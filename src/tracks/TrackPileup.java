package tracks;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.ListIterator;
import java.util.Map;
import java.util.TreeMap;

import org.apache.commons.lang3.StringUtils;

import coloring.Config;
import coloring.ConfigKey;
import exceptions.InvalidColourException;
import exceptions.InvalidGenomicCoordsException;
import exceptions.InvalidRecordException;
import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import samTextViewer.GenomicCoords;
import samTextViewer.Utils;

/** Collect depth and coverage info over a region of a chromosome
 * */
public class TrackPileup extends TrackWiggles {

	/** Key: Position of the locus. Value: information at this position. 
	 * NB: This is a hashmap so positions are not returned in order. The method(s)
	 * getDepth(), etc should return a TreeMap for positions to be sorted. Here we use
	 * HashMap because is faster to build. 
	 * */
	private Map<Integer, Locus> loci;
	private List<ScreenWiggleLocusInfo> screenWiggleLocusInfoList= new ArrayList<ScreenWiggleLocusInfo>();
	private long nRecsInWindow;
	private long alnRecCnt;
	
	/*        C O N S T R U C T O R         */

	/** Initialize pileup with a list of empty loci at the region given by chrom:from-to
	 * At the start, the loci map is empty as no information is provided yet. 
	 * Positions without info (i.e. 0 coverage) are not present in the map at all 
	 * (to save memory and cpu). 
	 * @throws IOException 
	 * @throws SQLException 
	 * @throws InvalidRecordException 
	 * @throws InvalidGenomicCoordsException 
	 * @throws ClassNotFoundException 
	 * */
	protected TrackPileup(String bam, GenomicCoords gc) throws IOException, ClassNotFoundException, InvalidGenomicCoordsException, InvalidRecordException, SQLException {

		if(!Utils.bamHasIndex(bam)){
			File temp= File.createTempFile("asciigenome.", ".bam");
			Utils.sortAndIndexSamOrBam(bam, temp.getAbsolutePath(), true);
			this.setWorkFilename(temp.getAbsolutePath());
		} else {
			this.setWorkFilename(bam);
		}
		this.setFilename(bam);
		this.setGc(gc);
		this.alnRecCnt= Utils.getAlignedReadCount(this.getWorkFilename());
		this.setLastModified();
		
	}

	/*       M E T H O D S        */
	
	@Override
	public void update() throws InvalidGenomicCoordsException, IOException{
		
		// this.updateWorkFileName();
		
		if(this.getyMaxLines() == 0){
			return;
		}
		
		SamReader samReader= Utils.getSamReader(this.getWorkFilename());
		List<Boolean> passFilter= this.filterReads(samReader);
		this.nRecsInWindow= 0;
		for(boolean x : passFilter){ // The count of reads in window is the count of reads passing filters
			if(x){
				this.nRecsInWindow++;
			}
		}
		samReader= Utils.getSamReader(this.getWorkFilename());
		Iterator<SAMRecord> sam= samReader.query(this.getGc().getChrom(), this.getGc().getFrom(), this.getGc().getTo(), false);
		
		this.loci= new TreeMap<Integer, Locus>();
		ListIterator<Boolean> pass = passFilter.listIterator();
		while(sam.hasNext()){
			SAMRecord rec= sam.next();
			if( pass.next() ){
					this.add(rec);
			}
		}
		List<Double> screenScores= this.prepareScreenScores();
		this.setScreenScores(screenScores);
	}

//	private void updateWorkFileName() {
//		
//		long modified;
//		UrlValidator urlValidator = new UrlValidator();
//		if(urlValidator.isValid(this.getFilename())){
//			modified= this.getLastModified();
//		} else {
//			modified= new File(this.getFilename()).lastModified();
//		}
//		if(this.getLastModified() < modified){
//
//		}
//	}

	private List<Double> prepareScreenScores() throws InvalidGenomicCoordsException, IOException{
		// We need to walk along the genomic window spanned by the current coordinates and 
		// collect depth. Depth as to be binned into screen scores.
		int userWindowSize= this.getGc().getUserWindowSize();
		
		this.screenWiggleLocusInfoList.clear();
		for(int i= 0; i < userWindowSize; i++){
			this.screenWiggleLocusInfoList.add(new ScreenWiggleLocusInfo());
		}		
		
		List<Double> mapping= this.getGc().getMapping();
		Map<Integer, Integer> depthMap = this.getDepth();
		for( int refPos : depthMap.keySet()){
			int screenIdx= Utils.getIndexOfclosestValue(refPos, mapping);
			ScreenWiggleLocusInfo sloc = this.screenWiggleLocusInfoList.get(screenIdx);
			int depth= depthMap.get(refPos);
			sloc.increment(depth);
		}

		List<Double> screenScores= new ArrayList<Double>();
		for(ScreenWiggleLocusInfo screenLocusInfo: this.screenWiggleLocusInfoList){
			double score= screenLocusInfo.getMeanScore();
			screenScores.add(score);
		}
		return screenScores;
	}

	@Override
	protected List<Double> getScreenScores(){
		return this.screenScores;
	}
	
	/** Update the map of loci with the information in this record. 
	 * */
	private void add(SAMRecord samRecord){
		
		// Is this read forward or reverse? First or second in pair?
		boolean isFirstOFPair= ! samRecord.getFirstOfPairFlag();
		boolean isReverse= samRecord.getReadNegativeStrandFlag();
		
		List<AlignmentBlock> alnBlocks = samRecord.getAlignmentBlocks();
		if(alnBlocks.size() == 0){
			// Nothing to be done. This may happen with e.g. fully clipped reads
			return; 
		}
		int readPos= alnBlocks.get(0).getReadStart() - 1; // Leftmost position where the read starts aligning.
		for(AlignmentBlock block : alnBlocks){
			
			if(block.getReferenceStart() > this.getGc().getTo()){
				// This block is completely to the right of the user's coordinates: Not useful
				continue; 
			}
			if((block.getReferenceStart() + block.getLength() - 1) < this.getGc().getFrom()){
				// This block is completely to the left of the user's coordinates.
				// So no useful information here.
				continue;
			}
			
			int refPos= block.getReferenceStart() - 1; // -1 because we increment immediately in the loop
			for(int i= 0; i < block.getLength(); i++){
				// Iterate through each base of this block
				refPos++;
				readPos++;
				if(refPos < this.getGc().getFrom()){
					continue; // Skip this position. Part of this block is outside user's coordinates. 
				}
				if(refPos > this.getGc().getTo()){
					break; // We have passed the user's endpoint, no need to process this record anymore
				}
				// What read base do we have at this position?
				char base= samRecord.getReadBases().length == 0 ? 'N' : (char) samRecord.getReadBases()[readPos-1];				
								
				// Start collecting info	
				if( ! this.loci.containsKey(refPos)){
					// Add this position to the map.
					this.loci.put(refPos, new Locus(this.getGc().getChrom(), refPos));
				}
				this.loci.get(refPos).add(base, isReverse, isFirstOFPair);
			}
		}
		// Now we need to increment counts corresponding to deletions in the reference
		List<int[]>deletedBlocks= this.getRefPositionOfDeletedBlocks(samRecord);
		for(int[] block : deletedBlocks){
			for(int refPos= block[0]; refPos <= block[1]; refPos++){
				if(refPos < this.getGc().getFrom() || refPos > this.getGc().getTo()){
					continue; // Position is outside user's coordinates.
				}
				if( ! this.loci.containsKey(refPos)){
					// Add this position to the map.
					this.loci.put(refPos, new Locus(this.getGc().getChrom(), refPos));
				}
				this.loci.get(refPos).add('D', isReverse, isFirstOFPair);
			}
		}
	}
	
	/**
	 * MEMO: Deletion does not consume read bases. It consumes reference bases:
	 * ref  NNNNNNNN
	 * read NNN---NN
	 * */
	private List<int[]> getRefPositionOfDeletedBlocks(SAMRecord samRecord){
		
		int ndel= StringUtils.countMatches(samRecord.getCigarString().toUpperCase(), 'D');
		List<int[]> deletedBlocks= new ArrayList<int[]>(ndel); // int array contains start and end of deletion
		if( ndel == 0 ){
			return deletedBlocks; // No deletetions, nothing to be done.
		}
		
		// 1S 2M 3D 4M 5D 6M 7S
		// Walk along the cigar string. When an operator consumes read bases, advance the position tracker
		// When you hit a deletion record the start and end of the deleted block in ref coordinates.
		List<CigarElement> cigarOps = samRecord.getCigar().getCigarElements();
		int readTracker= 0;
		int[] deletedBlock= new int[2]; // Length 2 because it contains start and end of deleted block
		for(CigarElement op : cigarOps){
			if(op.getOperator().equals(CigarOperator.DELETION)){
				deletedBlock[0]= samRecord.getReferencePositionAtReadPosition(readTracker) + 1;
				deletedBlock[1]= samRecord.getReferencePositionAtReadPosition(readTracker + 1) - 1;
				deletedBlocks.add(deletedBlock);
				deletedBlock= new int[2];
				if(deletedBlocks.size() == ndel){
					// We have accumulated all the deleted blocks. No need to process this read any more.
					return deletedBlocks; 
				}
			}
			if(op.getOperator().consumesReadBases()){
				readTracker += op.getLength();
			}
		}
		return deletedBlocks;
	}

	/** Depth at each position. Key: reference position. Value: depth. 
	 * */
	protected Map<Integer, Integer> getDepth(){
		// Important: Use TreeMap to have position sorted. 
		Map<Integer, Integer> depth= new TreeMap<Integer, Integer>();
		for(int pos : this.loci.keySet()){
			Locus loc= this.loci.get(pos);
			depth.put(pos, loc.getDepth());
		}
		return depth;
	}
	
	private char[] getConsensusSequence() throws IOException {
		
		// We could get the refseq from genomicCoords but maybe safer to extract it again from scratch.
		byte[] refSeq= null;
		if(this.getGc().getFastaFile() != null){
			IndexedFastaSequenceFile faSeqFile = new IndexedFastaSequenceFile(new File(this.getGc().getFastaFile()));
			refSeq= faSeqFile.getSubsequenceAt(this.getGc().getChrom(), this.getGc().getFrom(), this.getGc().getTo()).getBases();
			faSeqFile.close();
		}
		char[] consensusSequence= new char[this.getGc().getTo() - this.getGc().getFrom() + 1];
		int i= 0;
		for(int pos= this.getGc().getFrom(); pos <= this.getGc().getTo(); pos++){
			char consensus= ' '; // Empty char assuming there is no coverage.
			if(this.loci.containsKey(pos)){ // If containsKey then the position has coverage.
				Locus loc= this.loci.get(pos);
				consensus= loc.getConsensus();
				if(refSeq != null){
					char ref= Character.toUpperCase((char) refSeq[loc.pos - this.getGc().getFrom()]);
					if(ref == Character.toUpperCase(consensus)){
						consensus= '=';
					}
				}
			}
			consensusSequence[i]= consensus;
			i++;
		}
		return consensusSequence;
	}

	public String getPrintableConsensusSequence() throws IOException, InvalidGenomicCoordsException, InvalidColourException{
		if( ! this.getGc().isSingleBaseResolution || this.isBisulf()){
			return "";
		}
		
		if(new String(this.getConsensusSequence()).trim().isEmpty()){
			return ""; // If there is no coverage at all
		}
		String faSeqStr= "";
		for(char base : this.getConsensusSequence()){
			
			if(this.isNoFormat()){
				faSeqStr += base;
			} else { 
				faSeqStr += "\033[48;5;" + Config.get256Color(ConfigKey.background) + ";38;5;";
				     if(base == 'A') { faSeqStr += Config.get256Color(ConfigKey.seq_a);} 
				else if(base == 'C') { faSeqStr += Config.get256Color(ConfigKey.seq_c);} 
				else if(base == 'G') { faSeqStr += Config.get256Color(ConfigKey.seq_g);} 
				else if(base == 'T') { faSeqStr += Config.get256Color(ConfigKey.seq_t);} 
				else { faSeqStr += Config.get256Color(ConfigKey.seq_other); }
				faSeqStr += "m" + base + "\033[0m\033[38;5;0;48;5;" + Config.get256Color(ConfigKey.background) + "m"; // Clear formatting and fg to black and bg to white;
			}
		}
		return faSeqStr + "\n";
	}
	
	@Override
	public String getTitle() throws InvalidColourException, InvalidGenomicCoordsException, IOException{
		
		if(this.isHideTitle()){
			return "";
		}
		
		Double[] range = Utils.range(this.getScreenScores());
		
		String rpmTag= "";
		if(this.isRpm()){
			range[0]= (range[0] / this.alnRecCnt) * 1000000.0;
			range[1]= (range[1] / this.alnRecCnt) * 1000000.0;
			rpmTag= "; rpm";
		}
		Double[] rounded= Utils.roundToSignificantDigits(range[0], range[1], 2);
		

		String ymin= this.getYLimitMin().isNaN() ? "auto" : this.getYLimitMin().toString();
		String ymax= this.getYLimitMax().isNaN() ? "auto" : this.getYLimitMax().toString();

		String samtools= "";
		if( ! (this.get_F_flag() == Track.F_FLAG) ){
			samtools += " -F " + this.get_F_flag();
		}
		if( ! (this.get_f_flag() == Track.f_FLAG) ){
			samtools += " -f " + this.get_f_flag();
		}
		if( ! (this.getMapq() == Track.MAPQ) ){
			samtools += " -q " + this.getMapq();
		}
		if( ! samtools.isEmpty()){
			samtools= "; samtools" + samtools;
		}
		String awk= "";
		if(!this.getAwk().isEmpty()){
			awk= "; awk:on";
		}
		
		String xtitle= this.getTrackTag() 
				+ "; ylim[" + ymin + " " + ymax + "]" 
				+ "; range[" + rounded[0] + " " + rounded[1] + "]"
				+ "; Reads: " + this.nRecsInWindow + "/" + this.alnRecCnt
				+ samtools 
				+ rpmTag
				+ awk;
		// xtitle= Utils.padEndMultiLine(xtitle, this.getGc().getUserWindowSize());
		return this.formatTitle(xtitle) + "\n";
	}
 
	@Override
	public String getAwk(){
		return this.awk;
	}
	
	@Override
	public void setAwk(String awk) throws ClassNotFoundException, IOException, InvalidGenomicCoordsException, InvalidRecordException, SQLException {
		this.awk= awk;
		this.update();
	}

	@Override
	public void setRpm(boolean rpm){
		this.rpm= rpm;
	}
}
