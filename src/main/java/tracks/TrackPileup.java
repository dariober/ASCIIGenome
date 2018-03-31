package tracks;

import java.io.File;
import java.io.IOException;
import java.net.MalformedURLException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.ListIterator;
import java.util.Map;
import java.util.Random;
import java.util.regex.Pattern;

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
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.IntervalTree;
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
	private Map<String, Map<Integer, Locus>> loci= new HashMap<String, Map<Integer, Locus>>();
	@SuppressWarnings("rawtypes")
	private Map<String, IntervalTree> zeroDepthIntervals= new HashMap<String, IntervalTree>(); 
	
	private List<ScreenWiggleLocusInfo> screenWiggleLocusInfoList= new ArrayList<ScreenWiggleLocusInfo>();
	private long alnRecCnt= -1;
	
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

		this.setTrackFormat(TrackFormat.BAM);
		
		if(!Utils.bamHasIndex(bam)){
			File temp= Utils.createTempFile(".asciigenome.", ".bam");
			Utils.sortAndIndexSamOrBam(bam, temp.getAbsolutePath(), true);
			this.setWorkFilename(temp.getAbsolutePath());
		} else {
			this.setWorkFilename(bam);
		}
		this.setFilename(bam);
		this.setGc(gc);
		// this.alnRecCnt= Utils.getAlignedReadCount(this.getWorkFilename());
		this.setLastModified();
	}

	/*                  F I L T E R S           */
	@Override
	void setSamRecordFilter(List<SamRecordFilter> samRecordFilter) throws MalformedURLException, ClassNotFoundException, IOException, InvalidGenomicCoordsException, InvalidRecordException, SQLException {
		this.clearCache();
		this.getFeatureFilter().setSamRecordFilter(samRecordFilter);
		this.update();
	}
	
	@Override
	public void setShowHideRegex(Pattern showRegex, Pattern hideRegex) throws InvalidGenomicCoordsException, IOException, ClassNotFoundException, InvalidRecordException, SQLException{
		this.clearCache();
		this.getFeatureFilter().setShowHideRegex(showRegex, hideRegex);
		this.update();
	}
		
	@Override
	public void setAwk(String awk) throws ClassNotFoundException, IOException, InvalidGenomicCoordsException, InvalidRecordException, SQLException {
		this.clearCache();
		this.getFeatureFilter().setAwk(awk);
		this.update();
	}

	@Override 
	public String getAwk(){
		// MEMO: You need to override TrackWiggles not Tracks!
		return this.getFeatureFilter().getAwk();
	};

    @Override
    public void setVariantReadInInterval(String chrom, int from, int to, boolean variantOnly) throws MalformedURLException, ClassNotFoundException, IOException, InvalidGenomicCoordsException, InvalidRecordException, SQLException{
    	this.clearCache();
        super.setVariantReadInInterval(chrom, from, to, variantOnly);
    }
	
	/*       M E T H O D S        */
	
	@SuppressWarnings({ "rawtypes", "unchecked" })
	@Override
	public void update() throws InvalidGenomicCoordsException, IOException{
		
		if(this.getyMaxLines() == 0){
			return;
		}
		String chrom= this.getGc().getChrom();
		
		if(! this.loci.containsKey(chrom)){
			this.loci.put(chrom, new HashMap<Integer, Locus>());
		}
		if(! this.zeroDepthIntervals.containsKey(chrom)){
			this.zeroDepthIntervals.put(chrom, new IntervalTree());
		}

		// Check cache is not growing too much
		if(this.loci.get(chrom).keySet().size() > 500000){
			this.loci.get(chrom).clear();
			this.zeroDepthIntervals.get(chrom).clear();
		}
		
		// Find the positions that we haven't visited before:
		List<Integer> missingPos= new ArrayList<Integer>();
		for(int pos= this.getGc().getFrom(); pos <= this.getGc().getTo(); pos++){
			if( ! this.loci.get(chrom).containsKey(pos) &&  
				! this.zeroDepthIntervals.get(chrom).overlappers(pos, pos).hasNext()){
				missingPos.add(pos);
			}
		}
		for(List<Integer> gap : mergePositionsInIntervals(missingPos)){

			int qryFrom= gap.get(0);
			int qryTo= gap.get(1);
			
			SamReader samReader= Utils.getSamReader(this.getWorkFilename());
			List<Boolean> passFilter= this.filterReads(samReader, chrom, qryFrom, qryTo);
			samReader= Utils.getSamReader(this.getWorkFilename());
			
			Iterator<SAMRecord> sam= samReader.query(chrom, qryFrom, qryTo, false);
	
			ListIterator<Boolean> pass = passFilter.listIterator();
			while(sam.hasNext()){
				SAMRecord rec= sam.next();
				if(pass.next()){
					this.add(rec, qryFrom, qryTo, this.loci.get(chrom));
				}
			}
			
			// Now add the loci that have been collected in this last update
			List<Integer> zeroDepthPos= new ArrayList<Integer>();
			for(int pos= qryFrom; pos <= qryTo; pos++){
				if( ! this.loci.get(chrom).containsKey(pos)){
					zeroDepthPos.add(pos);
				}
			}
			List<List<Integer>> zeroDepthGaps = mergePositionsInIntervals(zeroDepthPos);
			for(List<Integer> z : zeroDepthGaps){
				this.zeroDepthIntervals.get(chrom).put(z.get(0), z.get(1), null);
			}
		}
		List<Float> screenScores= this.prepareScreenScores();
		this.setScreenScores(screenScores);
	}

	private List<Float> prepareScreenScores() throws InvalidGenomicCoordsException, IOException{
		// We need to walk along the genomic window spanned by the current coordinates and 
		// collect depth. Depth as to be binned into screen scores.
		int userWindowSize= this.getGc().getUserWindowSize();
		
		this.screenWiggleLocusInfoList.clear();
		for(int i= 0; i < userWindowSize; i++){
			this.screenWiggleLocusInfoList.add(new ScreenWiggleLocusInfo());
		}		
		
		List<Double> mapping= this.getGc().getMapping();
		Map<Integer, Integer> depthMap = this.getDepth(this.getGc().getChrom(), this.getGc().getFrom(), this.getGc().getTo());
		// Winsorise here:
//		if(this.getWinsorizeMultiple() > 0){
//			Map<Integer, Float> depthMapWins= depthMap; 
//			List<Float> xwin= Utils.winsorise(depthMap.entrySet(), this.getWinsorizeMultiple());
//			for x in depthMapWins:
//				//
//		}
		
		for(int refPos : depthMap.keySet()){
			int screenIdx= Utils.getIndexOfclosestValue(refPos, mapping);
			ScreenWiggleLocusInfo sloc = this.screenWiggleLocusInfoList.get(screenIdx);
			int depth= depthMap.get(refPos);
			sloc.increment(depth);
		}

		List<Float> screenScores= new ArrayList<Float>();
		for(ScreenWiggleLocusInfo screenLocusInfo: this.screenWiggleLocusInfoList){
			float score= screenLocusInfo.getMeanScore();
			screenScores.add(score);
		}
		return screenScores;
	}

	@Override
	protected List<Float> getScreenScores(){
		return this.screenScores;
	}
	
	/** Update the given map of loci with the information in this record. Only consider positions
	 * between qryFrom and qryTo. 
	 * */
	private void add(SAMRecord samRecord, int qryFrom, int qryTo, Map<Integer, Locus> accumulator){
		
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
			
			if(block.getReferenceStart() > qryTo){
				// This block is completely to the right of the user's coordinates: Not useful
				continue; 
			}
			if((block.getReferenceStart() + block.getLength() - 1) < qryFrom){
				// This block is completely to the left of the user's coordinates.
				// So no useful information here.
				continue;
			}
			
			int refPos= block.getReferenceStart() - 1; // -1 because we increment immediately in the loop
			for(int i= 0; i < block.getLength(); i++){
				// Iterate through each base of this block
				refPos++;
				readPos++;
				if(refPos < qryFrom){
					continue; // Skip this position. Part of this block is outside user's coordinates. 
				}
				if(refPos > qryTo){
					break; // We have passed the user's endpoint, no need to process this record anymore
				}
				// What read base do we have at this position?
				char base= samRecord.getReadBases().length == 0 ? 'N' : (char) samRecord.getReadBases()[readPos-1];				

				// Start collecting info	
				if( ! accumulator.containsKey(refPos)){
					// Add this position to the map.
					accumulator.put(refPos, new Locus(this.getGc().getChrom(), refPos));
				}
				accumulator.get(refPos).add(base, isReverse, isFirstOFPair);
			}
		}
		// Now we need to increment counts corresponding to deletions in the reference
		List<int[]>deletedBlocks= this.getRefPositionOfDeletedBlocks(samRecord);
		for(int[] block : deletedBlocks){
			for(int refPos= block[0]; refPos <= block[1]; refPos++){
				if(refPos < qryFrom || refPos > qryTo){
					continue; // Position is outside user's coordinates.
				}
				if( ! accumulator.containsKey(refPos)){
					// Add this position to the map.
					accumulator.put(refPos, new Locus(this.getGc().getChrom(), refPos));
				}
				accumulator.get(refPos).add('D', isReverse, isFirstOFPair);
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
	 * @throws IOException 
	 * */
	protected Map<Integer, Integer> getDepth(String chrom, int from, int to) throws IOException{
		
		List<Integer> posInWindow= new ArrayList<Integer>();  
		for(int pos : this.loci.get(chrom).keySet()){
			if(pos < from || pos > to){
				// Because of locus caching, some loci maybe outside the current window.
				continue;
			} else {
				posInWindow.add(pos);
			}
		}

		double samplingRate= (200000.0) / posInWindow.size();
		Random rand = new Random();
		// Important: Use have positions returned sorted. 
		Map<Integer, Integer> depth= new LinkedHashMap<Integer, Integer>();
		for(int pos : posInWindow){
			int posDepth= this.loci.get(chrom).get(pos).getDepth();
			if(rand.nextFloat() < samplingRate){
				depth.put(pos, posDepth);			
			}
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
			if(this.loci.get(this.getGc().getChrom()).containsKey(pos)){ // If containsKey then the position has coverage.
				Locus loc= this.loci.get(this.getGc().getChrom()).get(pos);
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
	protected String getTitleForActiveFilters() {
		List<String> title= new ArrayList<String>();
		if( ! this.getAwk().equals(Filter.DEFAULT_AWK.getValue())){
			title.add("awk");
		}
		if( ! this.getShowRegex().pattern().equals(Filter.DEFAULT_SHOW_REGEX.getValue()) || ! this.getHideRegex().pattern().equals(Filter.DEFAULT_HIDE_REGEX.getValue())){
			title.add("grep");
		}
		if( this.get_f_flag() != Integer.valueOf(Filter.DEFAULT_f_FLAG.getValue()) || 
			this.get_F_flag() != Integer.valueOf(Filter.DEFAULT_F_FLAG.getValue())){
			title.add("bit-flag");
		}
		if(this.getMapq() != Integer.valueOf(Filter.DEFAULT_MAPQ.getValue())){
			title.add("mapq");
		}
		if( ! this.getFeatureFilter().getVariantChrom().equals(Filter.DEFAULT_VARIANT_CHROM.getValue())){
			title.add("var-read");
		}
		if(title.size() > 0){
			return "; filters: " + title.toString(); 
		} else {
			return "";	
		}
	}
	
	@Override
	public String getTitle() throws InvalidColourException, InvalidGenomicCoordsException, IOException{
		
		if(this.isHideTitle()){
			return "";
		}
		
		Float[] range = Utils.range(this.getScreenScores());
		
		String rpmTag= "";
		if(this.isRpm()){
			if(this.alnRecCnt == -1){
				this.alnRecCnt= Utils.getAlignedReadCount(this.getWorkFilename());
			}
			range[0]= (float) ((range[0] / this.alnRecCnt) * 1000000.0);
			range[1]= (float) ((range[1] / this.alnRecCnt) * 1000000.0);
			rpmTag= "; rpm";
		}
		String[] rounded= Utils.roundToSignificantDigits(range[0], range[1], 2);
		

		String ymin= this.getYLimitMin().isNaN() ? "auto" : this.getYLimitMin().toString();
		String ymax= this.getYLimitMax().isNaN() ? "auto" : this.getYLimitMax().toString();

		String libsize= "";
		if(this.alnRecCnt != -1){
			libsize= "; lib size: " + this.alnRecCnt;
		}
		String xtitle= this.getTrackTag() 
				+ "; ylim[" + ymin + " " + ymax + "]" 
				+ "; range[" + rounded[0] + " " + rounded[1] + "]"
				+ libsize
				+ rpmTag
				+ this.getTitleForActiveFilters();
		return this.formatTitle(xtitle) + "\n";
	}

	@Override
	public void setRpm(boolean rpm){
		this.rpm= rpm;
	}

	/**Merge the *sorted* list of positions into intervals of consecutive ints.
	 * */
	protected static List<List<Integer>> mergePositionsInIntervals(List<Integer> positions){
		List<Integer> v= new ArrayList<Integer>();
		v.add(null); v.add(null); 
        int from= -1;
        int to= -1;
        List<List<Integer>> intervals= new ArrayList<List<Integer>>();
		for(int i : positions){
		    if(v.get(0) == null){
		        from= i;
		        to= i;
		        v.set(0, from); v.set(1, to);
		    } else if(i == to){
		    	continue; // You have duplicate positions.
		    } else if(i < to){
		    	System.err.println("Positions are not sorted: " + i + " after " + to);
		    	throw new RuntimeException();
		    } else if (i == (to + 1)){
		    	to= i;
		    	v.set(0, from);
		    	v.set(1, i);
		    } else {
		    	intervals.add(v);
		        from= i;
		        to= i;
		        v= new ArrayList<Integer>();
		        v.add(from); v.add(to);
			}
		}
		if(v.get(0) != null){
		    intervals.add(v);
		} else {
			intervals.clear();
		}
		return intervals;
	}

	private void clearCache(){
		this.loci.clear(); // clear cached positions
		this.zeroDepthIntervals.clear();
	}
	
	@Override
	public void reload() throws InvalidGenomicCoordsException, IOException, ClassNotFoundException, InvalidRecordException, SQLException{
		if( ! Files.isSameFile(Paths.get(this.getWorkFilename()), Paths.get(this.getFilename()))){
			TrackPileup tr= new TrackPileup(this.getFilename(), this.getGc());
			String fname= this.getWorkFilename();
			Files.move(Paths.get(tr.getWorkFilename()), Paths.get(fname), java.nio.file.StandardCopyOption.REPLACE_EXISTING);
			Files.move(Paths.get(tr.getWorkFilename().replaceAll("\\.bam$", ".bai")), Paths.get(fname.replaceAll("\\.bam$", ".bai")), java.nio.file.StandardCopyOption.REPLACE_EXISTING);
		}
		this.clearCache();
		this.update();
	}
}
