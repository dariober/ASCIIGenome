package tracks;

import java.io.File;
import java.io.IOException;
import java.net.MalformedURLException;
import java.net.URL;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import org.apache.commons.lang3.StringUtils;
import org.apache.commons.validator.routines.UrlValidator;

import com.google.common.base.Joiner;

import coloring.Config;
import coloring.ConfigKey;
import coloring.Xterm256;
import exceptions.InvalidColourException;
import exceptions.InvalidGenomicCoordsException;
import exceptions.InvalidRecordException;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import samTextViewer.GenomicCoords;
import samTextViewer.SamLocusIterator;
import samTextViewer.Utils;

@Deprecated
public class TrackCoverage extends Track {

	/* A t t r i b u t e s */
	
	private List<ScreenLocusInfo> screenLocusInfoList= new ArrayList<ScreenLocusInfo>(); 
	/** Library size as extracted from querying bam index */
	private long alnRecCnt= -1;
	/** N. records in window */
	private long nRecsInWindow= -1;
	private SamReader samReader;

	/* C o n s t r u c t o r */

	/**
	 * Construct coverage track from bam alignment in the provided interval. 
	 * Loci will be sampled according to the size of the interval and the size of the printable screen. 
	 * @param bam Input bam file
	 * @param gc Interval to sample positions from
	 * @param windowSize The size of the screen in number of characters.
	 * @param bs Should loci be parsed also as BS-Seq data? 
	 * @throws IOException 
	 * @throws InvalidGenomicCoordsException 
	 * @throws SQLException 
	 * @throws InvalidRecordException 
	 * @throws ClassNotFoundException 
	 */
	public TrackCoverage(String bam, GenomicCoords gc, boolean bs) throws IOException, InvalidGenomicCoordsException, ClassNotFoundException, InvalidRecordException, SQLException{
		this.setFilename(bam);
		this.setWorkFilename(bam);
		this.setBisulf(bs);
		this.alnRecCnt= Utils.getAlignedReadCount(bam);
		/*  ------------------------------------------------------ */
		/* This chunk prepares SamReader from local bam or URL bam */
		UrlValidator urlValidator = new UrlValidator();
		SamReaderFactory srf=SamReaderFactory.make();
		srf.validationStringency(ValidationStringency.SILENT);
		if(urlValidator.isValid(this.getWorkFilename())){
			this.samReader = srf.open(SamInputResource.of(new URL(bam)).index(new URL(bam + ".bai")));
		} else {
			this.samReader= srf.open(new File(bam));
		}
		/*  ------------------------------------------------------ */
		this.setGc(gc);
	}
	
	/* M e t h o d s */
	@Override
	protected void update() throws IOException, InvalidGenomicCoordsException{
		
		this.screenLocusInfoList= new ArrayList<ScreenLocusInfo>();
		if(this.getGc().getGenomicWindowSize() < this.MAX_REGION_SIZE * 5){
			
			IntervalList il= new IntervalList(samReader.getFileHeader());
			il.add(new Interval(this.getGc().getChrom(), this.getGc().getFrom(), this.getGc().getTo()));
			SamLocusIterator samLocIter= new SamLocusIterator(samReader, il, true);
			samLocIter.setSamFilters(this.getSamRecordFilter());
			Iterator<samTextViewer.SamLocusIterator.LocusInfo> iter= samLocIter.iterator();
			
			int userWindowSize= this.getGc().getUserWindowSize();
			
			for(int i= 0; i < this.getGc().getMapping(userWindowSize).size(); i++){
				this.screenLocusInfoList.add(new ScreenLocusInfo());	
			}
			
			List<Double> mapping = this.getGc().getMapping(userWindowSize);

			while(iter.hasNext()){

				samTextViewer.SamLocusIterator.LocusInfo locusInfo= iter.next();
				int screenPos= Utils.getIndexOfclosestValue(locusInfo.getPosition(), mapping);
				byte refBase= '\0';
				if(this.getGc().getRefSeq() != null){
					refBase= this.getGc().getRefSeq()[screenPos];
				}
				this.screenLocusInfoList.get(screenPos).increment(locusInfo, refBase);
			}


			this.nRecsInWindow= Utils.countReadsInWindow(this.getWorkFilename(), this.getGc(), this.getSamRecordFilter());
			samLocIter.close();
		}

		ArrayList<Double> yValues = new ArrayList<Double>();
		for(ScreenLocusInfo x : this.screenLocusInfoList){
			yValues.add(x.getMeanDepth());
		}
		this.setScreenScores(yValues);
			
		if(this.isRpm()){
			this.rpm();
		}		
	}

	/** Convert screen scores to RPM. Raw values overwritten.
	 * */
	private void rpm(){
		ArrayList<Double> rpm = new ArrayList<Double>();	
		// long libSize= Utils.getAlignedReadCount(new File(this.getFilename()));
		for(Double y : this.getScreenScores()){
			rpm.add(y / this.alnRecCnt * 1000000.0);
		}		
		this.setScreenScores(rpm);
	} 
	
	@Override
	protected void updateToRPM(){
		try {
			this.update(); // It shouldn't be necessary to go for full update!!!
		} catch (IOException | InvalidGenomicCoordsException e) {
			e.printStackTrace();
		}
	}
	
	
	protected void setSamRecordFilter(List<SamRecordFilter> samRecordFilter) throws MalformedURLException, ClassNotFoundException, IOException, InvalidGenomicCoordsException, InvalidRecordException, SQLException {
		this.samRecordFilter = samRecordFilter;
		this.update();
	}

	
	/**
	 * Printable coverage track. The height of the track in lines is `yMaxLines`.
	 * @param screenToGenomeMap List of genomic positions corresponding to each column on screen.
	 * @param yMaxLines
	 * @param rpm Should read counts be normalized by library size as Reads Per Million
	 * @return HashMapwith with keys/values the printable characteristics of the track. 
	 * @throws IOException 
	 * @throws InvalidGenomicCoordsException 
	 * @throws InvalidColourException 
	 */
	@Override
	public String printToScreen() throws InvalidGenomicCoordsException, IOException, InvalidColourException{
		 //This method should not do any computation like RPM etc. Just print stuff.
		final Xterm256 xterm256= new Xterm256();
		if(this.getyMaxLines() == 0){
			return "";
		} else if(this.screenLocusInfoList.size() == 0){
			if(this.getGc().getGenomicWindowSize() >= this.MAX_REGION_SIZE * 5){
				return "Track not shown: Window is too large";
			}
			return "";
		}
		
		TextProfile textProfile= new TextProfile(this.getScreenScores(), this.getyMaxLines(), this.getYLimitMin(), this.getYLimitMax());
		ArrayList<String> lineStrings= new ArrayList<String>();
		for(int i= (textProfile.getProfile().size() - 1); i >= 0; i--){
			List<String> xl= textProfile.getProfile().get(i);
			lineStrings.add(StringUtils.join(xl, ""));
		}
		String printable= Joiner.on("\n").join(lineStrings); 
		if(!this.isNoFormat()){
			printable= "\033[48;5;"
			+ Config.get256Color(ConfigKey.background)
			+ ";38;5;"
			+ xterm256.colorNameToXterm256(this.getTitleColour())
			+ "m"
			+ printable;
		}
		return printable;
		
	}
	
    /* S e t t e r s   and   G e t t e r s */
    
    /* This method makes sense calling only after having set the profile. Typically after */
    public List<ScreenLocusInfo> getScreenLocusInfoList(){
    	return screenLocusInfoList;
    }

	@Override
	public String getTitle() throws InvalidColourException, InvalidGenomicCoordsException, IOException{
		
		if(this.isHideTitle()){
			return "";
		}
		
		Double[] range = Utils.range(this.getScreenScores());
		Double[] rounded= Utils.roundToSignificantDigits(range[0], range[1], 2);
		String rpmTag= this.isRpm() ? "; rpm" : "";

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
		String xtitle= this.getTrackTag() 
				+ "; ylim[" + ymin + " " + ymax + "]" 
				+ "; range[" + rounded[0] + " " + rounded[1] + "]"
				+ "; Recs here/all: " + this.nRecsInWindow + "/" + this.alnRecCnt
				+ samtools 
				+ rpmTag;
		// xtitle= Utils.padEndMultiLine(xtitle, this.getGc().getUserWindowSize());
		return this.formatTitle(xtitle) + "\n";
	}

	
	private char[] getConsensusSequence() throws IOException {
		
		IntervalList il= new IntervalList(this.samReader.getFileHeader());
		il.add(new Interval(this.getGc().getChrom(), this.getGc().getFrom(), this.getGc().getTo()));
		SamLocusIterator samLocIter= new SamLocusIterator(this.samReader, il, true);
		samLocIter.setSamFilters(this.getSamRecordFilter());
		Iterator<samTextViewer.SamLocusIterator.LocusInfo> iter= samLocIter.iterator();
	
		// We could get the refseq from genomicCoords but maybe safer to extract it again from scratch.
		byte[] refSeq= null;
		if(this.getGc().getFastaFile() != null){
			IndexedFastaSequenceFile faSeqFile = new IndexedFastaSequenceFile(new File(this.getGc().getFastaFile()));
			refSeq= faSeqFile.getSubsequenceAt(this.getGc().getChrom(), this.getGc().getFrom(), this.getGc().getTo()).getBases();
			faSeqFile.close();
		}
		
		char[] consensusSequence= new char[this.getGc().getTo() - this.getGc().getFrom() + 1];
		int i= 0;
		while(iter.hasNext()){
			samTextViewer.SamLocusIterator.LocusInfo locusInfo= iter.next();
			char ref= '.';
			if(refSeq != null){
				ref= Character.toUpperCase((char) refSeq[locusInfo.getPosition() - this.getGc().getFrom()]);
			}
			if( ! this.isBisulf()){
				consensusSequence[i]= (new PileupLocus_DEPRECATE_ME(locusInfo, ref)).getConsensus();
			} else {
				// consensusSequence[i]= (new MethylLocus(locusInfo, ref)).getConsensus();
			}
			i++;
		}
		samLocIter.close();
		return consensusSequence;
	}

	@Override
	public String getPrintableConsensusSequence() throws IOException, InvalidGenomicCoordsException, InvalidColourException{
		if( ! this.getGc().isSingleBaseResolution || this.isBisulf()){
			return "";
		}
		String faSeqStr= "";
		boolean allEmpty= true;
		for(char base : this.getConsensusSequence()){
			
			if(base != ' '){ // Switch to check if entire sequence doesn't have any coverage. 
				allEmpty= false;
			}
			
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
		if(allEmpty){
			return "";
		} 
		return faSeqStr + "\n";
	}
	
	/** This method is not really used. For each position print ACGT counts. 
	 * */
	protected List<PileupLocus_DEPRECATE_ME> getPileupList() throws IOException {
		
		IntervalList il= new IntervalList(this.samReader.getFileHeader());
		il.add(new Interval(this.getGc().getChrom(), this.getGc().getFrom(), this.getGc().getTo()));
		SamLocusIterator samLocIter= new SamLocusIterator(this.samReader, il, true);
		samLocIter.setSamFilters(this.getSamRecordFilter());
		Iterator<samTextViewer.SamLocusIterator.LocusInfo> iter= samLocIter.iterator();
	
		// We could get the refseq from genomicCoords but maybe safer to extract it again from scratch.
		byte[] refSeq= null;
		if(this.getGc().getFastaFile() != null){
			IndexedFastaSequenceFile faSeqFile = new IndexedFastaSequenceFile(new File(this.getGc().getFastaFile()));
			refSeq= faSeqFile.getSubsequenceAt(this.getGc().getChrom(), this.getGc().getFrom(), this.getGc().getTo()).getBases();
			faSeqFile.close();
		}
		
		List<PileupLocus_DEPRECATE_ME> pileup= new ArrayList<PileupLocus_DEPRECATE_ME>(); 
		while(iter.hasNext()){			
			samTextViewer.SamLocusIterator.LocusInfo locusInfo= iter.next();
			char ref= '.';
			if(refSeq != null){
				ref= Character.toUpperCase((char) refSeq[locusInfo.getPosition() - this.getGc().getFrom()]);
			}
			pileup.add(new PileupLocus_DEPRECATE_ME(locusInfo, ref));
		}
		samLocIter.close();
		return pileup;
	}
	
	@Override
	public List<String> printPileupList(){
		List<String> plist= new ArrayList<String>();
		try {
			for(PileupLocus_DEPRECATE_ME x : this.getPileupList()){
				plist.add(x.toString() + "\n");
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
		return plist;
	}

	@Override
	public String printLines(){
		return "";
	}
//	public String printFeaturesToFile() throws IOException, InvalidGenomicCoordsException, InvalidColourException {
//		return "";
//	}

	@Override
	protected List<String> getRecordsAsStrings() {
		return new ArrayList<String>();
	}

	@Override
	public void setAwk(String awk) throws ClassNotFoundException, IOException, InvalidGenomicCoordsException,
			InvalidRecordException, SQLException {
		// 
	}

	@Override
	public String getAwk() {
		// TODO Auto-generated method stub
		return null;
	}
}
