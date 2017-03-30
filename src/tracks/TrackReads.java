package tracks;

import java.io.File;
import java.io.IOException;
import java.net.URISyntaxException;
import java.net.URL;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Random;

import org.apache.commons.io.FileUtils;
import org.apache.commons.lang3.StringUtils;
import org.apache.commons.validator.routines.UrlValidator;

import exceptions.InvalidColourException;
import exceptions.InvalidGenomicCoordsException;
import exceptions.InvalidRecordException;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.filter.AggregateFilter;
import samTextViewer.GenomicCoords;
import samTextViewer.Main;
import samTextViewer.Utils;

/**
 * Prepare track for printing aligned reads 
 * @author berald01
 *
 */
public class TrackReads extends Track{

	private List<List<TextRead>> readStack;
	// private boolean bisulf= false;
	private boolean withReadName= false;
	private static int MAX_READS_STACK= 2000;
	private long nRecsInWindow= -1;
	/* C o n s t r u c t o r s */
	/**
	 * Create read track
	 * @param sam Bam file to extract reads from
	 * @param gc Interval to consider
	 * @bs Should the track be displayed as BS-Seq data
	 * @param maxReadsStack Accumulate at most this many reads.
	 * @throws IOException 
	 * @throws InvalidGenomicCoordsException 
	 * @throws SQLException 
	 * @throws InvalidRecordException 
	 * @throws ClassNotFoundException 
	 */
	public TrackReads(String bam, GenomicCoords gc) throws IOException, InvalidGenomicCoordsException, ClassNotFoundException, InvalidRecordException, SQLException{

		if(!Utils.bamHasIndex(bam)){
			System.err.println("\nAlignment file " + bam + " has no index.\n");
			throw new RuntimeException();
		}
		this.setFilename(bam);
		this.setWorkFilename(bam);
		this.setGc(gc);

	}
	
	protected TrackReads(){};
	
	/* M e t h o d s */
	
	public void update() throws InvalidGenomicCoordsException, IOException{
		
		this.readStack= new ArrayList<List<TextRead>>();
		if(this.getGc().getGenomicWindowSize() < this.MAX_REGION_SIZE){

			/*  ------------------------------------------------------ */
			/* This chunk prepares SamReader from local bam or URL bam */
			UrlValidator urlValidator = new UrlValidator();
			SamReaderFactory srf=SamReaderFactory.make();
			srf.validationStringency(ValidationStringency.SILENT);
			SamReader samReader;
			if(urlValidator.isValid(this.getWorkFilename())){
				samReader = srf.open(SamInputResource.of(new URL(this.getWorkFilename())).index(new URL(this.getWorkFilename() + ".bai")));
			} else {
				samReader= srf.open(new File(this.getWorkFilename()));
			}
			/*  ------------------------------------------------------ */
			
			this.nRecsInWindow= Utils.countReadsInWindow(this.getWorkFilename(), this.getGc(), this.getSamRecordFilter());
			float probSample= (float) TrackReads.MAX_READS_STACK / this.nRecsInWindow;
			
			Iterator<SAMRecord> sam= samReader.query(this.getGc().getChrom(), this.getGc().getFrom(), this.getGc().getTo(), false);
			List<TextRead> textReads= new ArrayList<TextRead>();
			AggregateFilter aggregateFilter= new AggregateFilter(this.getSamRecordFilter());
			
			while(sam.hasNext() && textReads.size() < TrackReads.MAX_READS_STACK){
	
				SAMRecord rec= sam.next();

				Boolean passAwk= true; // Utils.passAwkFilter(rec.getSAMString(), this.getAwk());
				
				if( !rec.getReadUnmappedFlag() && !aggregateFilter.filterOut(rec) && passAwk){
					Random rand = new Random();
					if(rand.nextFloat() < probSample){ // Downsampler
						TextRead tr= new TextRead(rec, this.getGc());
						textReads.add(tr);
					}
				}
			}
			this.readStack= stackReads(textReads);
		} else {
			this.nRecsInWindow= -1;
		}
	}
	
	/** 
	 * Printable track on screen. This is what should be called by Main 
	 * @throws InvalidGenomicCoordsException 
	 * @throws InvalidColourException */
	@Override
	public String printToScreen() throws InvalidGenomicCoordsException, InvalidColourException{
		
		int yMaxLines= (this.getyMaxLines() < 0) ? Integer.MAX_VALUE : this.getyMaxLines();;
		
		// If there are more lines (inner lists) than desired lines of output (yMaxLines), get a representative sample
		List<Double> keep= new ArrayList<Double>();
		if(this.readStack.size() == 0){
			return  "";
		} else if(this.readStack.size() > yMaxLines){
			keep= Utils.seqFromToLenOut(0, this.readStack.size()-1, yMaxLines);
		} else {
			keep= Utils.seqFromToLenOut(0, this.readStack.size()-1, this.readStack.size());
		}
		StringBuilder printable= new StringBuilder();
		for(Double idx : keep){
			List<TextRead> line= this.readStack.get((int)Math.rint(idx));
			try {
				printable.append(linePrinter(line, this.bisulf, this.isNoFormat(), withReadName));
				printable.append("\n");
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
		return printable.toString().replaceAll("\n$", "");
	}
	
	/**		
	 * Put in the same list reads that will go in the same line of text 
	 * Example Input, a list of TextRead's:
	 AAAAAAAAAAAA
	  CCCCCCCCCCCC
	              TTTTTTTTTTT
	                           GGGGGGGGGGG
	                                   AAAAAAAA

	 * Output, each line is a list of TextRead:
     [AAAAAAAAAAAA TTTTTTTTTTT  GGGGGGGGGGG]       
	 [ CCCCCCCCCCCC                     AAAAAAAA]
	 * @throws IOException 
	 * @throws InvalidGenomicCoordsException 
	 */
	private List<List<TextRead>> stackReads(List<TextRead> textReads) throws InvalidGenomicCoordsException, IOException{
		
		List<List<TextRead>> listOfLines= new ArrayList<List<TextRead>>();
		if(textReads.size() == 0){
			return listOfLines;
		}
		List<TextRead> line= new ArrayList<TextRead>();
		line.add(textReads.get(0)); 
		textReads.remove(0);
		listOfLines.add(line);
		int gap= (this.getGc().getBpPerScreenColumn() > 1) ? 0 : 1; // If reads are very compressed, do not add space between adjacent ones.
		while(true){
			ArrayList<TextRead> trToRemove= new ArrayList<TextRead>();
			// Find a read in input whose start is greater then end of current
			for(int i=0; i < textReads.size(); i++){
				TextRead tr= textReads.get(i);
				if(tr.getTextStart() > line.get(line.size()-1).getTextEnd()+gap){ // +2 because we want some space between adjacent reads
					listOfLines.get(listOfLines.size()-1).add(tr); // Append to the last line. 
					trToRemove.add(tr);
				}
			} // At the end of the loop you have put in line as many reads as you can. 
			for(TextRead tr : trToRemove){ 
				textReads.remove(textReads.indexOf(tr));
			}
			// Create a new line, add the first textRead in list
			if(textReads.size() > 0){
				line= new ArrayList<TextRead>();
				line.add(textReads.get(0));
				listOfLines.add(line);
				textReads.remove(0);
			} else {
				break;
			}
		}
		return listOfLines;
	}
	
	/** Prepare a printable string of each output line. 
	 * @param textReads List reads to print out on the same line.
	 * @param noFormat Do not format reads.
	 * @return
	 * @throws IOException 
	 * @throws InvalidGenomicCoordsException 
	 * @throws InvalidColourException 
	 */
	private String linePrinter(List<TextRead> textReads, boolean bs, boolean noFormat, boolean withReadName) throws IOException, InvalidGenomicCoordsException, InvalidColourException{
		StringBuilder sb= new StringBuilder();

		int curPos= 0; // Position on the line, needed to pad with blanks btw reads.
		double bpPerScreenColumn= this.getGc().getBpPerScreenColumn();
		for(TextRead tr : textReads){
			String line= StringUtils.repeat(" ", (tr.getTextStart()-1) - curPos);
			sb.append(line);
			String printableRead= tr.getPrintableTextRead(bs, noFormat, withReadName, bpPerScreenColumn);
			sb.append(printableRead);
			curPos= tr.getTextEnd();
		}
		return sb.toString();
	}

	
	@Override
	public String getTitle() throws InvalidColourException, InvalidGenomicCoordsException, IOException{
		
		if(this.isHideTitle()){
			return "";
		}
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
		String title= this.getTrackTag()
				+ samtools;
		
		//String xtitle= Utils.padEndMultiLine(title, this.getGc().getUserWindowSize());
		return this.formatTitle(title) + "\n";
	}
	
	@Override
	protected List<String> getRecordsAsStrings() {
		List<String> featureList= new ArrayList<String>();
		
		for(List<TextRead> x : this.readStack){
			for(TextRead txr : x ){
				featureList.add(txr.getSamRecord().getSAMString());
			}
		}
		return featureList;
	}

	
	/**Print raw features under track. 
	 * windowSize size the number of characters before clipping occurs. This is 
	 * typically the window size for plotting. windowSize is used only by CLIP mode.  
	 * @throws IOException 
	 * @throws InvalidGenomicCoordsException 
	 * */
//	private String printFeatures() throws InvalidGenomicCoordsException, IOException{
//
//		int windowSize= this.getGc().getUserWindowSize();
//		if(this.getPrintMode().equals(PrintRawLine.FULL)){
//			windowSize= Integer.MAX_VALUE;
//		} else if(this.getPrintMode().equals(PrintRawLine.CLIP)){
//			// Keep windowSize as it is
//		} else {
//			return "";
//		} 
//		
//		List<String> featureList= new ArrayList<String>();
//		
//		int count= this.getPrintRawLineCount();
//		String omitString= "";
//		outerloop:
//		for(List<TextRead> x : this.readStack){
//			for(TextRead txr : x ){
//				featureList.add(txr.getSamRecord().getSAMString());
//				count--;
//				if(count == 0){
//					int omitted= (int) (this.nRecsInWindow - this.getPrintRawLineCount());
//					if(omitted > 0){
//						omitString= "[" + omitted + "/"  + this.nRecsInWindow + " features omitted]";
//					}
//					break outerloop;
//				}
//			}
//		}
//		List<String> tabList= Utils.tabulateList(featureList);
//		StringBuilder sb= new StringBuilder();
//		if( ! omitString.isEmpty()){
//			sb.append(omitString + "\n");
//		}
//		for(String x : tabList){
//			if(x.length() > windowSize){
//				x= x.substring(0, windowSize);
//			}			
//			sb.append(x + "\n");
//		}
//		return sb.toString(); // NB: Leave last trailing \n
//	}
//	
//	@Override
//	/** Write the features in interval to file by appending to existing file. 
//	 * If the file to write to null or empty, return the data that would be
//	 * written as string.
//	 * printFeaturesToFile aims at reproducing the behavior of Linux cat: print to file, possibly appending or to stdout. 
//	 * */
//	public String printFeaturesToFile() throws IOException, InvalidGenomicCoordsException, InvalidColourException {
//		
//		if(this.getExportFile() == null || this.getExportFile().isEmpty()){
//			if(this.isNoFormat()){
//				return this.printFeatures();
//			} else {
//				return "\033[38;5;" + Config.get256Color(ConfigKey.foreground) + 
//						";48;5;" + Config.get256Color(ConfigKey.background) + "m" + this.printFeatures();				
//			}
//		}
//		
//		BufferedWriter wr= null;
//		try{
//			wr = new BufferedWriter(new FileWriter(this.getExportFile(), true));
//			for(List<TextRead> x : this.readStack){
//				for(TextRead txr : x){
//					wr.write(txr.getSamRecord().getSAMString() + "\n");
//				}
//			}
//			wr.close();
//		} catch(IOException e){
//			System.err.println("Cannot write to " + this.getExportFile());
//			throw e;
//		}
//		return "";
//	}

	
	/* S e t t e r s   and   G e t t e r s */

	@Override
	public void setAwk(String awk) throws ClassNotFoundException, IOException, InvalidGenomicCoordsException, InvalidRecordException, SQLException {
		
		String awkFunc= ""; 
		try {
			awkFunc= FileUtils.readFileToString(new File(Main.class.getResource("/functions.awk").toURI()));
		} catch (URISyntaxException e) {

		}

		this.awk= awk;
		this.update();
	}
	
	@Override
	public String getAwk(){
		return this.awk;
	}

	
	public boolean isWithReadName() {
		return withReadName;
	}

	public void setWithReadName(boolean withReadName) {
		this.withReadName = withReadName;
	}
	
}
