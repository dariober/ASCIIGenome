package tracks;

import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.ListIterator;
import java.util.Random;

import org.apache.commons.lang3.StringUtils;

import com.google.common.hash.Hashing;

import coloring.Config;
import coloring.ConfigKey;
import exceptions.InvalidColourException;
import exceptions.InvalidGenomicCoordsException;
import exceptions.InvalidRecordException;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import samTextViewer.GenomicCoords;
import samTextViewer.Utils;

/**
 * Prepare track for printing aligned reads 
 * @author berald01
 *
 */
public class TrackReads extends Track{

	private List<List<TextRead>> readStack;
	private boolean withReadName= false;
	private long nRecsInWindow= -1;
	private int userWindowSize;
	
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
		
		if(this.getyMaxLines() == 0){
			return;
		}
		
		this.userWindowSize= this.getGc().getUserWindowSize();
		
		this.readStack= new ArrayList<List<TextRead>>();
		if(this.getGc().getGenomicWindowSize() < this.MAX_REGION_SIZE){

			List<TextRead> textReads= new ArrayList<TextRead>();
			
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
			
			float probSample= Float.parseFloat(Config.get(ConfigKey.max_reads_in_stack)) / this.nRecsInWindow;
			
			// Add this random String to the read name so different screenshot will generate 
			// different samples. 
			String rndOffset= Integer.toString(new Random().nextInt());
			                                                            
			ListIterator<Boolean> pass = passFilter.listIterator();
			while(sam.hasNext() && textReads.size() < Float.parseFloat(Config.get(ConfigKey.max_reads_in_stack))){
				SAMRecord rec= sam.next();
				if( pass.next() ){
					// Get read name (template name in fact), w/o the suffixes and what comes after the first blank.
					String templ_name= rec.getReadName().replaceAll(" .*", "").replaceAll("/1$|/2$", "");
					long v= Hashing.md5().hashBytes((templ_name + rndOffset).getBytes()).asLong();
					Random rand = new Random(v);
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
//			for(TextRead tr : line){
//				tr.setIdeogram(tr.getIdeogram(true, false), false);
//				this.changeFeatureColor(this.getColorForRegex()); // For each TextRead in this interval set colour according to regex in List<Argument>
//			}
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
		int gap= (this.getGc().isSingleBaseResolution) ? 1 : 0; // If reads are very compressed, do not add space between adjacent ones.
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
		for(TextRead tr : textReads){
			String line= StringUtils.repeat(" ", (tr.getTextStart()-1) - curPos);
			sb.append(line);
			String printableRead= tr.getPrintableTextRead(bs, noFormat, withReadName);
			sb.append(printableRead);
			curPos= tr.getTextEnd();
		}
		int nchars= Utils.stripAnsiCodes(sb.toString()).length();
		if(nchars < this.userWindowSize){
			sb.append(StringUtils.repeat(' ', this.userWindowSize - nchars));
		}
		return this.highlightMidCharacter(sb.toString()); // highlightMidCharacter(fmtLine);
	}

	/** Find the mid character and add some formatting to highlight it.
	 * */
	private String highlightMidCharacter(String fmtLine){
		if(this.isNoFormat()){
			return fmtLine;
		}
		List<Integer> idx = Utils.indexOfCharsOnFormattedLine(fmtLine);
		if(idx.size() <= 7){
			return fmtLine;
		}
		int mid= idx.get(this.userWindowSize/2);
		char midChar= fmtLine.charAt(mid);
		if(midChar != ' '){
			String hLine= fmtLine.substring(0, mid) + "\033[1;7m" + midChar + "\033[21;27m" + fmtLine.substring(mid+1);
			return hLine;
		}
		return fmtLine;
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
		
		String awk= "";
		if(!this.getAwk().isEmpty()){
			awk= "; awk:on";
		}
		
		String xtitle= this.getTrackTag() 
				+ "; Reads: " + this.nRecsInWindow + "/" + Utils.getAlignedReadCount(this.getWorkFilename()) 
				+ samtools 
				+ awk; 
		return this.formatTitle(xtitle) + "\n";
	}
//	public String getTitle() throws InvalidColourException, InvalidGenomicCoordsException, IOException{
//		
//		if(this.isHideTitle()){
//			return "";
//		}
//		String samtools= "";
//		if( ! (this.get_F_flag() == Track.F_FLAG) ){
//			samtools += " -F " + this.get_F_flag();
//		}
//		if( ! (this.get_f_flag() == Track.f_FLAG) ){
//			samtools += " -f " + this.get_f_flag();
//		}
//		if( ! (this.getMapq() == Track.MAPQ) ){
//			samtools += " -q " + this.getMapq();
//		}
//		if( ! samtools.isEmpty()){
//			samtools= "; samtools" + samtools;
//		}
//		String title= this.getTrackTag()
//				+ samtools;
//		
//		//String xtitle= Utils.padEndMultiLine(title, this.getGc().getUserWindowSize());
//		return this.formatTitle(title) + "\n";
//	}
	
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

	/* S e t t e r s   and   G e t t e r s */

	@Override
	public void setAwk(String awk) throws ClassNotFoundException, IOException, InvalidGenomicCoordsException, InvalidRecordException, SQLException {
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
