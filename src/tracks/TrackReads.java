package tracks;

import java.io.File;
import java.io.IOException;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Random;

import org.apache.commons.lang3.StringUtils;
import org.apache.commons.validator.routines.UrlValidator;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.filter.AggregateFilter;
import htsjdk.samtools.filter.SamRecordFilter;
import samTextViewer.GenomicCoords;
import samTextViewer.Utils;

/**
 * Prepare track for printing aligned reads 
 * @author berald01
 *
 */
public class TrackReads extends Track{

	private List<List<TextRead>> readStack;
	private boolean bs= false;
	private boolean withReadName= false;
	private int maxReadStack;
	
	/* C o n s t r u c t o r s */
	/**
	 * Create read track
	 * @param sam Bam file to extract reads from
	 * @param gc Interval to consider
	 * @bs Should the track be displayed as BS-Seq data
	 * @param maxReadsStack Accumulate at most this many reads.
	 * @throws IOException 
	 */
	public TrackReads(String bam, GenomicCoords gc, List<SamRecordFilter> filters, int maxReadStack) throws IOException{

		if(!Utils.bamHasIndex(bam)){
			System.err.println("\nAlignment file " + bam + " has no index.\n");
			throw new RuntimeException();
		}
		this.setFilename(bam);
		this.setGc(gc);
		this.setFilters(filters);
		this.maxReadStack= maxReadStack;
		this.update();

	} 
	
	/* M e t h o d s */
	
	public void update() throws MalformedURLException{
		
		this.readStack= new ArrayList<List<TextRead>>();
		if(this.getGc().getGenomicWindowSize() < this.MAX_REGION_SIZE){

			/*  ------------------------------------------------------ */
			/* This chunk prepares SamReader from local bam or URL bam */
			UrlValidator urlValidator = new UrlValidator();
			SamReaderFactory srf=SamReaderFactory.make();
			srf.validationStringency(ValidationStringency.SILENT);
			SamReader samReader;
			if(urlValidator.isValid(this.getFilename())){
				samReader = srf.open(SamInputResource.of(new URL(this.getFilename())).index(new URL(this.getFilename() + ".bai")));
			} else {
				samReader= srf.open(new File(this.getFilename()));
			}
			/*  ------------------------------------------------------ */
			
			long cnt= countReadsInWindow(this.getFilename(), this.getGc(), this.getFilters());
			float probSample= (float) this.maxReadStack / cnt;
			
			Iterator<SAMRecord> sam= samReader.query(this.getGc().getChrom(), this.getGc().getFrom(), this.getGc().getTo(), false);
			List<TextRead> textReads= new ArrayList<TextRead>();
			AggregateFilter aggregateFilter= new AggregateFilter(this.getFilters());
			
			while(sam.hasNext() && textReads.size() < this.maxReadStack){
	
				SAMRecord rec= sam.next();
				if( !aggregateFilter.filterOut(rec) ){
					Random rand = new Random();
					if(rand.nextFloat() < probSample){ // Downsampler
						TextRead tr= new TextRead(rec, this.getGc());
						textReads.add(tr);
					}
				}
			}
			this.readStack= stackReads(textReads);
		} 
	}
	
	/** 
	 * Printable track on screen. This is what should be called by Main */
	@Override
	public String printToScreen(){
		
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
		String printable= "";
		for(Double idx : keep){
			List<TextRead> line= this.readStack.get((int)Math.rint(idx));
			printable += linePrinter(line, bs, this.isNoFormat(), withReadName) + "\n";
		}
		return printable.replaceAll("\n$", "");
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
	 */
	private List<List<TextRead>> stackReads(List<TextRead> textReads){
		
		List<List<TextRead>> listOfLines= new ArrayList<List<TextRead>>();
		if(textReads.size() == 0){
			return listOfLines;
		}
		List<TextRead> line= new ArrayList<TextRead>();
		line.add(textReads.get(0)); 
		textReads.remove(0);
		listOfLines.add(line);
		while(true){
			ArrayList<TextRead> trToRemove= new ArrayList<TextRead>();
			// Find a read in input whose start is greater then end of current
			for(int i=0; i < textReads.size(); i++){
				TextRead tr= textReads.get(i);
				int gap= (this.getGc().getBpPerScreenColumn() > 1) ? 0 : 1; // If reads are very compressed, do not add space between adjacent ones.
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
	 */
	private String linePrinter(List<TextRead> textReads, boolean bs, boolean noFormat, boolean withReadName){
		StringBuilder sb= new StringBuilder();

		int curPos= 0; // Position on the line, needed to pad with blanks btw reads.
		for(TextRead tr : textReads){
			sb.append(StringUtils.repeat(" ", (tr.getTextStart()-1) - curPos));
			String printableRead= tr.getPrintableTextRead(bs, noFormat, withReadName);
			sb.append(printableRead);
			curPos= tr.getTextEnd();
		}
		return sb.toString();
	}

	
	/**
	 * Count reads in interval using the given filters.
	 * NB: Make sure you use the same filter as in readAndStackSAMRecords();
	 * @param bam
	 * @param gc
	 * @param filters List of filters to apply
	 * @return
	 * @throws MalformedURLException 
	 */
	private long countReadsInWindow(String bam, GenomicCoords gc, List<SamRecordFilter> filters) throws MalformedURLException {

		/*  ------------------------------------------------------ */
		/* This chunk prepares SamReader from local bam or URL bam */
		UrlValidator urlValidator = new UrlValidator();
		SamReaderFactory srf=SamReaderFactory.make();
		srf.validationStringency(ValidationStringency.SILENT);
		SamReader samReader;
		if(urlValidator.isValid(bam)){
			samReader = srf.open(SamInputResource.of(new URL(bam)).index(new URL(bam + ".bai")));
		} else {
			samReader= srf.open(new File(bam));
		}
		/*  ------------------------------------------------------ */
		
		long cnt= 0;
		
		//SamReaderFactory srf=SamReaderFactory.make();
		//srf.validationStringency(ValidationStringency.SILENT);
		//SamReader samReader = srf.open(new File(bam));

		Iterator<SAMRecord> sam= samReader.query(gc.getChrom(), gc.getFrom(), gc.getTo(), false);
		AggregateFilter aggregateFilter= new AggregateFilter(filters);
		while(sam.hasNext()){
			SAMRecord rec= sam.next();
			if( !aggregateFilter.filterOut(rec) ){
				cnt++;
			}
		}
		try {
			samReader.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return cnt;
	}
	
	/* S e t t e r s   and   G e t t e r s */
	
	public boolean isBs() {
		return bs;
	}

	public void setBs(boolean bs) {
		this.bs = bs;
	}

	public boolean isWithReadName() {
		return withReadName;
	}

	public void setWithReadName(boolean withReadName) {
		this.withReadName = withReadName;
	}
	
}