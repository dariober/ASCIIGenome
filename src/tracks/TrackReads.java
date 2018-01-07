package tracks;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;
import java.util.ListIterator;
import java.util.Random;
import java.util.regex.Pattern;

import org.apache.commons.lang3.StringUtils;

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

	private List<List<SamSequenceFragment>> readStack;
	// private boolean withReadName= false;
	private long nRecsInWindow= -1;
	private int userWindowSize;
	private List<Argument> colorForRegex= null;
	private long alnRecCnt= -1;
	
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
	}
	
	protected TrackReads(){
		this.setTrackFormat(TrackFormat.BAM);
	};
	
	/* M e t h o d s */
	
	public void update() throws InvalidGenomicCoordsException, IOException{

		if(this.getyMaxLines() == 0){
			return;
		}
		
		this.userWindowSize= this.getGc().getUserWindowSize();
		
		this.readStack= new ArrayList<List<SamSequenceFragment>>();
		if(this.getGc().getGenomicWindowSize() < this.MAX_REGION_SIZE){

			SamReader samReader= Utils.getSamReader(this.getWorkFilename());
			List<Boolean> passFilter= this.filterReads(samReader, this.getGc().getChrom(), this.getGc().getFrom(), this.getGc().getTo());

			this.nRecsInWindow= 0;
			for(boolean x : passFilter){ // The count of reads in window is the count of reads passing filters
				if(x){
					this.nRecsInWindow++;
				}
			}
			samReader= Utils.getSamReader(this.getWorkFilename());
			Iterator<SAMRecord> sam= samReader.query(this.getGc().getChrom(), this.getGc().getFrom(), this.getGc().getTo(), false);
			
			float max_reads= Float.parseFloat(Config.get(ConfigKey.max_reads_in_stack));
			float probSample= max_reads / this.nRecsInWindow;
			
			// Add this random String to the read name so different screenshot will generate 
			// different samples. 
			String rndOffset= Integer.toString(new Random().nextInt());

			List<TextRead> textReads= new ArrayList<TextRead>();
			ListIterator<Boolean> pass = passFilter.listIterator();
			while(sam.hasNext() && textReads.size() < max_reads){
				SAMRecord rec= sam.next();
				if( pass.next() ){
					String templ_name= Utils.templateNameFromSamReadName(rec.getReadName());
					long v= (templ_name + rndOffset).hashCode(); // Hashing.md5().hashBytes((templ_name + rndOffset).getBytes()).asLong();
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
			keep= Utils.seqFromToLenOut(0, yMaxLines, yMaxLines);
			// keep= Utils.seqFromToLenOut(0, this.readStack.size()-1, yMaxLines);
		} else {
			keep= Utils.seqFromToLenOut(0, this.readStack.size()-1, this.readStack.size());
		}
		StringBuilder printable= new StringBuilder();
		// this.changeFeatureColor(null);
		for(Double idx : keep){
			List<SamSequenceFragment> line= this.readStack.get((int)Math.rint(idx));
			try {
				printable.append(linePrinter(line, this.bisulf, this.isNoFormat()));
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
	private List<List<SamSequenceFragment>> stackReads(List<TextRead> textReads) throws InvalidGenomicCoordsException, IOException{
		
		List<List<SamSequenceFragment>> listOfLines= new ArrayList<List<SamSequenceFragment>>();
		if(textReads.size() == 0){
			return listOfLines;
		}

		List<SamSequenceFragment> fragments= this.makeFragments(textReads, this.getReadsAsPairs()); 
		List<SamSequenceFragment> line= new ArrayList<SamSequenceFragment>(); // All fragments going in the same line
		line.add(fragments.get(0)); 
		fragments.remove(0);
		listOfLines.add(line);
		int gap= (this.getGc().isSingleBaseResolution) ? 1 : 0; // If reads are very compressed, do not add space between adjacent ones.
		while(true){
			ArrayList<SamSequenceFragment> fragToRemove= new ArrayList<SamSequenceFragment>();
			// Find a fragment in input whose start is greater then end of current
			for(int i=0; i < fragments.size(); i++){
				SamSequenceFragment frag= fragments.get(i);
				if(frag.getTextStart() > line.get(line.size()-1).getTextEnd()+gap){ // +2 because we want some space between adjacent reads
					listOfLines.get(listOfLines.size()-1).add(frag); // Append to the last line. 
					fragToRemove.add(frag);
				}
			} // At the end of the loop you have put in line as many reads as you can. 
			for(SamSequenceFragment tr : fragToRemove){ 
				fragments.remove(fragments.indexOf(tr));
			}
			// Create a new line, add the first textRead in list
			if(fragments.size() > 0){
				line= new ArrayList<SamSequenceFragment>();
				line.add(fragments.get(0));
				listOfLines.add(line);
				fragments.remove(0);
			} else {
				break;
			}
		}
		return listOfLines;
	}
	
	/**Match reads in textReads list to return a list fragments. Fragments are
	 * returned sorted by start position. 
	 * @param paired If true try to match up read pairs in the same fragment. If false, each read
	 * is a fragment so the resulting list is effectively the same as the input.  
	 * */
	private List<SamSequenceFragment> makeFragments(List<TextRead> textReads, boolean asPair) {

		List<SamSequenceFragment> fragments= new ArrayList<SamSequenceFragment>();

		while(textReads.size() > 0){
			// Keep going until all reads have been moved to the list of fragments.
			TextRead tr= textReads.get(0);
			textReads.remove(0);
			if( ! asPair || ! tr.getSamRecord().getProperPairFlag()){
				SamSequenceFragment frag= new SamSequenceFragment(tr);
				if(! asPair){ 
					frag.setSingleton(true);
				}
				fragments.add(frag);
			}
			else {
				// Find the mate of this read, if present.
				TextRead mate= null;
				for(TextRead candidateMate : textReads){
					if(candidateMate.getSamRecord().getProperPairFlag() && 
					   Utils.equalReadNames(tr.getSamRecord().getReadName(), candidateMate.getSamRecord().getReadName()) &&
					   tr.getSamRecord().getAlignmentStart() == candidateMate.getSamRecord().getMateAlignmentStart()){
						mate= candidateMate;
						break;
					}
				} // After this loop either we have found a mate or not. Either way, create a fragment from a singleton or a pair.
				if(mate == null){
					fragments.add(new SamSequenceFragment(tr));
				} else {
					fragments.add(new SamSequenceFragment(tr, mate));
					textReads.remove(mate);
				}
			}
		}
		sortFragmentsByStartPosition(fragments); // This may be redundant.
		return fragments;
	}

    private void sortFragmentsByStartPosition(List<SamSequenceFragment> fragments){
    	  Collections.sort(fragments, new Comparator<SamSequenceFragment>() {
    	      @Override
    	      public int compare(final SamSequenceFragment frag1, final SamSequenceFragment frag2) {
    	          return Integer.compare(frag1.getLeftRead().getSamRecord().getAlignmentStart(), 
    	        		                 frag2.getLeftRead().getSamRecord().getAlignmentStart());
    	      }
    	  });    	
    }
	
	/** Prepare a printable string of each output line. 
	 * @param textReads List reads to print out on the same line.
	 * @param noFormat Do not format reads.
	 * @return
	 * @throws IOException 
	 * @throws InvalidGenomicCoordsException 
	 * @throws InvalidColourException 
	 */
	private String linePrinter(List<SamSequenceFragment> fragments, boolean bs, boolean noFormat) throws IOException, InvalidGenomicCoordsException, InvalidColourException{
		StringBuilder sb= new StringBuilder();

		int curPos= 0; // Position on the line, needed to pad with blanks btw reads.
		for(SamSequenceFragment frag : fragments){
			String line= StringUtils.repeat(" ", (frag.getTextStart()-1) - curPos); // Left pad with spaces
			sb.append(line);
			String printableRead= frag.getPrintableFragment(bs, noFormat);
			sb.append(printableRead);
			curPos= frag.getTextEnd();
		}
		int nchars= Utils.stripAnsiCodes(sb.toString()).length();
		if(nchars < this.userWindowSize){
			sb.append(StringUtils.repeat(' ', this.userWindowSize - nchars));
		}
		return this.highlightMidCharacter(sb.toString());
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
		
		String libsize= "";
		if(this.alnRecCnt != -1){
			libsize= "/" + this.alnRecCnt;
		}
		String xtitle= this.getTrackTag() 
				+ "; Reads: " + this.nRecsInWindow + libsize 
				+ this.getTitleForActiveFilters(); 
		return this.formatTitle(xtitle) + "\n";
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
	protected List<String> getRecordsAsStrings() {
		List<String> featureList= new ArrayList<String>();
		
		for(List<SamSequenceFragment> x : this.readStack){
			for(SamSequenceFragment frag : x){
				featureList.add(frag.getLeftRead().getSamRecord().getSAMString());
				if(frag.getRightRead() != null){
					featureList.add(frag.getRightRead().getSamRecord().getSAMString());	
				}
			}
		}
		return featureList;
	}

	/* S e t t e r s   and   G e t t e r s */
	
	@Override
	public boolean getReadsAsPairs(){
		return this.readsAsPairs;
	}
	
	@Override
	public void setReadsAsPairs(boolean readsAsPairs) throws InvalidGenomicCoordsException, IOException {
		this.readsAsPairs= readsAsPairs;
		this.update();
	}

	@Override
	public void changeFeatureColor(List<Argument> args){
		List<List<SamSequenceFragment>> stack = this.readStack;
		for(List<SamSequenceFragment> frags : stack){
			for(SamSequenceFragment frag : frags){
				List<TextRead> reads= new ArrayList<TextRead>();
				reads.add(frag.getLeftRead());
				if(frag.getRightRead() != null){
					reads.add(frag.getRightRead());	
				}
				for(TextRead tr : reads){
					for(Argument arg : args){
						String regex= arg.getKey();
						String color= arg.getArg();
						boolean matched= Pattern.compile(regex).matcher(tr.getSamRecord().getSAMString()).find();
						if(arg.isInvert()){
							matched= ! matched;
						}
						if(matched){
							// tr.getTextReadAsFeatureChars(this.isBisulf());
							// f.setFgColor(color);
						}
					}
				}
			}
		}
	}
	
	@Override
	protected void setColorForRegex(List<Argument> xcolorForRegex) {
		if(xcolorForRegex == null){
			this.colorForRegex= null;
			return;
		} else {
			if(this.colorForRegex == null){
				this.colorForRegex= new ArrayList<Argument>();
			}
			for(Argument p : xcolorForRegex){
				this.colorForRegex.add(p);
			}
		}
	}

	private List<Argument> getColorForRegex() {
		return this.colorForRegex;
	}

}
