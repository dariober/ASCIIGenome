package samTextViewer;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.net.MalformedURLException;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Locale;
import java.util.Set;

import org.apache.commons.lang3.StringUtils;

import coloring.Config;
import coloring.ConfigKey;
import exceptions.InvalidColourException;
import exceptions.InvalidCommandLineException;
import exceptions.InvalidGenomicCoordsException;
import faidx.Faidx;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;

/**
 * Class to set up the horizontal axis on screen. 
 * 
 * * Parse and check input genomic coordinates in input.
 * * Map genomic coordinates to screen coordinates.
 * * Display info about window: ruler, coordinates, bp/column
 * 
 * The actual interval can be smaller than windw size. 
 * @author berald01
 */
public class GenomicCoords implements Cloneable {
	
	public static final String SPACER= "-"; // Background character as spacer
	public static final String TICKED= "*"; // Char to use to mark region
	
	private String chrom;
	private Integer from;
	private Integer to;
	private SAMSequenceDictionary samSeqDict; // Can be null
	/** Size of the screen window 
	 * Only user getUserWindowSize() to access it after init at constructor.
	 */
	private String fastaFile= null;
	private String samSeqDictSource= null; // Source of the sequence dictionary. Can be path to file of fasta, bam, genome. Or genome tag e.g. hg19. 
	private byte[] refSeq= null;
	public boolean isSingleBaseResolution= false;
	private int terminalWidth;
	private List<Double> mapping;
	
	/* Constructors */
	public GenomicCoords(String region, int terminalWidth, SAMSequenceDictionary samSeqDict, String fastaFile, boolean verbose) throws InvalidGenomicCoordsException, IOException{
		
		this.setTerminalWidth(terminalWidth);
		
		GenomicCoords xgc= parseStringToGenomicCoords(region);
		
		this.chrom= xgc.getChrom();
		this.from= xgc.getFrom();
		this.to= xgc.getTo();
		
		if(from == null){ 
			from= 1; 
		}
		
		if(to == null){ 
			to= from + this.getTerminalWidth() - 1;
		}
		
		// Check valid input
		if(chrom == null || (to != null && from == null) || (from > to) || from < 1 || to < 1){
			System.err.println("Got: " + chrom + ":" + from + "-" + to);
			InvalidGenomicCoordsException e = new InvalidGenomicCoordsException();
			throw e;
		}
		if(samSeqDict != null && samSeqDict.size() > 0){ // If dict is present, check against it
			if(samSeqDict.getSequence(chrom) == null){
				if(verbose){
					System.err.println("\nCannot find chromosome '" + chrom + "' in sequence dictionary.");
				}
				InvalidGenomicCoordsException e = new InvalidGenomicCoordsException();
				throw e;
			}
		}
		this.samSeqDict= samSeqDict;
		if(this.samSeqDict != null && this.samSeqDict.size() > 0){
			correctCoordsAgainstSeqDict(samSeqDict);
		}
		if(fastaFile != null){
			this.setSamSeqDictFromFasta(fastaFile);
			this.setFastaFile(fastaFile);
		}
		
		this.update();
	}
	
	public GenomicCoords(String region, int terminalWidth, SAMSequenceDictionary samSeqDict, String fastaFile) throws InvalidGenomicCoordsException, IOException{
		this(region, terminalWidth, samSeqDict, fastaFile, true);
	}
	
	GenomicCoords(int terminalWidth) throws InvalidGenomicCoordsException, IOException{ 
		this.terminalWidth= terminalWidth;
	};

	/**Update a bunch of fields when coordinates change*/
	private void update() throws InvalidGenomicCoordsException, IOException {
		// this.setTerminalWindowSize();
		this.setSingleBaseResolution(); // True if one text character corresponds to 1 bp
		this.setRefSeq();
		this.mapping= this.seqFromToLenOut(this.getTerminalWidth());
	}
	
	/* Methods */
	
	/** Set genome dictionary and fasta file ref if available. See 
	 * GenomicCoords.getSamSeqDictFromAnyFile() for available inputs.
	 * @param includeGenomeFile: Should the input data be treated as a genome file?
	 * Set to true only if the input can be a genome file. Other files (bed, vcf, gff) 
	 * look like valid genome file and this can result in wring dictionary.  
	 * */
	public void setGenome(List<String> input, boolean includeGenomeFile) throws IOException {
		
		List<String> cleanList= new ArrayList<String>(); 
		for(String x : input){
			if(x != null  && ! x.trim().isEmpty()){
				cleanList.add(Utils.tildeToHomeDir(x));
			}
		}		
		if(cleanList.size() == 0){
			return;
		}		

		// Set Dictionary
		this.setSamSeqDictFromAnySource(cleanList, includeGenomeFile);
		// Try to set fasta sequence
		for(String x : cleanList){
			boolean done= true;
			try{
				if(new File(x + ".fai").exists()){
					this.setFastaFile(x);
				} else {
					throw new FileNotFoundException();
				}
//				IndexedFastaSequenceFile fa= new IndexedFastaSequenceFile(new File(x));
//				this.setFastaFile(x);
//				fa.close();
			} catch(FileNotFoundException e){
				try {
					new Faidx(new File(x));
					(new File(x + ".fai")).deleteOnExit();
					this.setFastaFile(x);
				} catch (Exception e1) {
					done= false;
				}
			}
			if(done){
				break;
			}
		}
	}
	
	protected void setRefSeq() throws IOException, InvalidGenomicCoordsException{
		if(this.fastaFile == null ||  ! this.isSingleBaseResolution){
			this.refSeq= null;
			return;
		}
		this.refSeq= this.getSequenceFromFasta();
	}
	
	public byte[] getRefSeq() throws IOException, InvalidGenomicCoordsException {
		if(this.refSeq == null){
			this.setRefSeq();
		}
		return this.refSeq;
	}
	
	
	public byte[] getSequenceFromFasta() throws IOException{
		IndexedFastaSequenceFile faSeqFile = null;
		try {
			faSeqFile = new IndexedFastaSequenceFile(new File(this.fastaFile));
			try{
				byte[] seq= faSeqFile.getSubsequenceAt(this.chrom, this.from, this.to).getBases();
				faSeqFile.close();
				return seq;
			} catch (NullPointerException e){
				System.err.println("Cannot fetch sequence " + this.chrom + ":" + this.from + "-" + this.to +
						" for fasta file " + this.fastaFile);
				e.printStackTrace();
			}
			faSeqFile.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		return null;	
	}
	
	/**
	 * Parse string to return coordinates. This method simply populates the fields chrom, from, to by 
	 * parsing the input string.  
	 * This object can't be used as such as there is no check for valid input. It should be used only by the constructor.
	 * @return
	 * @throws InvalidGenomicCoordsException 
	 * @throws IOException 
	 */
	private GenomicCoords parseStringToGenomicCoords(String x) throws InvalidGenomicCoordsException, IOException{

		Integer from= 1; // Default start/end coords.
		Integer to= this.getTerminalWidth();
		
		GenomicCoords xgc= new GenomicCoords(this.getTerminalWidth());
		
		if(x == null || x.isEmpty()){
			x= "Undefined_contig";
		}
		
		x= x.trim();
		int nsep= StringUtils.countMatches(x, ":");
		if(nsep == 0){ // Only chrom present. It will not handle well chrom names containing ':'
			xgc.chrom= x.trim();
			xgc.from= from;
			xgc.to= to;
			if(xgc.samSeqDict != null && xgc.to > samSeqDict.getSequence(xgc.chrom).getSequenceLength()){
				xgc.to= samSeqDict.getSequence(xgc.chrom).getSequenceLength(); 
			}
		} else {

			xgc.chrom= StringUtils.substringBeforeLast(x, ":").trim();
			
			// Strip chromosome name, remove commas from integers.
			String fromTo= StringUtils.substringAfterLast(x, ":").replaceAll(",", "").trim();
			nsep= StringUtils.countMatches(fromTo, "-");
			if(nsep == 0){ // Only start position given
				xgc.from= Integer.parseInt(StringUtils.substringBefore(fromTo, "-").trim());
				xgc.to= xgc.from + this.getTerminalWidth() - 1;
			} else if(nsep == 1){ // From and To positions given.
				xgc.from= Integer.parseInt(StringUtils.substringBefore(fromTo, "-").trim());
				xgc.to= Integer.parseInt(StringUtils.substringAfter(fromTo, "-").trim());
			} else {
				InvalidGenomicCoordsException e = new InvalidGenomicCoordsException();
				System.err.println("\nUnexpected format for region " + x + "\n");
				e.printStackTrace();
				throw e;
			}
		}
		return xgc;
	} 

	public void correctCoordsAgainstSeqDict(SAMSequenceDictionary samSeqDict) throws InvalidGenomicCoordsException, IOException{

		if(samSeqDict == null || samSeqDict.size() == 0){
			// Just check start pos
			if (this.from <=0 ){
				this.from= 1;
			}			
			return;
		}
		
		if(this.chrom == null){ // Nothing to do
			return;
		}
		if(samSeqDict.getSequence(this.chrom) == null){ // Not found: Nullify everything
			this.chrom= null;
			this.from= null;
			this.to= null;
			return;
		} 
		// Reset min coords
		if( this.from != null && this.from < 1) {
			this.from= 1;
		}
		// Reset max coords
		if( this.from != null && this.from > samSeqDict.getSequence(this.chrom).getSequenceLength() ) {
			this.from= samSeqDict.getSequence(this.chrom).getSequenceLength() - this.getGenomicWindowSize() + 1;
			if(this.from <= 0){
				this.from= 1;
			}
			this.to= this.from + this.getGenomicWindowSize() - 1;
			if(this.to > samSeqDict.getSequence(this.chrom).getSequenceLength()){
				this.to= samSeqDict.getSequence(this.chrom).getSequenceLength();
			}
		}
		if( this.to != null && this.to > samSeqDict.getSequence(this.chrom).getSequenceLength() ) {			
			this.to= samSeqDict.getSequence(this.chrom).getSequenceLength();
		}
	}
	
	public String toString(){
		int range= this.to - this.from + 1;
		return this.chrom + ":" + this.from + "-" + this.to + "; " + NumberFormat.getNumberInstance(Locale.UK).format(range) + " bp";
	}
	
	/** Return current position in the form chrom:start-end */
	public String toStringRegion(){
		return this.getChrom() + ":" + this.getFrom() + "-" + this.getTo();
	}

	/** Get midpoint of genomic interval 
	 * */
	private int getMidpoint(){
		int range= this.to - this.from + 1;
		if(range % 2 == 1){
			range--;
		}
		int midpoint= range / 2 + this.from;
		return midpoint;
	}
	
	/**
	 * Rescale coords to extend them as in zooming-in/-out
	 * @throws IOException 
	 * @throws InvalidGenomicCoordsException 
	 */
	public void zoomOut() throws IOException, InvalidGenomicCoordsException{
		// If window size is 1 you need to extend it otherwise zoom will have no effect!
		if((this.to - this.from) == 0){
			if((this.from - 1) > 0){
				// Try to extend left by 1 bp:
				this.from -= 1;
			} else {
				// Else extend right
				this.to += 1; // But what if you have a chrom of 1bp?!
			}
		}
		
		int zoom= 1;
		// * Get size of window (to - from + 1)
		int range= this.to - this.from + 1;
		if(range % 2 == 1){
			range--;
		}
		int midpoint= this.getMidpoint();

		// Extend midpoint right		
		long zoomTo= midpoint + ((long)range * (long)zoom);
		
		if(zoomTo >= Integer.MAX_VALUE){
			System.err.println("Invalid 'to' coordinate to fetch " + zoomTo + " (integer overflow?)");
			zoomTo= Integer.MAX_VALUE;
		}
		this.to= (int)zoomTo;
		
		// * Extend midpoint left by window size x2 and check coords
		this.from= midpoint - (range * zoom);
		this.from= (this.from <= 0) ? 1 : this.from; 
		if(this.samSeqDict != null && this.samSeqDict.size() > 0){
			if(this.samSeqDict.getSequence(this.chrom).getSequenceLength() > 0){
				this.to= (this.to > this.samSeqDict.getSequence(this.chrom).getSequenceLength()) ? 
						this.samSeqDict.getSequence(this.chrom).getSequenceLength() : this.to;
			}
		}
		this.update();
	}

	/**
	 * Zoom into range. 
	 * @throws IOException 
	 * @throws InvalidGenomicCoordsException 
	 */
	public void zoomIn() throws IOException, InvalidGenomicCoordsException{
		float zoom= (float) (1/4.0);
		// * Get size of window (to - from + 1)
		int range= this.to - this.from + 1;
		if(range % 2 == 1){
			range--;
		}
		// * Get midpoint of range
		int midpoint= this.getMidpoint();
		int extendBy= (int) Math.rint(range * zoom);
		int newFrom= midpoint - extendBy;
		int newTo= midpoint + extendBy;
		if((newTo - newFrom + 1) < this.getUserWindowSize()){ // Reset new coords to be at least windowSize in span
			int diff= this.getUserWindowSize() - (newTo - newFrom + 1);
			if(diff % 2 == 0){
				newFrom -= diff/2;
				newTo += diff/2;
			} else {
				newFrom -= diff/2+1;
				newTo += diff/2;
			}
			// Check new coords are no larger then starting values
			newFrom= (newFrom < this.from) ? this.from : newFrom;
			newTo= (newTo > this.to) ? this.to : newTo;
		}
		this.from= newFrom;
		this.to= newTo;
		if(this.from > this.to){ // Not sure this can happen.
			this.to= this.from;
		}
		this.update();
	}
	
	/** Move coordinates to the left hand side of the current window
	 * @throws IOException 
	 * @throws InvalidGenomicCoordsException 
	 * */
	public void left() throws InvalidGenomicCoordsException, IOException{
		int w= this.getUserWindowSize();
		this.to= this.getMidpoint();
		if((this.to - this.from) < w){
			this.to += (w - (this.to - this.from) - 1);  
		}
		this.update();
	}
	
	/** Move coordinates to the right hand side of the current window
	 * @throws IOException 
	 * @throws InvalidGenomicCoordsException 
	 * */
	public void right() throws InvalidGenomicCoordsException, IOException{
		int w= this.getUserWindowSize();
		this.from= this.getMidpoint();
		if((this.to - this.from) < w){
			this.from -= (w - (this.to - this.from) - 1);  
		}
		this.update();
	}
	
	/**
	 * Same as R seq(from to, length.out). See also func in Utils
	 * @return
	 * @throws IOException 
	 * @throws InvalidGenomicCoordsException 
	 */
	//private List<Double> seqFromToLenOut() throws InvalidGenomicCoordsException, IOException {
	//	return seqFromToLenOut(this.getUserWindowSize());
	//}
	
	/**
	 * Same as R seq(from to, length.out). See also func in Utils
	 * @return
	 * @throws IOException 
	 * @throws InvalidGenomicCoordsException 
	 */
	private List<Double> seqFromToLenOut(int size) throws InvalidGenomicCoordsException, IOException {
		
		if(this.getFrom() == null || this.getTo() == null){
			return null;
		}
		
		List<Double> mapping= new ArrayList<Double>();
		
		if(this.from < 1 || this.from > this.to){
			System.err.println("Invalid genome coordinates: from " + this.from + " to " + this.to);
			try {
				throw new InvalidGenomicCoordsException();
			} catch (InvalidGenomicCoordsException e) {
				e.printStackTrace();
				System.exit(1);
			}
		}
		int span= this.to - this.from + 1;
		// If the genomic span is less then screen size, reduce screen size to.
		// If genomic span == screenSize then you have a mapping one to one.
		if(span <= size){ 
			for(int i= this.from; i <= this.to; i++){
				mapping.add((double)i);
			}
			return mapping;
		}
		
		double step= ((double)span - 1)/(size - 1);
		mapping.add((double)this.from);
		for(int i= 1; i < size; i++){
			mapping.add((double)mapping.get(i-1)+step);
		}
		
		// First check last point is close enough to expectation. If so, replace last point with
		// exact desired.
		double diffTo= Math.abs(mapping.get(mapping.size() - 1) - this.to);
		if(diffTo > ((float)this.to * 0.001)){
			System.err.println("Error generating sequence:");
			System.err.println("Last point: " + mapping.get(mapping.size() - 1));
			System.err.println("To diff: " + diffTo);
			System.err.println("Step: " + step);
		} else {
			mapping.set(mapping.size()-1, (double)this.to);
		}
		
		double diffFrom= Math.abs(mapping.get(0) - this.from);		
		if(diffFrom > 0.01 || mapping.size() != size){
			System.err.println("Error generating sequence:");
			System.err.println("Expected size: " + size + "; Effective: " + mapping.size());
			System.err.println("From diff: " + diffFrom);
			System.exit(1);
		}
		return mapping;
	}	
	
	/**Take of care you are working with double precision here. Do not use this method to 
	 * for exact calculations. 
	 * */
	public double getBpPerScreenColumn() throws InvalidGenomicCoordsException, IOException{
		List<Double> mapping = seqFromToLenOut(this.getUserWindowSize());
		double bpPerScreenColumn= (to - from + 1) / (double)mapping.size();
		return bpPerScreenColumn;
	}
	
	/**
	 * Produce a string representation of the current position on the chromosome
	 * @param nDist: Distance between labels   
	 * @throws IOException 
	 * @throws InvalidGenomicCoordsException 
	 * @throws InvalidColourException 
	 * */
	public String getChromIdeogram(int nDist, boolean noFormat) throws InvalidGenomicCoordsException, IOException, InvalidColourException {
		
		if(this.samSeqDict == null || this.samSeqDict.size() == 0){
			return null;
		}
		List<Double> positionMap = null;
		try{
			positionMap = Utils.seqFromToLenOut(1, this.samSeqDict.getSequence(this.chrom).getSequenceLength(), 
					this.getUserWindowSize());
		} catch (NullPointerException e){
			throw new InvalidGenomicCoordsException();
		}
		// This code taken from printableRuler() above.
		String numberLine= "";
    	int prevLen= 0;
    	int j= 0;
		while(j < positionMap.size()){
			int num= (int)Math.rint(Utils.roundToSignificantFigures(positionMap.get(j), 2));
			String posMark= Utils.parseIntToMetricSuffix(num); // String.valueOf(num);
			if(j == 0){
				numberLine= posMark;
				j += posMark.length();
			} else if((numberLine.length() - prevLen) >= nDist){
				prevLen= numberLine.length();
				numberLine= numberLine + posMark;
				j += posMark.length();
			} else {
				numberLine= numberLine + SPACER;
				j++;
			}
		}
		List<String> map= new ArrayList<String>();
		for(int i= 0; i < numberLine.length(); i++){
			map.add(numberLine.charAt(i) +"");
		}
		
		// ------------------
		
		int fromTextPos= Utils.getIndexOfclosestValue(this.from, positionMap);
		int toTextPos= Utils.getIndexOfclosestValue(this.to, positionMap);

		boolean isFirst= true;
		int lastTick= -1;
		for(int i= fromTextPos; i <= toTextPos; i++){
			if(isFirst || map.get(i).equals(SPACER)){
				map.set(i, TICKED);
				isFirst= false;
			}
			lastTick= i;
		}
		map.set(lastTick, TICKED);
		String ideogram= StringUtils.join(map, "");
		if(ideogram.length() > this.getUserWindowSize()){
			ideogram= ideogram.substring(0, this.getUserWindowSize());
		}
		if(!noFormat){
			ideogram= "\033[48;5;" + Config.get256Color(ConfigKey.background) + ";38;5;" + Config.get256Color(ConfigKey.chrom_ideogram) + "m" + ideogram;
		}
		return ideogram;
	}
	
	/** For debugging only 
	 * @throws IOException 
	 * @throws InvalidGenomicCoordsException */
	public String toStringVerbose(int windowSize) throws InvalidGenomicCoordsException, IOException{
		List<Double> mapping = seqFromToLenOut(this.getUserWindowSize());
		String str= "Genome coords: " + from + "-" + to 
				+ "; screen width: " + mapping.size()
				+ "; scale: " + this.getBpPerScreenColumn() + " bp/column" 
				+ "; Mapping: " + mapping;
		str += "\n";
		str += this.toString();
		return str;
	}
	
	public String printableGenomicRuler(int markDist, boolean noFormat) throws InvalidGenomicCoordsException, IOException, InvalidColourException{
		List<Double> mapping = this.seqFromToLenOut(this.getUserWindowSize());
		String numberLine= this.printRulerFromList(mapping, markDist, 0);
		numberLine= numberLine.replaceFirst("^0 ", "1 "); // Force to start from 1
		numberLine= numberLine.substring(0, this.getUserWindowSize());
		if(!noFormat){
			numberLine= "\033[48;5;" + Config.get256Color(ConfigKey.background) + 
					";38;5;" + Config.get256Color(ConfigKey.ruler) +
					"m" + numberLine;
		}
    	return numberLine;
    }

	public String printablePercentRuler(int markDist, boolean noFormat) throws InvalidGenomicCoordsException, IOException, InvalidColourException{
		// List<Double> mapping =  Utils.seqFromToLenOut(1, this.getUserWindowSize(), this.getUserWindowSize());
		List<Double> mapping =  Utils.seqFromToLenOut(0, 1, this.getUserWindowSize());
		String numberLine= this.printRulerFromList(mapping, markDist, 2);
		numberLine= numberLine.substring(0, this.getUserWindowSize());
		if(!noFormat){
			numberLine= "\033[48;5;" + Config.get256Color(ConfigKey.background) + 
					";38;5;" + Config.get256Color(ConfigKey.ruler) +
					"m" + numberLine;
		}
    	return numberLine;
    }
		
	private String printRulerFromList(List<Double> marks, int markDist, int digits) throws InvalidGenomicCoordsException, IOException {
		int prevLen= 0;
    	int i= 0;
    	// First round numbers and see if we can round digits
    	List<String>rMarks= new ArrayList<String>();
    	markLoop:
    	for(int sfx : new int[]{1000000, 100000, 10000, 1000, 100, 10, 5, 1}){
    		rMarks.clear();
	    	for(Double mark : marks){
	    		double n = Utils.round(mark/sfx, digits) * sfx;
	    		String str= String.format("%." + digits + "f", n);
	    		str= str.replaceAll("^0\\.00", "0"); // Strip leading zero if any.
	    		str= str.replaceAll("^0\\.", "."); 
	    		str= str.replaceAll("\\.00$", ""); // Change x.00 to x. E.g. 1.00 -> 1
	    		rMarks.add(str);
	    	}
	    	Set<String> uniq= new HashSet<String>(rMarks);
	    	if(uniq.size() == marks.size()){
	    		// No duplicates after rounding
	    		break markLoop;
	    	}
    	}

    	StringBuilder numberLine= new StringBuilder();
		while(i < rMarks.size()){
			String label= rMarks.get(i);
			
			if(label.length() >= markDist){
				// Increase markDist if the number is bigger than the space itself
				markDist= label.length() + 1;
			}
			if(i == 0){
				numberLine.append(label);
				i += label.length();
			} else if((numberLine.length() - prevLen) >= markDist){
				prevLen= numberLine.length();
				numberLine.append(label);
				i += label.length();
			} else {
				numberLine.append(" ");
				i++;
			}
		}
		return numberLine.toString();
	}
	
	/** Ref sequence usable for print on screen. 
	 * @throws IOException 
	 * @throws InvalidGenomicCoordsException 
	 * @throws InvalidColourException */
	public String printableRefSeq(boolean noFormat) throws IOException, InvalidGenomicCoordsException, InvalidColourException{

		byte[] refSeq= this.getRefSeq();
		if(refSeq == null){
			return "";
		}
		
		if(noFormat){
			return new String(refSeq) + "\n";
		} else {
			String faSeqStr= "";
			for(byte c : refSeq){
				// For colour scheme see http://www.umass.edu/molvis/tutorials/dna/atgc.htm
				char base= (char) c;
				String prefix= "\033[48;5;" + Config.get256Color(ConfigKey.background) + ";38;5;";
				if(base == 'A' || base == 'a'){
					faSeqStr += prefix + Config.get256Color(ConfigKey.seq_a) + "m" + base;
				} else if(base == 'C' || base == 'c') {
					faSeqStr += prefix + Config.get256Color(ConfigKey.seq_c) + "m" + base;
				} else if(base == 'G' || base == 'g') {
					faSeqStr += prefix + Config.get256Color(ConfigKey.seq_g) + "m" + base;
				} else if(base == 'T' || base == 't') {
					faSeqStr += prefix + Config.get256Color(ConfigKey.seq_t) + "m" + base;
				} else {
					faSeqStr += prefix + Config.get256Color(ConfigKey.seq_other) + "m" + base;
				} 
			}
			return faSeqStr + "\n";
		}
	}

	/**
	 * Get a SAMSequenceDictionary by querying the input files, if there is any bam, or from the indexed fasta file or from genome files.
	 * @param testfiles List of input files
	 * @param fasta Reference sequence
	 * @param genome 
	 * @return
	 * @throws IOException 
	 */
	private boolean setSamSeqDictFromAnySource(List<String> testfiles, boolean includeGenomeFile) throws IOException{

		boolean isSet= false;
		for(String testfile : testfiles){ // Get sequence dict from bam, if any				
			try{
				isSet= this.setSamSeqDictFromBam(testfile);
				if(isSet){
					return isSet;
				}
			} catch(Exception e){
				//
			}
			try{
				isSet= this.setSamSeqDictFromFasta(testfile);
				if(isSet){
					return isSet;
				}
			} catch(Exception e){
				//
			}
			if(includeGenomeFile){
				try{
					isSet= this.setSamSeqDictFromGenomeFile(testfile);
					if(isSet){
						return isSet;
					}
				} catch(Exception e){
					//
				}
			}
		}
		// If we still haven't found a sequence dict, try looking for a VCF file.
		for(String testfile : testfiles){ 
			try{
				isSet= this.setSamSeqDictFromVCF(testfile);
				if(isSet){
					return isSet;
				}
			} catch(Exception e){
				//
			}
		}
		return isSet;
	}

	private boolean setSamSeqDictFromVCF(String vcf) throws MalformedURLException {

		SAMSequenceDictionary samSeqDict= Utils.getVCFHeader(vcf).getSequenceDictionary();
		if(samSeqDict != null){
			this.setSamSeqDictSource(new File(vcf).getAbsolutePath());
			this.setSamSeqDict(samSeqDict);
			return true;
		}
		return false;

	}

	
	private boolean setSamSeqDictFromFasta(String fasta) throws IOException{
		
		// IndexedFastaSequenceFile fa= null;
		
		try{
			if(new File(fasta + ".fai").exists() && ! new File(fasta + ".fai").isDirectory()){
				//
			} else {
				throw new FileNotFoundException(); 
			}
			// fa= new IndexedFastaSequenceFile(new File(fasta));
		} catch(FileNotFoundException e){
			try {
				new Faidx(new File(fasta));
				(new File(fasta + ".fai")).deleteOnExit();
			} catch (Exception e1) {
				//
			}
		}		
		
		BufferedReader br= new BufferedReader(new FileReader(new File(fasta + ".fai")));
		SAMSequenceDictionary seqDict= new SAMSequenceDictionary(); // null;
		while(true){
			String line= br.readLine();
			if(line == null){
				break;
			}
			SAMSequenceRecord ssqRec= new SAMSequenceRecord(
					line.split("\t")[0], 
					Integer.parseInt(line.split("\t")[1]));
			seqDict.addSequence(ssqRec);
		}
		br.close();
//			fa.close();
		this.setSamSeqDictSource(new File(fasta).getAbsolutePath());
		this.setSamSeqDict(seqDict);
		return true;
//		}
//		fa.close();
//		return false;
	}
	
	private boolean setSamSeqDictFromBam(String bamfile) {

		/*  ------------------------------------------------------ */
		/* This chunk prepares SamReader from local bam            */
		SamReaderFactory srf=SamReaderFactory.make();
		srf.validationStringency(ValidationStringency.SILENT);
		SamReader samReader;
		samReader= srf.open(new File(bamfile));
		/*  ------------------------------------------------------ */
		
		SAMSequenceDictionary seqDict = samReader.getFileHeader().getSequenceDictionary();
		if(seqDict != null && !seqDict.isEmpty()){
			this.setSamSeqDictSource(new File(bamfile).getAbsolutePath());
			this.setSamSeqDict(seqDict);
			return true;
		}
		return false;
	}
	
	/** Get SamSequenceDictionary either from local file or from built-in resources.
	 * If reading from resources, "genome" is the tag before '.genome'. E.g. 
	 * 'hg19' will read file hg19.genome  
	 * */
	private boolean setSamSeqDictFromGenomeFile(String genome) throws IOException {

		SAMSequenceDictionary samSeqDict= new SAMSequenceDictionary();

		BufferedReader reader=null;
		try{
			// Attempt to read from resource
			InputStream res= Main.class.getResourceAsStream("/genomes/" + genome + ".genome");
			reader= new BufferedReader(new InputStreamReader(res));
		} catch (NullPointerException e){
			try{
				// Read from local file
				reader= new BufferedReader(new FileReader(new File(genome)));
			} catch (FileNotFoundException ex){
				return false; 
			}
		}
		
		String line = null;
		while ((line = reader.readLine()) != null) {
			if(line.trim().isEmpty() || line.trim().startsWith("#")){
				continue;
			}
			String[] chromLine= line.split("\t");
			SAMSequenceRecord sequenceRecord = new SAMSequenceRecord(chromLine[0], Integer.parseInt(chromLine[1]));
			samSeqDict.addSequence(sequenceRecord);
		}
		reader.close();
		this.setSamSeqDictSource(genome);
		this.setSamSeqDict(samSeqDict);
		return true;
	}

	public boolean equalCoords(GenomicCoords other){
		return this.chrom.equals(other.chrom) && this.from.equals(other.from) && this.to.equals(other.to); 
	}
	
	/** True if all fileds in this object equla those in the other
	 * @throws IOException 
	 * @throws InvalidGenomicCoordsException */
	public boolean equalCoordsAndWindowSize(GenomicCoords other) throws InvalidGenomicCoordsException, IOException{
		return this.equalCoords(other) && 
				this.getUserWindowSize() == other.getUserWindowSize();
	}
	
	public Object clone() {
	//shallow copy
		try {
			return super.clone();
		} catch (CloneNotSupportedException e) {
			return null;
		}
	}

	public void centerAndExtendGenomicCoords(GenomicCoords gc, int size, double slop) throws InvalidGenomicCoordsException, IOException {
		
		if(size <= 0){
			System.err.println("Invalid feature size. Must be > 0, got " + size);
			throw new InvalidGenomicCoordsException();
		}

		if(slop > 0){
			double center= (size/2.0) + gc.getFrom();
			gc.from= (int)Math.rint(center - (size * slop));
			gc.to= (int)Math.rint(center + (size * slop));
		}
		if(slop == 0){
			// Decrease gc.from so that the start of the feature is in the middle of the screen
			int newFrom= gc.from - (int)Utils.getTerminalWidth() / 2;
			newFrom= newFrom < 1 ? 1 : newFrom;
			int newTo= newFrom + Utils.getTerminalWidth() - 1;
			gc.from= newFrom;
			gc.to= newTo;
		}
		if(((gc.to - gc.from)+1) < gc.getUserWindowSize()){
			int span= (gc.to - gc.from);
			int extendBy= (int)Math.rint((gc.getUserWindowSize() / 2.0) - (span / 2.0));
			gc.from -= extendBy;
			gc.to += extendBy;
		}
		gc.correctCoordsAgainstSeqDict(samSeqDict);
		this.update();
	}

	/** Reset window size according to current terminal screen. 
	 * If the user reshapes the terminal window size or the font size, 
	 * detect the new size and add it to the history. 
	 * */
	private int getTerminalWidth() {
		return this.terminalWidth;
	}	
	
	/* Getters and setters */

	public List<Double> getMapping() {
		return this.mapping;
	}
	
	/** Map using this.getUserWindowSize() as window size. Consider using 
	 * getMapping(int size) to avoid computing the terminal width for each call. */
	//public List<Double> getMapping() throws InvalidGenomicCoordsException, IOException {
	//	return seqFromToLenOut(this.getUserWindowSize());
	//}	
	
	public String getChrom() {
		return chrom;
	}

	public Integer getFrom() {
		return from;
	}

	public Integer getTo() {
		return to;
	}
	
	/** Width of the terminal screen window in number of characters. 
	 * Not to be confused with the genomic window size (as in bp) 
	 * @throws IOException 
	 * @throws InvalidGenomicCoordsException */
	public int getUserWindowSize() throws InvalidGenomicCoordsException, IOException{

		int userWindowSize= this.getTerminalWidth();
		if(userWindowSize > this.getGenomicWindowSize()){
			return this.getGenomicWindowSize();
		} else {
			return userWindowSize;
		}
	}
	
	/** Size of genomic interval. Can be smaller than windowSize set by user. 
	 * Not to be confused with userWindowSize. 
	 * @throws IOException 
	 * @throws InvalidGenomicCoordsException */
	public Integer getGenomicWindowSize() throws InvalidGenomicCoordsException {
		if(this.getTo() == null || this.getFrom() == null){
			return null;
		}
		return this.to - this.from + 1;
	}

	public SAMSequenceDictionary getSamSeqDict(){
		return this.samSeqDict;
	}
	
	public void setSamSeqDict(SAMSequenceDictionary samSeqDict){
		this.samSeqDict= samSeqDict;
	}
	
	public String getFastaFile(){
		return this.fastaFile;
	}

	protected void setFastaFile(String fastaFile) {
		this.fastaFile = fastaFile;
	}

	public String getSamSeqDictSource() {
		if(this.getFastaFile() != null){
			return this.getFastaFile();
		}
		return samSeqDictSource;
	}

	public void setSamSeqDictSource(String samSeqDictSource) {
		this.samSeqDictSource = samSeqDictSource;
	}

	/** Extend genomic coordinates left and right by given bases. 
	 * refpoint= "mid": The new coordinates are given by the midpoint of the current one +/- left and right.
	 * refpoint= "window": The new coordinates are the given by the current window extended left and right. 
	 * @throws IOException 
	 * @throws InvalidGenomicCoordsException 
	 * */
	private void extend(int left, int right, String refpoint) throws InvalidGenomicCoordsException, IOException {
		int newFrom= 0;
		int newTo= 0;
		if(refpoint.equals("mid")){
			int mid= this.getMidpoint();
			newFrom= mid - left;
			newTo= mid + right - 1; // -1 so the window size becomes exactly right-left
		}
		if(refpoint.equals("window")){
			newFrom= this.getFrom() - left;
			newTo= this.getTo() + right;
		}
		if(newFrom > newTo){
			int tmp= newFrom;
			newFrom= newTo;
			newTo= tmp;
		}
		this.from= newFrom;
		this.to= newTo;
		this.correctCoordsAgainstSeqDict(this.getSamSeqDict());
		this.update();
	}
	
	/** Apply this.extend() after having parsed cmd line args.  
	 * @throws InvalidCommandLineException 
	 * @throws IOException 
	 * @throws InvalidGenomicCoordsException 
	 * */
	public void cmdInputExtend(List<String> cmdInput) throws InvalidCommandLineException, InvalidGenomicCoordsException, IOException{
		List<String> args= new ArrayList<String>(cmdInput);
		args.remove(0); // Remove cmd name
		
		String refpoint= "window";
		if(args.contains("window")){
			args.remove("window");
		}
		if(args.contains("mid")){
			refpoint= "mid";
			args.remove("mid");
		} 
		if(args.size() == 0){
			throw new InvalidCommandLineException();
		}
		int left= Integer.parseInt(args.get(0));
		int right= left;
		if(args.size() > 1){
			right= Integer.parseInt(args.get(1));
		}
		this.extend(left, right, refpoint);
	}

	public void setTerminalWidth(int terminalWidth) throws InvalidGenomicCoordsException, IOException {
		this.terminalWidth= terminalWidth;
		this.update();
	}

	protected void setSingleBaseResolution() throws InvalidGenomicCoordsException, IOException {
		if(this.getGenomicWindowSize() == null){
			this.isSingleBaseResolution= false;
			return;
		}
		if(this.getUserWindowSize() == this.getGenomicWindowSize()){
			this.isSingleBaseResolution= true;
		} else {
			this.isSingleBaseResolution= false;
		}		
	}

}
