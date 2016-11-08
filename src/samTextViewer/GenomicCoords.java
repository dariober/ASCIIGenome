package samTextViewer;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.sql.SQLException;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;
import org.apache.commons.lang3.StringUtils;
import exceptions.InvalidGenomicCoordsException;
import exceptions.InvalidRecordException;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import tracks.TrackWiggles;

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
	
	/* Constructors */
	public GenomicCoords(String region, SAMSequenceDictionary samSeqDict, String fastaFile) throws InvalidGenomicCoordsException, IOException{

		GenomicCoords xgc= parseStringToGenomicCoords(region);

		this.chrom= xgc.getChrom();
		this.from= xgc.getFrom();
		this.to= xgc.getTo();
		
		if(from == null){ 
			from= 1; 
		}
		
		if(to == null){ 
			to= from + this.getTerminalWindowSize() - 1; 
		}
		
		// Check valid input
		if(chrom == null || (to != null && from == null) || (from > to) || from < 1 || to < 1){
			System.err.println("Got: " + chrom + ":" + from + "-" + to);
			InvalidGenomicCoordsException e = new InvalidGenomicCoordsException();
			throw e;
		}
		if(samSeqDict != null && samSeqDict.size() > 0){ // If dict is present, check against it
			if(samSeqDict.getSequence(chrom) == null){
				System.err.println("\nCannot find chromosome '" + chrom + "' in sequence dictionary.");
				InvalidGenomicCoordsException e = new InvalidGenomicCoordsException();
				throw e;
			}
		}
		this.samSeqDict= samSeqDict;
		if(this.samSeqDict != null && this.samSeqDict.size() > 0){
			correctCoordsAgainstSeqDict(samSeqDict);
		}
		this.fastaFile= fastaFile;
		
	}
		
	GenomicCoords(){ 

	};
	
	/* Methods */
	
	/** Set genome dictionary and fasta file ref if available. See 
	 * GenomicCoords.getSamSeqDictFromAnyFile() for available inputs.
	 * */
	public void setGenome(List<String> input) throws IOException {
		
		List<String> cleanList= new ArrayList<String>(); 
		for(String x : input){
			if(x != null  && ! x.trim().isEmpty()){
				cleanList.add(x);
			}
		}		
		if(cleanList.size() == 0){
			return;
		}		

		// Set Dictionary
		this.setSamSeqDictFromAnySource(cleanList);
		
		// Try to set fasta sequence
		for(String x : cleanList){
			try{
				IndexedFastaSequenceFile fa= new IndexedFastaSequenceFile(new File(x));
				if(fa.isIndexed()){
					this.setFastaFile(x);
				}
				fa.close();
			} catch(FileNotFoundException e){
				//
			}
		}
	}

	
	/** Get sequence from fasta, but only if it can fit the screen. Null otherwise. 
	 * @throws InvalidGenomicCoordsException 
	 * */
	public byte[] getRefSeq() throws IOException, InvalidGenomicCoordsException {
		if(this.fastaFile == null || this.getBpPerScreenColumn() > 1){
			return null;
		}
		return this.getSequenceFromFasta(); // this.refSeq;
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
		Integer to= this.getTerminalWindowSize();
		
		GenomicCoords xgc= new GenomicCoords();
		
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
				xgc.to= xgc.from + this.getTerminalWindowSize() - 1;
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
	
	/** Retunr current position in the form chrom:start-end */
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
	 */
	public void zoomOut() throws IOException{
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
		// this.setRefSeq();
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
		// this.setRefSeq();
	}
	
	/** Move coordinates to the left hand side of the current window
	 * @throws IOException 
	 * @throws InvalidGenomicCoordsException 
	 * */
	public void left() throws InvalidGenomicCoordsException, IOException{
		int w= this.getUserWindowSize();
		this.to= this.getMidpoint();
		if((this.to - this.from) < w){
			this.to += (w - (this.to - this.from));  
		}
	}
	
	public void right() throws InvalidGenomicCoordsException, IOException{
		int w= this.getUserWindowSize();
		this.from= this.getMidpoint();
		if((this.to - this.from) < w){
			this.from -= (w - (this.to - this.from));  
		}
		
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
		
		List<Double> mapping= new ArrayList<Double>();
		
		if(from < 1 || from > to){
			System.err.println("Invalid genome coordinates: from " + from + " to " + to);
			try {
				throw new InvalidGenomicCoordsException();
			} catch (InvalidGenomicCoordsException e) {
				e.printStackTrace();
				System.exit(1);
			}
		}
		int span= to - from + 1;
		// If the genomic span is less then screen size, reduce screen size to.
		// If genomic span == screenSize then you have a mapping one to one.
		if(span <= size){ 
			for(int i= from; i <= to; i++){
				mapping.add((double)i);
			}
			return mapping;
		}
		
		double step= ((double)span - 1)/(size - 1);
		mapping.add((double)from);
		for(int i= 1; i < size; i++){
			mapping.add((double)mapping.get(i-1)+step);
		}
		
		// First check last point is close enough to expectation. If so, replace last point with
		// exact desired.
		double diffTo= Math.abs(mapping.get(mapping.size() - 1) - to);
		if(diffTo > ((float)to * 0.001)){
			System.err.println("Error generating sequence:");
			System.err.println("Last point: " + mapping.get(mapping.size() - 1));
			System.err.println("To diff: " + diffTo);
			System.err.println("Step: " + step);
		} else {
			mapping.set(mapping.size()-1, (double)to);
		}
		
		double diffFrom= Math.abs(mapping.get(0) - from);		
		if(diffFrom > 0.01 || mapping.size() != size){
			System.err.println("Error generating sequence:");
			System.err.println("Expected size: " + size + "; Effective: " + mapping.size());
			System.err.println("From diff: " + diffFrom);
			System.exit(1);
		}
		return mapping;
	}	
	
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
	 * */
	public String getChromIdeogram(int nDist, boolean noFormat) throws InvalidGenomicCoordsException, IOException {
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
			ideogram= "\033[30m" + ideogram + "\033[48;5;231m";
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
	
	public String printableRuler(int markDist, boolean noFormat) throws InvalidGenomicCoordsException, IOException{
		List<Double> mapping = this.seqFromToLenOut(this.getUserWindowSize());
    	String numberLine= "";
    	int prevLen= 0;
    	int i= 0;
		while(i < mapping.size()){
			String posMark= String.valueOf(Math.round(mapping.get(i)));
			// Increase markDist if the number is bigger than the space itself
			if(posMark.length() >= markDist){
				markDist= posMark.length() + 1;
			}
			if(i == 0){
				numberLine= posMark;
				i += posMark.length();
			} else if((numberLine.length() - prevLen) >= markDist){
				prevLen= numberLine.length();
				numberLine= numberLine + posMark;
				i += posMark.length();
			} else {
				numberLine= numberLine + " ";
				i++;
			}
		}
		numberLine= numberLine.substring(0, this.getUserWindowSize());
		if(!noFormat){
			numberLine= "\033[30m" + numberLine + "\033[48;5;231m";
		}
    	return numberLine;
    }
	
	/** Ref sequence usable for print on screen. 
	 * @throws IOException 
	 * @throws InvalidGenomicCoordsException */
	public String printableRefSeq(boolean noFormat) throws IOException, InvalidGenomicCoordsException{

		if(this.fastaFile == null || this.getBpPerScreenColumn() > 1){
			return "";
		}
		
		byte[] refSeq= this.getRefSeq();
		
		if(noFormat){
			return new String(refSeq) + "\n";
		} else {
			String faSeqStr= "";
			for(byte c : refSeq){
				// For colour scheme see http://www.umass.edu/molvis/tutorials/dna/atgc.htm
				char base= (char) c;
				if(base == 'A' || base == 'a'){
					faSeqStr += "\033[34m" + base + "\033[48;5;231m";
				} else if(base == 'C' || base == 'c') {
					faSeqStr += "\033[31m" + base + "\033[48;5;231m";
				} else if(base == 'G' || base == 'g') {
					faSeqStr += "\033[32m" + base + "\033[48;5;231m";
				} else if(base == 'T' || base == 't') {
					faSeqStr += "\033[33m" + base + "\033[48;5;231m";
				} else {
					faSeqStr += "\033[30m" + base + "\033[48;5;231m";;
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
	private boolean setSamSeqDictFromAnySource(List<String> testfiles) throws IOException{

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
			try{
				isSet= this.setSamSeqDictFromGenomeFile(testfile);
				if(isSet){
					return isSet;
				}
			} catch(Exception e){
				//
			}
		}
		return isSet;
	}

	private boolean setSamSeqDictFromFasta(String fasta) throws IOException{
		
		IndexedFastaSequenceFile fa= new IndexedFastaSequenceFile(new File(fasta));
		if(fa.isIndexed()){
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
			fa.close();
			this.setSamSeqDictSource(new File(fasta).getAbsolutePath());
			this.setSamSeqDict(seqDict);
			return true;
		}
		fa.close();
		return false;
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
				// System.err.println("\nGenome file not found: " + genome);
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

	/** Return GC profile in region by sliding window of given step size 
	 * @throws IOException 
	 * @throws InvalidRecordException 
	 * @throws InvalidGenomicCoordsException 
	 * @throws SQLException 
	 * @throws ClassNotFoundException 
	 * */
	public TrackWiggles getGCProfile() throws IOException, InvalidRecordException, InvalidGenomicCoordsException, ClassNotFoundException, SQLException {
		if(this.fastaFile == null){
			return null; 
		}
		
		// Read fasta sequence in a single string. 
		// Really there is no need to put the whole thing in memory, but it gives easier implementation.
		IndexedFastaSequenceFile faSeqFile = null;
		faSeqFile = new IndexedFastaSequenceFile(new File(this.fastaFile));
		String fa = new String(faSeqFile.getSubsequenceAt(this.chrom, this.from, this.to).getBases());
		faSeqFile.close();
		
		// Determine step size:
		int step= (int) (this.getBpPerScreenColumn() < 10 ? 10 : Math.rint(this.getBpPerScreenColumn()));
		
		// * Write a temp bedgraph file with coordinates given by step
		File temp = File.createTempFile("gcProfile.", ".tmp.bedgraph");
		temp.deleteOnExit();
		BufferedWriter bw= new BufferedWriter(new FileWriter(temp)); 
		double ymin= 100;
		double ymax= 0;
		for(int xfrom= 0; xfrom < fa.length(); xfrom += step){
			int xto= xfrom + step;
			// Look ahead: If the next round is going to be the last one and the step
			// is not a multiple of the sequence length, merge this round with the next in a 
			// longer sequence. This is to avoid the last chunk to be too skewed in percentage.
			if(((xto+step) > fa.length())){
				xto= fa.length();
			}
			String chunk= fa.substring(xfrom, xto).toUpperCase();

			int gcCnt= StringUtils.countMatches(chunk, 'C') + StringUtils.countMatches(chunk, 'G');
			float pct= (float) ((float)gcCnt / chunk.length() * 100.0);
			if(pct < ymin){
				ymin= pct;
			} 
			if(pct > ymax){
				ymax= pct;
			} 
			// Shouldn't use + to concatenate...
			String line= this.chrom + "\t" + (this.from + xfrom - 1) + "\t" + (this.from + xto - 1) + "\t" + pct;
			bw.write(line + "\n");
			if(xto >= fa.length()){
				break;
			}
		}
		bw.close();
		// * Read this bedGraph via TrackWiggles 
		TrackWiggles cgWiggle= new TrackWiggles(temp.getAbsolutePath(), this, 4);
		temp.delete(); // Track has been created - remove al associated tmp files.
		new File(cgWiggle.getFilename()).delete();
		new File(cgWiggle.getFilename() + ".tbi").delete();

//		cgWiggle.setTrackTag(GenomicCoords.gcProfileFileTag);
		cgWiggle.setYLimitMin(0);
		cgWiggle.setYLimitMax(100);
		
		return cgWiggle;
	}


	public void centerAndExtendGenomicCoords(GenomicCoords gc, int size, double slop) throws InvalidGenomicCoordsException, IOException {
		
		if(size <= 0){
			System.err.println("Invalid feature size. Must be > 0, got " + size);
			throw new InvalidGenomicCoordsException();
		}
				
		double center= (size/2.0) + gc.getFrom();
		gc.from= (int)Math.rint(center - (size * slop));
		gc.to= (int)Math.rint(center + (size * slop));
		
		if(((gc.to - gc.from)+1) < gc.getUserWindowSize()){
			int span= (gc.to - gc.from);
			int extendBy= (int)Math.rint((gc.getUserWindowSize() / 2.0) - (span / 2.0));
			gc.from -= extendBy;
			gc.to += extendBy;
		}
		gc.correctCoordsAgainstSeqDict(samSeqDict);
		
	}

	/** Reset window size according to current terminal screen. 
	 * If the user reshapes the terminal window size or the font size, 
	 * detect the new size and add it to the history. 
	 * */
	private int getTerminalWindowSize() throws InvalidGenomicCoordsException, IOException{
		return jline.TerminalFactory.get().getWidth() - 1;
	}	

	
	/* Getters and setters */

	public List<Double> getMapping(int size) throws InvalidGenomicCoordsException, IOException {
		return seqFromToLenOut(size);
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

		int userWindowSize= this.getTerminalWindowSize();
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
	public int getGenomicWindowSize() throws InvalidGenomicCoordsException, IOException{
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
		this.setSamSeqDictSource(fastaFile);
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
}
