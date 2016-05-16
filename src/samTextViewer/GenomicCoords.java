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
import java.net.URL;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;

import org.apache.commons.lang3.StringUtils;
import org.apache.commons.validator.routines.UrlValidator;

import exceptions.InvalidGenomicCoordsException;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import tracks.TrackFormat;
//import readWriteBAMUtils.ReadWriteBAMUtils;
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
	
	private String chrom;
	private Integer from;
	private Integer to;
	private SAMSequenceDictionary samSeqDict; // Can be null
	private int windowSize; // Size of the screen window
	private byte[] refSeq;
	private String fastaFile= null;
	private String gcProfileFileTag= "CG_percent";
	
	/* Constructors */
	public GenomicCoords(String chrom, Integer from, Integer to, SAMSequenceDictionary samSeqDict, int windowSize, String fastaFile) 
			throws InvalidGenomicCoordsException, IOException{
		
		if(from == null){ from= 1; }
		
		if(to == null){ to= from + windowSize - 1; }
		
		// Check valid input
		if(chrom == null || (to != null && from == null) || (from > to) || from < 1 || to < 1){
			System.err.println("Got: " + chrom + ":" + from + "-" + to);
			InvalidGenomicCoordsException e = new InvalidGenomicCoordsException();
			throw e;
		}
		if(samSeqDict != null && samSeqDict.size() > 0){ // If dict is present, check against it
			if(samSeqDict.getSequence(chrom) == null){
				System.err.println("\nCannot find chromosome '" + chrom + "' in sequence dictionary.\n");
				InvalidGenomicCoordsException e = new InvalidGenomicCoordsException();
				throw e;
			}
		}
		if(windowSize < 1){
			InvalidGenomicCoordsException e = new InvalidGenomicCoordsException();
			throw e;			
		}
		this.chrom= chrom;
		this.from= from;
		this.to= to;
		this.windowSize= windowSize;
		this.samSeqDict= samSeqDict;
		if(this.samSeqDict != null && this.samSeqDict.size() > 0){
			correctCoordsAgainstSeqDict(samSeqDict);
		}
		this.fastaFile= fastaFile;
		this.setRefSeq();
	}
	
	public GenomicCoords(String region, SAMSequenceDictionary samSeqDict, int windowSize, String fastaFile) throws InvalidGenomicCoordsException, IOException{
		GenomicCoords xgc= parseStringToGenomicCoords(region, windowSize);
		GenomicCoords gc= new GenomicCoords(xgc.chrom, xgc.from, xgc.to, samSeqDict, windowSize, fastaFile);
		this.chrom= gc.chrom;
		this.from= gc.from;
		this.to= gc.to;
		this.windowSize= gc.windowSize;
		this.samSeqDict= gc.samSeqDict;
		this.fastaFile= fastaFile;
		this.setRefSeq();
	}
		
	private GenomicCoords(){ };
	
	/* Methods */
	
	/** Extract sequence from fasta file. Return null if fastaFile is null or the genomic span is larger than 
	 * windowSize 
	 * @throws IOException */
	private void setRefSeq() throws IOException{
		if(this.fastaFile == null || this.getBpPerScreenColumn() > 1){
			this.refSeq= null;
			return;
		}
		IndexedFastaSequenceFile faSeqFile = null;
		try {
			faSeqFile = new IndexedFastaSequenceFile(new File(this.fastaFile));
			try{
				this.refSeq= faSeqFile.getSubsequenceAt(this.chrom, this.from, this.to).getBases();
			} catch (NullPointerException e){
				System.err.println("Cannot fetch sequence " + this.chrom + ":" + this.from + "-" + this.to +
						" for fasta file " + this.fastaFile);
				e.printStackTrace();
			}
			faSeqFile.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
	}
	
	public byte[] getRefSeq(){
		return this.refSeq;
	}
	/**
	 * Parse string to return coordinates. This method simply populates the fields chrom, from, to by 
	 * parsing the input string.  
	 * This object can't be used as such as there is no check for valid input. It should be used only by the constructor.
	 * @return
	 * @throws InvalidGenomicCoordsException 
	 */
	private GenomicCoords parseStringToGenomicCoords(String x, int windowSize) throws InvalidGenomicCoordsException{

		Integer from= 1; // Default start/end coords.
		Integer to= 1;
		
		GenomicCoords xgc= new GenomicCoords();
		
		x= x.trim();
		int nsep= StringUtils.countMatches(x, ":");
		if(nsep == 0){ // Only chrom present. It will not handle well chrom names containing ':'
			xgc.chrom= x.trim();
			xgc.from= from;
			xgc.to= to + windowSize - 1;
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
				xgc.to= xgc.from + windowSize -1;
			} else if(nsep == 1){
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

	public void correctCoordsAgainstSeqDict(SAMSequenceDictionary samSeqDict){
		
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
	
	private int getMidpoint(){
		int range= this.to - this.from + 1;
		if(range % 2 == 1){
			range--;
		}
		// * Get midpoint of genomic interval
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
		
		// * Extend midpoint left by window size x2 and check coords
		this.from= midpoint - (range * zoom);
		this.from= (this.from <= 0) ? 1 : this.from; 
		
		// Extend midpoint right
		this.to= midpoint + (range * zoom);
		if(this.samSeqDict != null && this.samSeqDict.size() > 0){
			if(this.samSeqDict.getSequence(this.chrom).getSequenceLength() > 0){
				this.to= (this.to > this.samSeqDict.getSequence(this.chrom).getSequenceLength()) ? 
						this.samSeqDict.getSequence(this.chrom).getSequenceLength() : this.to;
			}
		}
		this.setRefSeq();
	}

	/**
	 * Zoom into range. 
	 * @throws IOException 
	 */
	public void zoomIn() throws IOException{
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
		if((newTo - newFrom + 1) < windowSize){ // Reset new coords to be at least windowSize in span
			int diff= windowSize - (newTo - newFrom + 1);
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
		this.setRefSeq();
	}
	
	/**
	 * Same as R seq(from to, length.out). See also func in Utils
	 * @param from
	 * @param to
	 * @param lengthOut Length of sequence, effectively the desired screen width.
	 * @return
	 */
	private List<Double> seqFromToLenOut(){
		
		List<Double> mapping= new ArrayList<Double>();
		
		if(from < 1 || from > to){
			System.err.println("Invalid genome coordinates: from " + from + " to " + to);
			System.exit(1);
		}
		int span= to - from + 1;
		// If the genomic span is less then screen size, reduce screen size to.
		// If genomic span == screenSize then you have a mapping one to one.
		if(span <= this.windowSize){ 
			for(int i= from; i <= to; i++){
				mapping.add((double)i);
			}
			return mapping;
		}
		
		double step= ((double)span - 1)/(this.windowSize - 1);
		mapping.add((double)from);
		for(int i= 1; i < this.windowSize; i++){
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
		if(diffFrom > 0.01 || mapping.size() != this.windowSize){
			System.err.println("Error generating sequence:");
			System.err.println("Expected size: " + this.windowSize + "; Effective: " + mapping.size());
			System.err.println("From diff: " + diffFrom);
			System.exit(1);
		}
		return mapping;
	}
	
	public double getBpPerScreenColumn(){
		List<Double> mapping = seqFromToLenOut();
		double bpPerScreenColumn= (to - from + 1) / (double)mapping.size();
		return bpPerScreenColumn;
	}
	
	/**
	 * Produce a string representation of the current position on the chromosome
	 * @param nDist: Distance between labels   
	 * */
	public String getChromIdeogram(int nDist) {
		String SPACER= "-"; // Background character as spacer
		String TICKED= "*"; // Char to use to mark region
		if(this.samSeqDict == null || this.samSeqDict.size() == 0){
			return null;
		}
		List<Double> positionMap = Utils.seqFromToLenOut(1, this.samSeqDict.getSequence(this.chrom).getSequenceLength(), this.windowSize);
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
		return String.join("", map);
	}
	
	/** For debugging only */
	public String toStringVerbose(int windowSize){
		List<Double> mapping = seqFromToLenOut();
		String str= "Genome coords: " + from + "-" + to 
				+ "; screen width: " + mapping.size()
				+ "; scale: " + this.getBpPerScreenColumn() + " bp/column" 
				+ "; Mapping: " + mapping;
		str += "\n";
		str += this.toString();
		return str;
	}
	
	public String printableRuler(int markDist){
		List<Double> mapping = seqFromToLenOut();
    	String numberLine= "";
    	int prevLen= 0;
    	int i= 0;
		while(i < mapping.size()){
			String posMark= String.valueOf(Math.round(mapping.get(i)));
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
    	return numberLine;	
    }
	
	/** Ref sequence usable for print on screen. */
	public String printableRefSeq(boolean noFormat){

		if(this.refSeq == null){
			return "";
		} else if(noFormat){
			return new String(this.refSeq) + "\n";
		} else {
			String faSeqStr= "";
			for(byte c : this.refSeq){
				// For colour scheme see http://www.umass.edu/molvis/tutorials/dna/atgc.htm
				char base= (char) c;
				if(base == 'A' || base == 'a'){
					faSeqStr += "\033[107;34m" + base + "\033[0m";
				} else if(base == 'C' || base == 'c') {
					faSeqStr += "\033[107;31m" + base + "\033[0m";
				} else if(base == 'G' || base == 'g') {
					faSeqStr += "\033[107;32m" + base + "\033[0m";
				} else if(base == 'T' || base == 't') {
					faSeqStr += "\033[107;33m" + base + "\033[0m";
				} else {
					faSeqStr += base;
				} 
			}
			return faSeqStr + "\n";
		}
	}

	/**
	 * Get a SAMSequenceDictionary by querying the input files, if there is any bam, or from the indexed fasta file or from genome files.
	 * @param insam List of input files
	 * @param fasta Reference sequence
	 * @param genome 
	 * @return
	 * @throws IOException 
	 */
	public static SAMSequenceDictionary getSamSeqDictFromAnyFile(List<String> insam, String fasta, String genome) throws IOException{

		SAMSequenceDictionary seqDict= new SAMSequenceDictionary(); // null;

		if(insam != null){
			for(String x : insam){ // Get sequence dict from bam, if any				
				if(Utils.getFileTypeFromName(x).equals(TrackFormat.BAM)){
					
					/*  ------------------------------------------------------ */
					/* This chunk prepares SamReader from local bam or URL bam */
					UrlValidator urlValidator = new UrlValidator();
					SamReaderFactory srf=SamReaderFactory.make();
					srf.validationStringency(ValidationStringency.SILENT);
					SamReader samReader;
					if(urlValidator.isValid(x)){
						samReader = srf.open(SamInputResource.of(new URL(x)).index(new URL(x + ".bai")));
					} else {
						samReader= srf.open(new File(x));
					}
					/*  ------------------------------------------------------ */
					
					//SamReaderFactory srf=SamReaderFactory.make();
					//srf.validationStringency(ValidationStringency.SILENT);
					//SamReader samReader= srf.open(new File(x));
					
					seqDict= samReader.getFileHeader().getSequenceDictionary();
					if(!seqDict.isEmpty()){
						return seqDict;
					}
				}
			}
		}
		if(fasta != null){ // Try getting from fasta
			IndexedFastaSequenceFile fa= new IndexedFastaSequenceFile(new File(fasta));
			if(fa.isIndexed()){
				BufferedReader br= new BufferedReader(new FileReader(new File(fasta + ".fai")));
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
			}
			fa.close();
		}
		if(genome != null && !genome.isEmpty()){ // Try genome file as last option
			seqDict= getSamSeqDictFromGenomeFile(genome);
			return seqDict;
		}
		return seqDict;
	}
	
	/** Get SamSequenceDictionary either from local file or from built-in resources.
	 * If reading from resources, "genome" is the tag before '.genome'. E.g. 
	 * 'hg19' will read file hg19.genome  
	 * */
	private static SAMSequenceDictionary getSamSeqDictFromGenomeFile(String genome) throws IOException {

		SAMSequenceDictionary samSeqDict= new SAMSequenceDictionary();

		if(Utils.getFileTypeFromName(genome).equals(TrackFormat.BAM)){
			// If genome is a bam file, get header from there.
			SamReaderFactory srf=SamReaderFactory.make();
			srf.validationStringency(ValidationStringency.SILENT);
			SamReader samReader= srf.open(new File(genome));
			samSeqDict= samReader.getFileHeader().getSequenceDictionary();
			if(!samSeqDict.isEmpty()){
				return samSeqDict;
			}
		}
		
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
				System.err.println("\nGenome file not found: " + genome + "\n");
				ex.printStackTrace();
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
		return samSeqDict;
	}

	public boolean equalCoords(GenomicCoords other){
		return this.chrom.equals(other.chrom) && this.from.equals(other.from) && this.to.equals(other.to); 
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
	 * */
	public TrackWiggles getGCProfile() throws IOException {
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
			// is not a mulitple of the sequence length, merge thi sround with the next in a 
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
		temp.delete();
		cgWiggle.setFileTag(this.gcProfileFileTag);
		cgWiggle.setyMaxLines(5);
		cgWiggle.setYLimitMin(0);
		cgWiggle.setYLimitMax(100);
		// cgWiggle.setYmin(Math.round(ymin * 10.0)/10.0);
		// cgWiggle.setYmax(Math.round(ymax * 10.0)/10.0); 
		// cgWiggle.update();
		return cgWiggle;
	}

	public void centerAndExtendGenomicCoords(GenomicCoords gc, int size, double slop) throws InvalidGenomicCoordsException {
		
		if(size <= 0){
			System.err.println("Invalid feature size. Must be > 0, got " + size);
			throw new InvalidGenomicCoordsException();
		}
		
		double center= (size/2.0) + gc.getFrom();
		gc.from= (int)Math.rint(center - (size * slop));
		gc.to= (int)Math.rint(center + (size * slop));
		
		if(((gc.to - gc.from)+1) < gc.windowSize){
			int span= (gc.to - gc.from);
			int extendBy= (int)Math.rint((gc.windowSize / 2.0) - (span / 2.0));
			gc.from -= extendBy;
			gc.to += extendBy;
		}
		gc.correctCoordsAgainstSeqDict(samSeqDict);
	}

	
	/* Getters and setters */
	
	public List<Double> getMapping() {
		return seqFromToLenOut();
	}	
	
	public String getChrom() {
		return chrom;
	}

	public Integer getFrom() {
		return from;
	}

	public Integer getTo() {
		return to;
	}
	
	/** Size of the window as set in input this.windowSize */
	public int getUserWindowSize(){
		return this.windowSize;
	}
	
	/** Size of genomic interval. Can be smaller than windowSize set by user. */
	public int getGenomicWindowSize(){
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

	public String getGcProfileFileTag(){
		return this.gcProfileFileTag;
	}

}