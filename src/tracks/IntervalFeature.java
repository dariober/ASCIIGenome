package tracks;

import java.util.HashMap;
import java.util.List;

import org.apache.commons.lang3.math.NumberUtils;
import org.apache.commons.lang3.text.StrTokenizer;
import org.biojava.nbio.genome.parsers.gff.Feature;
import org.biojava.nbio.genome.parsers.gff.Location;

import samTextViewer.Utils;

import com.google.common.base.Splitter;
import com.google.common.collect.Lists;

/**
 * Class to hold bed or gtf features. Behaviour should be similar to pybedtools Interval.
 * Feature coords are 1-based. The first ten bases of the chrom have from-to = 1-10
 * @author berald01
 *
 */
public class IntervalFeature implements Comparable<IntervalFeature>{

	// When reading bed files, we expect fields to be in this order.
	private String chrom;       // Required
	private int from;           // Required. NB 1 based also for bed files.
	private int to;             // Required 

	private float score= Float.NaN;
	private char strand= '.'; 
	private String source= "."; // Gtf specific
	private String feature= "."; // Gtf specific

	private String raw; // Raw input string exactly as read from file.
	private TrackFormat format= TrackFormat.BED;
	
	/** Start position of feature in screen coordinates. 
	 * -1 if the feature is not part of the screenshot. */
	private int screenFrom= -1;
	private int screenTo= -1;
	
	/* C o n s t r u c t o r s */
		
	/**
	 * Create an IntervalFeature from a String. Typically this string is a bed or gtf line read from file.
	 * @param line
	 * @param type Format of this line
	 */
	public IntervalFeature(String line, TrackFormat type){
		if(type.equals(TrackFormat.BED) || type.equals(TrackFormat.BEDGRAPH)){
			intervalFeatureFromBedLine(line);
			this.format= TrackFormat.BED;
		} else if(type.equals(TrackFormat.GFF)){
			intervalFeatureFromGtfLine(line);
			this.format= TrackFormat.GFF;
		} else {
			System.err.println("Format " + type + " not supported");
			throw new RuntimeException();
		}
	}
	
	//public IntervalFeature(IntervalFeature intervalFeature){
	//	
	//}
	
	/* M e t h o d s */
	
	private IntervalFeature intervalFeatureFromBedLine (String bedLine){
		
		bedLine= bedLine.replace("\n", "");
		bedLine= bedLine.replace("\r", "");
		this.raw= bedLine;
		
		List<String> bedList = Lists.newArrayList(Splitter.on("\t").split(bedLine));
		if(bedList.size() < 3){
			throw new RuntimeException("intervalFeatureFromBedLine: Invalid bed line:\n" + bedList);
		}
		this.chrom= bedList.get(0).trim();
		this.from= Integer.parseInt(bedList.get(1)) + 1; // Make it 1-based
		this.to= Integer.parseInt(bedList.get(2));

		// Process optional fields
		if(bedList.size() > 3){
			this.source= bedList.get(3);
		}
		if(bedList.size() > 4){
			if(NumberUtils.isNumber(bedList.get(4))){ // NB: Returns false if leading or trailing spaces are present.
				this.score= Float.valueOf(bedList.get(4));
			}
		}
		if(bedList.size() > 5){
			if(bedList.get(5).equals("+")){
				this.strand= '+';
			} else if(bedList.get(5).equals("-")){
				this.strand= '-';
			} else {
				this.strand= '.';
			}
		}
		this.validateIntervalFeature();
		return this;
	}
	
	private IntervalFeature intervalFeatureFromGtfLine (String gtfLine){
		//chr1    unknown exon    11874   12227   .       +       .       gene_id "DDX11L1"; transcript_id "NR_046018_1"; gene_name "DDX11L1"; tss_id "TSS14523";
		gtfLine= gtfLine.replace("\n", "");
		gtfLine= gtfLine.replace("\r", "");
		this.raw= gtfLine;

		List<String> gtfList = Lists.newArrayList(Splitter.on("\t").split(gtfLine));
		this.chrom= gtfList.get(0).trim();
		this.source= gtfList.get(1).trim();
		this.feature= gtfList.get(2).trim();
		this.from= Integer.parseInt(gtfList.get(3));
		this.to= Integer.parseInt(gtfList.get(4));
		try{
			this.score= Float.parseFloat(gtfList.get(5));
		} catch (NumberFormatException e){
			this.score= Float.NaN;
		}
		
		// Strand
		if(gtfList.get(6).trim().length() == 0){
			this.strand= '.';
		} else {
			char strand= gtfList.get(6).trim().charAt(0);	
			if(strand == '+' || strand == '-'){
				this.strand= strand;
			} else {
				this.strand= '.';
			}
		}
		//if(gtfList.get(7).trim().length() == 0){
		//	this.frame= '.';
		//} else {
		//	this.frame= gtfList.get(7).trim().charAt(0);
		//}
		//this.attribute= gtfList.get(8).trim();
		this.validateIntervalFeature();
		return this;
	}
	
	/**
	 * A bunch of checks to make sure feature is ok.
	 * */
	private void validateIntervalFeature(){

		if(!chrom.trim().equals(chrom)){
			System.err.println("Chrom name must not start or end with whitespaces. Got '" + chrom + "'");
			System.exit(1);
		}
		
		if(from < 1 || to < 1 || (from > to)){
			System.err.println("Invalid coordinates: " + from + " " + to);
			System.exit(1);
		}
		if(this.strand != '+' && this.strand != '-' && this.strand != '.'){
			System.err.println("Invalid strand char " + this.strand);
			System.exit(1);			
		}		
	}
	 
	/** 
	 * Map interval to screen coordinates using the provided ruler.
	 * @param rulerMap List typically obtained from Ruler_TO_BE_DEPRECTED.mapping of length equal to the screen 
	 * width mapping genome coords to screen coords.
	 * */
	public void mapToScreen(List<Double> rulerMap) {

		/*        |============| <- ruler
		 *   ===                  ===  <- Interval(s) 
		 */	
		if((this.from < rulerMap.get(0) && this.to < rulerMap.get(0)) ||
				(this.from > rulerMap.get(rulerMap.size()-1)) && this.to > rulerMap.get(rulerMap.size()-1)){
			this.screenFrom= -1;
			this.screenTo= -1;
			return;
		}
		
		/*
		 * Feature fully contains ruler map?
		 *        |============| <- ruler
		 *   ===================== <- Interval 
		 */
		if(this.from <= rulerMap.get(0) && this.to >= rulerMap.get(rulerMap.size()-1)){
			this.screenFrom= 0;
			this.screenTo= rulerMap.size()-1;
			return;
		}
		
		// Feature is all or partially contained
		screenFrom= Utils.getIndexOfclosestValue(this.from, rulerMap);
		screenTo= Utils.getIndexOfclosestValue(this.to, rulerMap);
		/*        |============|      <- ruler
		 *   ========   ===    =====  <- Interval(s) 
		 */	
		 if(screenFrom == -1){
			 screenFrom= 0;
		 }
		 if(screenTo == -1){
			 screenTo= rulerMap.size()-1;
		 }
		 if(screenFrom == -1 || screenTo == -1){
			 System.err.println("Unexpected mapping of features to ruler.");
			 System.exit(1);
		 }
	}	

//	public String toGtfString(){
//		String scoreStr= (Float.isNaN(this.score)) ? "." : String.valueOf(this.score);
//		String feature= this.chrom + "\t" + 
//					    this.source + "\t" +
//					    this.feature + "\t" +
//					    this.from + "\t" +
//					    this.to + "\t" +
//					    scoreStr + "\t" +
//					    this.strand + "\t" +
//					    this.frame + "\t" +
//					    this.attribute;
//		return feature;
//	}
	
	/* For debugging only */
	public String toString(){
		String feature= this.chrom + ":" + this.from + "-" + this.to + ", " 
				+ this.source + ", " 
				+ this.score + ", " 
				+ this.strand;
		feature += "\nScreen coords from, to: " + this.screenFrom + ", " + this.screenTo;
		return feature;
	}
	
	/** Return true x has the same coordinates of this object */
	public boolean equalCoords(IntervalFeature x){
		return (this.chrom.equals(x.chrom) && this.from == x.from && this.to == x.to);
	}
	
	protected String getAttribute(String attributeName){
		// * Get attribute field,
		if(this.format != TrackFormat.GFF){
			return null;
		}
		String[] line= this.raw.split("\t");
		if(line.length < 9){
			return null;
		}
		Location location= new Location(Integer.parseInt(line[3]), Integer.parseInt(line[4]));
		double score= line[5].equals(".") ? Double.NaN : Double.parseDouble(line[5]);
		int frame= line[7].equals(".") ? -1 : Integer.parseInt(line[7]);
		Feature gff= new Feature(line[0], line[1], line[2], location, score, frame, line[8]);
		String x= gff.getAttribute(attributeName);
		
		// String attributes= this.raw.split("\t")[8];
		// StrTokenizer str= new StrTokenizer();
		// str.setQuoteChar('\'');
		// List<String> tokens= str.getTokenList();
		return x; 
	}
	
	//protected Object clone() {
	//	try {
	//		return super.clone();
	//	} catch (CloneNotSupportedException e) {
	//		return null;
	//	}
	//}
	
	/*   S e t t e r s   and   G e t t e r s   */
	
	public String getChrom() {
		return chrom;
	}

	public int getFrom() {
		return from;
	}

	public int getTo() {
		return to;
	}

	public String getName() {
		return source;
	}

	public float getScore() {
		return score;
	}

	public char getStrand() {
		return strand;
	}

	public int getScreenFrom() {
		return screenFrom;
	}

	public int getScreenTo() {
		return screenTo;
	}
	
	public String getSource() {
		return source;
	}

	public void setSource(String source) {
		this.source = source;
	}

	public String getFeature() {
		return feature;
	}

	public void setFeature(String feature) {
		this.feature = feature;
	}

	public String getRaw() {
		return raw;
	}
	
	@Override
	/**
	 * Sort by chrom, start, end. Should be the same as Unix `sort -k1,1 -k2,2n -k3,3n` 
	 */
	public int compareTo(IntervalFeature other) {
		
		int i= this.chrom.compareTo(other.chrom);
		    if (i != 0) return i;
		
		i = this.from - other.from;
		    if (i != 0) return i;

	    i = this.to - other.to;
		    if (i != 0) return i;

		return i;
	}

	protected String assignTextToFeature(boolean noFormat) {

		/* Map GTF features to characters. Forward capital LETTERS, reverse small letters  
		 * Feature names are case insensitive */
		HashMap<String, Character> fwdFeature= new HashMap<String, Character>();
		HashMap<String, Character> revFeature= new HashMap<String, Character>();
		HashMap<String, Character> unstrFeature= new HashMap<String, Character>();
		fwdFeature.put("exon",        'E'); revFeature.put("exon",        'e');
		fwdFeature.put("cds", 	      'C'); revFeature.put("cds",         'c');
		fwdFeature.put("start_codon", 'A'); revFeature.put("start_codon", 'a');
		fwdFeature.put("stop_codon",  'Z'); revFeature.put("stop_codon",  'z');
		fwdFeature.put("utr",         'U'); revFeature.put("utr",         'u');
		fwdFeature.put("3utr",        'U'); revFeature.put("3utr",        'u');
		fwdFeature.put("5utr",        'W'); revFeature.put("3utr",        'w');
		fwdFeature.put("gene",		  'G'); revFeature.put("gene",        'g');
		fwdFeature.put("transcript",  'T'); revFeature.put("transcript",  't');
		fwdFeature.put("mrna",        'M'); revFeature.put("mrna",        'm');
		fwdFeature.put("trna",        'X'); revFeature.put("trna",        'x');
		fwdFeature.put("rrna", 		  'R'); revFeature.put("rrna",        'r');
		fwdFeature.put("mirna",       'I'); revFeature.put("mirna",       'i');
		fwdFeature.put("ncrna",       'L'); revFeature.put("ncrna",       'l');
		fwdFeature.put("lncrna",      'L'); revFeature.put("lncrna",      'l');
		fwdFeature.put("sirna",       'S'); revFeature.put("sirna",       's');
		fwdFeature.put("pirna",       'P'); revFeature.put("pirna",       'p');
		fwdFeature.put("snorna",      'O'); revFeature.put("snorna",      'o');
		
		// For feature with strand not available, use forward encoding, unless feature unknown
		unstrFeature.putAll(fwdFeature);
		fwdFeature.put("other", '>'); 		
		revFeature.put("other", '<');
		unstrFeature.put("other", '|');
		
		HashMap<Character, HashMap<String, Character>> featureToTextCharDict= 
					new HashMap<Character, HashMap<String, Character>>();
		featureToTextCharDict.put('+', fwdFeature);
		featureToTextCharDict.put('-', revFeature);
		featureToTextCharDict.put('.', unstrFeature);

		// Get feature strand
		char strand= '.'; // Default for NA
		if(this.strand == '+'){
			strand= '+';
		} else if(this.strand == '-'){
			strand= '-';
		}
		// Get feature type
		String feature= this.getFeature().toLowerCase();
		if(!featureToTextCharDict.get('.').containsKey(feature)){
			feature= "other"; // Feature type NA or not found in dict
		}

		// Now you have the right char to be used for this feature type and strand.
		char text= featureToTextCharDict.get(strand).get(feature);

		// Add formatting if required
		if(noFormat){
			return Character.toString(text);
		} else if(strand == '+') {
			return "\033[48;5;147;38;5;240m" + text + "\033[0m";
		} else if(strand == '-') {
			return "\033[48;5;225;38;5;240m" + text + "\033[0m";
		} else {
			return "\033[48;5;250;38;5;240m" + text + "\033[0m";
		}			
	}
}
