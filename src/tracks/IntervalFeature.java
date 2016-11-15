package tracks;

import java.util.HashMap;
import java.util.List;

import org.apache.commons.lang3.math.NumberUtils;
import org.biojava.nbio.genome.parsers.gff.Feature;
import org.biojava.nbio.genome.parsers.gff.Location;

import samTextViewer.Utils;

import com.google.common.base.Splitter;
import com.google.common.collect.Lists;

import exceptions.InvalidGenomicCoordsException;

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

	private String raw; // Raw input string exactly as read from source.
	private TrackFormat format= TrackFormat.BED;
	/** Name to be displayed to the user */
	private String name= ".";
	/** Use this attribute to as key to assign the name field */
	private String gtfAttributeForName= null;
	
	/** Start position of feature in screen coordinates. 
	 * -1 if the feature is not part of the screenshot. */
	private int screenFrom= -1;
	private int screenTo= -1;
	
	/* C o n s t r u c t o r s */
		
	/**
	 * Create an IntervalFeature from a String. Typically this string is a bed or gtf line read from file.
	 * @param line
	 * @param type Format of this line
	 * @throws InvalidGenomicCoordsException 
	 */
	public IntervalFeature(String line, TrackFormat type) throws InvalidGenomicCoordsException{
		if(type.equals(TrackFormat.BED) || type.equals(TrackFormat.BEDGRAPH) || type.equals(TrackFormat.BIGBED)){
			this.intervalFeatureFromBedLine(line);
			this.format= TrackFormat.BED;
		} else if(type.equals(TrackFormat.GFF)){
			this.intervalFeatureFromGtfLine(line);
			this.format= TrackFormat.GFF;
		} else if(type.equals(TrackFormat.VCF)) {
			this.intervalFeatureFromVcfLine(line);
			this.format= TrackFormat.VCF;
		} else {
			System.err.println("Format " + type + " not supported");
			throw new RuntimeException();
		}
	}
	
	//public IntervalFeature(VariantContext variantContext){

	//	this.chrom= variantContext.getContig();
	//	this.from= variantContext.getStart();
	//	this.to= variantContext.getEnd();        
    //
	//	this.score= Float.NaN;
	//	
	//	this.raw= variantContext.toString();
	//	this.format= TrackFormat.VCF;
	//	
	//}
			
	/* M e t h o d s */
	
	private IntervalFeature intervalFeatureFromVcfLine(String vcfLine) throws InvalidGenomicCoordsException{
		vcfLine= vcfLine.replace("\n", "");
		vcfLine= vcfLine.replace("\r", "");
		this.raw= vcfLine;
		
		List<String> vcfList = Lists.newArrayList(Splitter.on("\t").split(vcfLine));
		if(vcfList.size() < 8){
			throw new RuntimeException("intervalFeatureFromVcfLine: Invalid vcf line:\n" + vcfList);
		}
				
		this.chrom= vcfList.get(0).trim();
		this.from= Integer.parseInt(vcfList.get(1)); // Make it 1-based

		// Feature coordinates are based on reference only. Insertions to the reference are indistinguishable from SNP (like IGV) 
		this.to= this.from + (vcfList.get(3).length()-1);
		this.name= vcfList.get(2);
		this.validateIntervalFeature();
		return this;
		
	}
	
	private IntervalFeature intervalFeatureFromBedLine (String bedLine) throws InvalidGenomicCoordsException{
		
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
			this.name= bedList.get(3);
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
	
	private IntervalFeature intervalFeatureFromGtfLine (String gtfLine) throws InvalidGenomicCoordsException{
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
		this.name= this.featureNameFromGTFAttribute(null);
		this.validateIntervalFeature();
		return this;
	}
	
	/**
	 * A bunch of checks to make sure feature is ok.
	 * @throws InvalidGenomicCoordsException 
	 * */
	private void validateIntervalFeature() throws InvalidGenomicCoordsException{

		if(!chrom.trim().equals(chrom)){
			System.err.println("Chrom name must not start or end with whitespaces. Got '" + chrom + "'");
			throw new InvalidGenomicCoordsException();
		}
		
		if(from < 1 || to < 1 || (from > to)){
			System.err.println("Invalid coordinates: " + from + " " + to);
			throw new InvalidGenomicCoordsException();
		}
		if(this.strand != '+' && this.strand != '-' && this.strand != '.'){
			System.err.println("Invalid strand char " + this.strand);
			throw new InvalidGenomicCoordsException();
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
			 throw new RuntimeException();
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
	
	/** Return true if x has the same coordinates of this object. Strand *not* taken into account */
	public boolean equals(IntervalFeature x){
		if(x == null){
			return false;
		}
		return (this.chrom.equals(x.chrom) && this.from == x.from && this.to == x.to);
	}

	/** Return true if x has the same coordinates of this object. Strand *is* taken into account */
	public boolean equalStranded(IntervalFeature x){
		if(x == null){
			return false;
		}
		return (this.chrom.equals(x.chrom) && this.from == x.from && this.to == x.to && this.strand == x.strand);
	}
	
	/** Returns the character to be used for this feature type and strand. E.g. 
	 * if the feature is exon on reverse strand return "e", if want it formatted return
	 *  "m31[e\033m" (or whatever formatting is needed). 
	 * */
	protected String assignTextToFeature(boolean noFormat) {
		
		if(this.format.equals(TrackFormat.VCF)){
			List<String> vcfList = Lists.newArrayList(Splitter.on("\t").split(this.raw));
			String text= FormatVCF.format(vcfList.get(3), vcfList.get(4), noFormat);
			return text;
		}
		
		// Get feature strand
		char strand= '.'; // Default for NA
		if(this.strand == '+'){
			strand= '+';
		} else if(this.strand == '-'){
			strand= '-';
		}
		// Get feature type
		String feature= this.getFeature().toLowerCase();
		HashMap<Character, HashMap<String, Character>> featureToTextCharDict = FormatGTF.getFeatureToTextCharDict();
		if(!featureToTextCharDict.get('.').containsKey(feature)){
			feature= "other"; // Feature type NA or not found in dict
		}
		
		// Now you have the right char to be used for this feature type and strand.
		char text= featureToTextCharDict.get(strand).get(feature);

		// Add formatting if required
		if(noFormat){
			return Character.toString(text);
		} else {
			return FormatGTF.format(text, strand);
		}
	}

	/** Get attribute value given key, e.g. transcript_id */
	protected String getAttribute(String attributeName){
		// * Get attribute field,
		if( ! this.format.equals(TrackFormat.GFF)){
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
		if(x != null ){
			x= x.trim();
		}
		return x; 
	}
	
	/**Return the name for the feature given the attributes.
	 * Here you decide what the user should see as name for different types of GTF/GFF attributes.
	 * @param attributeName: Attribute to use as key to get name from. E.g. gene_name. If null use 
	 * the default precedence rule to find one. If attributeName is not found, set name to '.' (missing).
	 * */
	private String featureNameFromGTFAttribute(String attributeName){

		String xname= this.name;
		if(attributeName != null){
			xname= this.getAttribute(attributeName);
			if(xname == null){
				xname= ".";
			}
		// Precedence to assign name 
		} else if (this.getAttribute("Name") != null){ xname= this.getAttribute("Name");
		
		} else if (this.getAttribute("ID") != null){ xname= this.getAttribute("ID");

		} else if (this.getAttribute("transcript_name") != null){ xname= this.getAttribute("transcript_name");
		
		} else if (this.getAttribute("transcript_id") != null){ xname= this.getAttribute("transcript_id");
		
		} else if (this.getAttribute("gene_name") != null){ xname= this.getAttribute("gene_name");
		
		} else if (this.getAttribute("gene_id") != null){ xname= this.getAttribute("gene_id");
		
		} else {
			// Leave as default
		}
		return xname;
	}
	
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

	/** Name be shown to the user */
	public String getName() {
		if(this.format.equals(TrackFormat.GFF)){
			return this.featureNameFromGTFAttribute(this.gtfAttributeForName);
		}
		return this.name;
	}

	public float getScore() {
		return score;
	}

	public char getStrand() {
		return strand;
	}

	public void setStrand(char strand) {
		this.strand= strand;
	}
	
	public int getScreenFrom() {
		return screenFrom;
	}
	public void setScreenFrom(int screenFrom) {
		this.screenFrom= screenFrom;
	}
	
	public int getScreenTo() {
		return screenTo;
	}
	public void setScreenTo(int screenTo) {
		this.screenTo= screenTo;
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

	public String getGtfAttributeForName() {
		return this.gtfAttributeForName;
	}

	public void setGtfAttributeForName(String gtfAttributeForName) {
		this.gtfAttributeForName = gtfAttributeForName;
	}
		
	@Override
	/**
	 * Sort by chrom, start, end, strand. 
	 */
	public int compareTo(IntervalFeature other) {
		
		int i= this.chrom.compareTo(other.chrom);
		    if (i != 0) return i;
		
		i = this.from - other.from;
		    if (i != 0) return i;

	    i = this.to - other.to;
		    if (i != 0) return i;

		i= Character.toString(this.strand).compareTo(Character.toString(other.strand));
		    if (i != 0) return i;
		    
		return i;
	}

//	public void setTo(int to) {
//		this.to= to;
//	}
//	public void setFrom(int from) {
//		this.from= from;
//	}
}
