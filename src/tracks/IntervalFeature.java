package tracks;

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

import org.apache.commons.lang3.math.NumberUtils;
import org.biojava.nbio.genome.parsers.gff.Feature;
import org.biojava.nbio.genome.parsers.gff.Location;

import com.google.common.base.Splitter;
import com.google.common.collect.Lists;

import exceptions.InvalidColourException;
import exceptions.InvalidGenomicCoordsException;
import samTextViewer.Utils;

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
	private TrackFormat trackFormat= TrackFormat.BED;
	/** Name to be displayed to the user */
	private String name= ".";
	/** Use this attribute to as key to assign the name field */
	private String gtfAttributeForName= null;
	
	/** Start position of feature in screen coordinates. 
	 * -1 if the feature is not part of the screenshot. */
	private int screenFrom= -1;
	private int screenTo= -1;
	final private String NAME_NA= "-na"; // String to use to set the name field to missing when retrieving feature name.
	
	// The feature as it would be represented on screen. Each String of the array
	// is a character to be printed on screen (e.g. "E") possibly formatted (e.g. "\033[m5;45E\033")
	private char[] ideogram;
	
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
			this.trackFormat= TrackFormat.BED;
		
		} else if(type.equals(TrackFormat.GFF) || type.equals(TrackFormat.GTF)){
			this.intervalFeatureFromGtfLine(line);
			this.trackFormat= TrackFormat.GFF;
					
		} else if(type.equals(TrackFormat.VCF)) {
			this.intervalFeatureFromVcfLine(line);
			this.trackFormat= TrackFormat.VCF;
		} else {
			System.err.println("Format " + type + " not supported");
			throw new RuntimeException();
		}
	}
	
	public IntervalFeature(String chrom, int from, int to, TrackFormat format) throws InvalidGenomicCoordsException{

		if(chrom == null || chrom.isEmpty()){
			throw new InvalidGenomicCoordsException();
		}
		if(from < 0 || to < 0){
			throw new InvalidGenomicCoordsException();
		}
		if(to < from){
			throw new InvalidGenomicCoordsException();
		}
		
		this.chrom= chrom;
		this.from= from;
		this.to= to;
		this.trackFormat= format;
	}
	
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

		// Feature coordinates are based on the longest operation, whether insertion or deletion
		// or snv (length 1). Note that representing insertions as strings of length > 1 is not
		// consistent with the genomic coordinates but it gives a better idea of the size of the 
		// insertion.
		int offset= vcfList.get(3).length() > vcfList.get(4).length() ? 
				vcfList.get(3).length() : 
				vcfList.get(4).length();
		
		this.to= this.from + (offset-1);
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
		this.name= this.getNameForIdeogram(null);
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
	
	/** Returns array of characters to be used for this feature type and strand.
	 * */
	private char[] makeIdeogram() {
		
		char text;
		if(this.trackFormat.equals(TrackFormat.VCF)){
			List<String> vcfList = Lists.newArrayList(Splitter.on("\t").split(this.raw));
			text= FormatVCF.textForVariant(vcfList.get(3), vcfList.get(4));

		} else {
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
			text= featureToTextCharDict.get(strand).get(feature);
		}
		char[] textArrayForScreen= new char[this.getScreenTo() - this.getScreenFrom() + 1];
		for(int i= 0; i < textArrayForScreen.length; i++){
			textArrayForScreen[i]= text;
		}
		return textArrayForScreen;
	}

	/** Returns a copy of input array with chars from feature name replaced */
	private char[] addNameToIdeogram(char[] xIdeoagram){
		
		if(this.getName() == null || this.getName().trim().isEmpty() || this.getName().trim().equals(".")){
			return xIdeoagram;
		}
		
		char[] ideogramWithName= Arrays.copyOf(xIdeoagram, xIdeoagram.length);
				
		// Find longest run of the same char to make space for the name
		int[] run= this.maxRun(ideogramWithName);
		
		// This is the amount of space (# chars) available to the name. Minus n because we leave some
		// space left and right.
		int space= run[1] - 1;
		if(space < 4){ // If the name has left than this many chars do not return it at all. 
			return ideogramWithName;
		}
						
		// Name to print, maybe shorter than the full name:
		String nameOnFeature= "_" + this.getName().trim() + "_";
		nameOnFeature= nameOnFeature.substring(0, Math.min(nameOnFeature.length(), space - 1));
		
		// This is where the name will start and end.
		int offset= (space - nameOnFeature.length())/2;
		int start= run[0] + offset + 1;
		int end= start + nameOnFeature.length();
		
		int j= 0;
		for(int i= (start); i < (end); i++){
			ideogramWithName[i]= nameOnFeature.charAt(j);
			j++;
		}
		ideogramWithName[end-1]= '_'; // Always terminate name with underscore 
		return ideogramWithName;
	}
	
	protected String[] makeIdeogramFormatted(boolean noFormat) throws InvalidColourException{

		char[] textArray;
		if(this.ideogram == null){
			// The array of chars has not been set from outside so create it:
			textArray= this.makeIdeogram();
		} else {
			textArray= this.ideogram;
		}
		
		textArray= this.addNameToIdeogram(textArray);
		
		String[] formattedTextForScreen= new String[textArray.length];
		for(int i=0 ; i < textArray.length; i++){
			char c= textArray[i];
			if(noFormat){
				formattedTextForScreen[i]= String.valueOf(c);
				
			} else if(this.trackFormat.equals(TrackFormat.GTF) 
					  || this.trackFormat.equals(TrackFormat.GFF)
					  || this.trackFormat.equals(TrackFormat.BED)){
				// You may want to set up formatting for each of these TrackFormats
				formattedTextForScreen[i]= FormatGTF.format(c, this.getStrand());
			
			} else if(this.trackFormat.equals(TrackFormat.VCF)){
				formattedTextForScreen[i]= FormatVCF.format(c);
			
			} else {
				System.out.println("I don't kknow how to format TrackFormat: " + this.trackFormat);
				throw new RuntimeException();
			}
		}
		return formattedTextForScreen;
	}
	
	/** Get attribute value given key, e.g. transcript_id */
	protected String getGFFValueFromKey(String attributeName){
		// * Get attribute field,
		if( ! (this.trackFormat.equals(TrackFormat.GFF) || this.trackFormat.equals(TrackFormat.GTF))){
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
	private String getNameForIdeogram(String attributeKey){

		String xname= this.name;
		
		if(attributeKey == null){

			// attribute name is null, assign default by Precedence to assign name 

			if (this.getGFFValueFromKey("Name") != null){ 
				xname= this.getGFFValueFromKey("Name");
			
			} else if (this.getGFFValueFromKey("ID") != null){ 
				xname= this.getGFFValueFromKey("ID");
	
			} else if (this.getGFFValueFromKey("transcript_name") != null){ 
				xname= this.getGFFValueFromKey("transcript_name");
			
			} else if (this.getGFFValueFromKey("transcript_id") != null){ 
				xname= this.getGFFValueFromKey("transcript_id");
			
			} else if (this.getGFFValueFromKey("gene_name") != null){ 
				xname= this.getGFFValueFromKey("gene_name");
			
			} else if (this.getGFFValueFromKey("gene_id") != null){ 
				xname= this.getGFFValueFromKey("gene_id");
			
			} else {
				// None of the above attributes found, leave as default
			}
			return xname;
		}
		
		if(attributeKey.equals(this.NAME_NA)){
			return ".";
		} 
		else if(! (this.trackFormat.equals(TrackFormat.GFF) || this.trackFormat.equals(TrackFormat.GTF))){
			// This is not a GFF/GTF file return default name without parsing attributes
			return xname;
		}				
		else {
			xname= this.getGFFValueFromKey(attributeKey);
			if(xname == null){
				xname= ".";
			}
			return xname;
		}
	}

	/**
	 * Find the start index and length of the longest run of identical characters. Examples
	 * AAAAAAzz   -> [0, 6]
	 * xAAAAAAzz  -> [1, 6]
	 * xxAAAAAAzz -> [2, 6]
	 * with ties the index of the first run found is returned.
	 * Modified from http://codereview.stackexchange.com/questions/75441/finding-the-length-of-the-largest-run
	 */
	private int[] maxRun(char[] seq) {
		
		if (seq == null || seq.length == 0) {
			int[] posMaxRun= {0, 0};
	    	return posMaxRun; 
	    }

		int maxStart= 0;
	    int maxRun = 1;
	    int currentRun = 1;

	    for (int i = 1; i < seq.length; i++) {
	        if (seq[i] == seq[i - 1]) {
	            currentRun++;
	            if (currentRun > maxRun) {
	            	maxRun = currentRun;
	            	maxStart= i - maxRun + 1;
	            }
	        } else {
	            currentRun = 1;
	        }
	    }
	    int[] posMaxRun= {maxStart, maxRun};
	    return posMaxRun;
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
		if(this.name != null &&  ! this.name.equals(".") && ! this.name.isEmpty() || this.raw == null){
			return this.name;
		}		
		return this.getNameForIdeogram(this.gtfAttributeForName);
	}

	public void setName(String name) {
		this.name= name;
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

	protected void setIdeogram(char[] ideogram) {

		if(ideogram == null){
			this.ideogram = null;
			return;
		}
		
		if((this.getScreenTo() - this.getScreenFrom() + 1) != ideogram.length){
			System.err.println("Length of text for screen (" + ideogram.length + ") "
					+ "does not equal feature length on screen from= " + this.getFrom() + " to= " + this.getScreenTo());
			throw new RuntimeException();
		}
		this.ideogram = ideogram;
	}
	
}
