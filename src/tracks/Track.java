package tracks;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import exceptions.InvalidColourException;
import exceptions.InvalidGenomicCoordsException;
import htsjdk.samtools.filter.AlignedFilter;
import htsjdk.samtools.filter.SamRecordFilter;
import samTextViewer.GenomicCoords;
import samTextViewer.Utils;

public class Track {

	private String title= "";
	protected int yMaxLines= 10;
	private String filename= "N/A"; // File name as given in input
	private String fileTag= "N/A"; // File name for title
	private List<Double> screenScores= new ArrayList<Double>();
	private GenomicCoords gc;
	private boolean noFormat= false; 
	private double yLimitMin= Double.NaN; // Same as R ylim()
	private double yLimitMax= Double.NaN;
	/** Max size of genomic region before the track shuts down to prevent excessive slow down */
	protected final int MAX_REGION_SIZE= 100000;   
	
	protected String titleColour= "default";
	protected boolean bisulf= false;

	private String gtfAttributeForName= null;
	/** Should features on with same coords be squashed into a single one? */
	// protected boolean squash= false;
	private PrintRawLine printMode= PrintRawLine.OFF;
	private FeatureDisplayMode featureDisplayMode= FeatureDisplayMode.EXPANDED;
	private int gap= 1;
	private boolean rpm= false;
	private int f_flag= 0;
	private int F_flag= 4;
	private int mapq= 0;
	// private boolean printPileup= false;
	private List<SamRecordFilter> samRecordFilter= new ArrayList<SamRecordFilter>(); 
	
	/* Min value of screen scores. Not to be confused with the y limit **/
	public double getMinScreenScores(){
		Double ymin= Double.NaN;
	 	for(Double x : this.screenScores){
	 		if(ymin.isNaN() && !x.isNaN()){
	 			ymin= x;
	 		} else if (x < ymin){
	 			ymin= x;
	 		}
	 	}
	 	return ymin;
	}
	
	/* Max value of screen scores. Not to be confused with the y limit **/
	public double getMaxScreenScores(){
		Double ymax= Double.NaN;
	 	for(Double x : this.screenScores){
	 		if(ymax.isNaN() && !x.isNaN()){
	 			ymax= x;
	 		} else if (x > ymax){
	 			ymax= x;
	 		}
	 	}
	 	return ymax;
	}
	
	/** Format the title string to add colour or return title as it is if
	 * no format is set.
	 * */
	protected String formatTitle(String title){
		if(this.isNoFormat()){
			return title;
		} else {
			int colourCode= Utils.ansiColorCodes().get(this.titleColour);
			return "\033[0;" + colourCode + "m" + title + "\033[0m";
		}
	}
	
	/* Printers */
	public String printToScreen() throws InvalidGenomicCoordsException{
		return null;
	}

	public String printFeatures(int windowSize){
		return "";
	};
	
	public String toString(){
		return this.getFilename();
	}
	
	/* Setters and getters */
	//public void setTitle(String title){
	//	this.title= title;
	//}
	public String getTitle(){
		return this.title;
	}
	public int getyMaxLines() {
		return yMaxLines;
	}
	public void setyMaxLines(int yMaxLines) {
		this.yMaxLines = yMaxLines;
	}
	public String getFilename() {
		return filename;
	}
	public void setFilename(String filename) {
		this.filename = filename;
	}

	public String getFileTag() { return fileTag; }
	
	public void setFileTag(String fileTag) { this.fileTag = fileTag; }
	
	protected List<Double> getScreenScores() {
		return screenScores;
	}
	protected void setScreenScores(List<Double> screenScores) {
		this.screenScores = screenScores;
	}

	public GenomicCoords getGc() {
		return gc;
	}

	public void setGc(GenomicCoords gc) {
		this.gc = gc;
	}

	public boolean isNoFormat() { return noFormat; }
	public void setNoFormat(boolean noFormat) { this.noFormat = noFormat; }

	public double getYLimitMin() { return yLimitMin;}
	public void setYLimitMin(double ymin) { this.yLimitMin = ymin;}

	public double getYLimitMax() { return yLimitMax; }
	public void setYLimitMax(double ymax) { this.yLimitMax = ymax; }

	/** Return filter making sure the AlignedFilter to discard unmapped is set.
	 * */
	public List<SamRecordFilter> getSamRecordFilter() { 
		AlignedFilter unmapped = new AlignedFilter(true);
		if(!this.samRecordFilter.contains(unmapped)){
			this.samRecordFilter.add(unmapped); // Unmapped reads are always discarded	
		}
		return this.samRecordFilter; 
	}

	protected void setSamRecordFilter(List<SamRecordFilter> samRecordFilter) {
		this.samRecordFilter = samRecordFilter;
	}
	
	public boolean isBisulf() { return this.bisulf; }
	public void setBisulf(boolean bisulf) { this.bisulf= bisulf; }

	public void setHideRegex(String hideRegex) { }
	public String getHideRegex() { return ""; }
	
	public void setShowRegex(String showRegex) { }
	public String getShowRegex() { return ""; }
	
	public String getTitleColour() {
		return this.titleColour;
	}

	public void setTitleColour(String colour) {
		if(!Utils.ansiColorCodes().containsKey(colour)){
			try {
				throw new InvalidColourException();
			} catch (InvalidColourException e) {
				// e.printStackTrace();
				System.err.println("\nGot invalid colour: " + colour + ". Resetting to default");
				System.err.println("Valid colours are: " + Utils.ansiColorCodes().keySet());
				colour= "blue";
			} 
		}
		this.titleColour = colour;
	}
	
	public String getGtfAttributeForName() {
		return this.gtfAttributeForName;
	}

	public void setGtfAttributeForName(String gtfAttributeForName) {
		this.gtfAttributeForName = gtfAttributeForName;
	}

	public PrintRawLine getPrintMode() {
		return printMode;
	}

	protected void setPrintMode(PrintRawLine printMode) {
		this.printMode = printMode;
	}

	public FeatureDisplayMode getFeatureDisplayMode() {
		return featureDisplayMode;
	}

	public void setFeatureDisplayMode(FeatureDisplayMode featureDisplayMode) {
		this.featureDisplayMode = featureDisplayMode;
	}

	protected int getGap() {
		return gap;
	}

	protected void setGap(int gap) {
		if(gap < 0){
			throw new RuntimeException("Cannot set gap < 0");
		}
		this.gap = gap;
	}

	public void setRpm(boolean rpm) {
		this.rpm = rpm;
	}
	public boolean isRpm(){
		return this.rpm;
	}

	//public void setPrintPileup(boolean printPileup) {
	//	this.printPileup = printPileup;
	//}
	//public boolean isPrintPileup(){
	//	return this.printPileup;
	//}

	
	/** This int is just a setting but is NOT translated to a filter! */
	protected int get_f_flag() {
		return f_flag;
	}

	/** This int is just a setting but is NOT translated to a filter! */
	protected void set_f_flag(int f_flag) {
		this.f_flag = f_flag;
	}

	/** This int is just a setting but is NOT translated to a filter! */
	protected int get_F_flag() {
		return F_flag;
	}

	/** This int is just a setting but is NOT translated to a filter! */
	protected void set_F_flag(int F_flag) {
		this.F_flag = F_flag;
	}

	/** This int is just a setting but is NOT translated to a filter! */
	public int getMapq() {
		return mapq;
	}

	/** This int is just a setting but is NOT translated to a filter! */
	protected void setMapq(int mapq) {
		this.mapq = mapq;
	}
	
	public List<String> printPileupList(){
		return new ArrayList<String>();
	}
	
	public String getPrintableConsensusSequence() throws IOException{
		return "";
	}
	
}

