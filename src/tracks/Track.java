package tracks;

import java.io.File;
import java.io.IOException;
import java.io.UnsupportedEncodingException;
import java.net.MalformedURLException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.io.FileUtils;
import org.apache.commons.validator.routines.UrlValidator;

import com.google.common.base.Joiner;

import coloring.Config;
import coloring.ConfigKey;
import coloring.Xterm256;
import exceptions.InvalidColourException;
import exceptions.InvalidGenomicCoordsException;
import exceptions.InvalidRecordException;
import htsjdk.samtools.filter.AlignedFilter;
import htsjdk.samtools.filter.SamRecordFilter;
import samTextViewer.GenomicCoords;
import samTextViewer.Main;
import samTextViewer.Utils;

// TODO: This class should be abstract 
public class Track {

	public static String awkFunc= "";
	
	static {
		  try {
		    awkFunc = FileUtils.readFileToString(new File(Main.class.getResource("/functions.awk").toURI()));
		  }
		  catch (Exception ex) {
		    /* Handle exception. */
		  }
		}
	
	private String title= "";
	protected int yMaxLines= 10;
	private String filename= "N/A"; // File name as given in input
	private String workFilename= "N/A"; // File actually used by ASCIIGenome. E.g. tmp tabix files 
	private String trackTag= "N/A"; // Tag name for title
	// private int id= 1;              // A unique identifier for the track. Changed when the track is added to a TrackSet. 
	private List<Double> screenScores= new ArrayList<Double>();
	private GenomicCoords gc;
	private boolean noFormat= false; 
	private double yLimitMin= Double.NaN; // Same as R ylim()
	private double yLimitMax= Double.NaN;
	/** Max size of genomic region before the track shuts down to prevent excessive slow down */
	protected final int MAX_REGION_SIZE= 100001;   
	
	protected String titleColour= null;
	protected boolean bisulf= false;

	private String gtfAttributeForName= null;
	/** Should features on with same coords be squashed into a single one? */
	private PrintRawLine printMode= PrintRawLine.OFF;
	private FeatureDisplayMode featureDisplayMode= FeatureDisplayMode.EXPANDED;
	private int gap= 1;
	private boolean rpm= false;
	protected static final int f_FLAG= 0; private int f_flag= f_FLAG;
	protected static final int F_FLAG= 4; private int F_flag= F_FLAG;
	protected static final int MAPQ= 0; private int mapq= MAPQ;
	protected List<SamRecordFilter> samRecordFilter= new ArrayList<SamRecordFilter>(); 
	private boolean hideTrack= false; 
	private boolean hideTitle= false;
	private TrackFormat trackFormat;
	protected String awk= ""; 
	
	/** A file to export track data
	 * */
	private String exportFile= null;
//	private boolean appendToExportFile= false;
	
	/** Format the title string to add colour or return title as it is if
	 * no format is set.
	 * @throws InvalidColourException 
	 * */
	protected String formatTitle(String title) throws InvalidColourException{

		if(this.isNoFormat()){
			return title;
		} else {
			int colourCode= Config.getColor(ConfigKey.title_colour);
			if(this.titleColour != null){
				colourCode= Xterm256.colorNameToXterm256(this.titleColour);
			}
			return "\033[48;5;" + Config.getColor(ConfigKey.background) + ";38;5;" + colourCode + "m" + title;
		}
	}
	
	/** Returns a string that parsed by `-exec` loads the current track with 
	 * the current settings (color, height, etc...).
	 * This method is currently quite approximative and it doesn't reproduce carefully all the settings.
	 * @throws IOException 
	 * @throws UnsupportedEncodingException 
	 * */
	public String settingsToString() throws UnsupportedEncodingException, IOException{
		String name= "^" + this.getTrackTag().replaceAll("(#|@)\\d+$", "");
		List<String> set= new ArrayList<String>();
		set.add("addTracks " + this.getFilename());
		set.add("colorTrack " + this.getTitleColour() + " " + name);
		set.add("trackHeight " + this.getyMaxLines() + " " + name);
		set.add("ylim " + this.getYLimitMin() + " " + this.getYLimitMax() + " " + name);
		set.add("samtools -q " + this.getMapq() + " -f " + this.get_f_flag() + " -F " + this.get_F_flag() + " " + name);
		set.add("grep -i " + this.getShowRegex() + " -e " + this.getHideRegex() + " " + name);
		if(this.isRpm()){
			set.add("rpm " + name);
		}
		if(this.isBisulf()){
			set.add("BSseq " + name);
		}
		if(this.isHideTitle()){
			set.add("hideTitle " + name);
		}
		return Joiner.on(" && ").join(set);
	}
	
	/* Printers */
	public String printToScreen() throws InvalidGenomicCoordsException, IOException, InvalidColourException{
		return null;
	}

//	public String printFeatures() throws InvalidGenomicCoordsException, IOException{
//		return "";
//	};
	
	/** Print track info - for debugging and development only.
	 * */
	public String toString(){
		return  "file name: " + this.getFilename() + 
				"; file type: " + Utils.getFileTypeFromName(this.getFilename()) +
				"; track tag: " + this.getTrackTag() +
				"; track class: " + this.getClass().getSimpleName();
	}
	
	/* Setters and getters */
	//public void setTitle(String title){
	//	this.title= title;
	//}
	public String getTitle() throws InvalidColourException, InvalidGenomicCoordsException, IOException{
		return this.title;
	}
	public int getyMaxLines() {
		return yMaxLines;
	}
	public void setyMaxLines(int yMaxLines) throws MalformedURLException, ClassNotFoundException, IOException, InvalidGenomicCoordsException, InvalidRecordException, SQLException {
		this.yMaxLines = yMaxLines;
		this.update();
	}
	public String getFilename() {
		return filename;
	}
	public void setFilename(String filename) {
		UrlValidator urlValidator = new UrlValidator();
		if(urlValidator.isValid(filename)){
			this.filename = filename;
		} else {
			this.filename = new File(filename).getAbsolutePath();
		}
	}
	
	public String getTrackTag() { 
		return trackTag; 
	}
	
	public void setTrackTag(String trackTag) { 
		this.trackTag = trackTag; 
	}
	
	protected List<Double> getScreenScores() {
		return screenScores;
	}
	protected void setScreenScores(List<Double> screenScores) {
		this.screenScores = screenScores;
	}

	public GenomicCoords getGc() {
		return gc;
	}

	/** Set the GenomicCoords object AND update the track by calling the update method.
	 * */
	public void setGc(GenomicCoords gc) throws MalformedURLException, ClassNotFoundException, IOException, InvalidGenomicCoordsException, InvalidRecordException, SQLException {
		this.gc = gc;
		this.update();
	}

	public boolean isNoFormat() { 
		return noFormat; 
	}
	
	public void setNoFormat(boolean noFormat) { 
		this.noFormat = noFormat; 
	}

	public Double getYLimitMin() { 
		return yLimitMin;
	}
	
	public void setYLimitMin(double ymin) { 
		this.yLimitMin = ymin;
	}

	public Double getYLimitMax() { 
		return yLimitMax; 
	}
	public void setYLimitMax(double ymax) { 
		this.yLimitMax = ymax; 
	}

	/** Return filter making sure the AlignedFilter to discard unmapped is set.
	 * */
	public List<SamRecordFilter> getSamRecordFilter() { 
		AlignedFilter unmapped = new AlignedFilter(true);
		if(!this.samRecordFilter.contains(unmapped)){
			this.samRecordFilter.add(unmapped); // Unmapped reads are always discarded	
		}
		return this.samRecordFilter; 
	}

	protected void setSamRecordFilter(List<SamRecordFilter> samRecordFilter) throws MalformedURLException, ClassNotFoundException, IOException, InvalidGenomicCoordsException, InvalidRecordException, SQLException {
		this.samRecordFilter = samRecordFilter;
		this.update();
	}
	
	public boolean isBisulf() { return this.bisulf; }
	public void setBisulf(boolean bisulf) { this.bisulf= bisulf; }

	public void setAwk(String awk) throws ClassNotFoundException, IOException, InvalidGenomicCoordsException, InvalidRecordException, SQLException {

	}
	public String getAwk(){
		return "";
	}
	
	public void setHideRegex(String hideRegex) throws ClassNotFoundException, IOException, InvalidGenomicCoordsException, InvalidRecordException, SQLException { 
	
	}
	
	public String getHideRegex() { 
		return ""; 
	}
	
	public void setShowRegex(String showRegex) throws ClassNotFoundException, IOException, InvalidGenomicCoordsException, InvalidRecordException, SQLException { 

	}
	
	public String getShowRegex() { 
		return ""; 
	}
	
	public String getTitleColour() {
		if(this.titleColour == null){
			return Config.get(ConfigKey.title_colour);
		}
		return this.titleColour;
	}

	public void setTitleColour(String colour) {
//		if(!Utils.xterm256ColorCodes().containsKey(colour)){
//			try {
//				throw new InvalidColourException();
//			} catch (InvalidColourException e) {
//				// e.printStackTrace();
//				System.err.println("\nGot invalid colour: " + colour + ". Resetting to default");
//				System.err.println("Valid colours are: " + Utils.xterm256ColorCodes().keySet());
//				colour= "blue";
//			} 
//		}
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

	public void setPrintMode(PrintRawLine printMode) {
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

	public void setRpm(boolean rpm) throws MalformedURLException, ClassNotFoundException, IOException, InvalidGenomicCoordsException, InvalidRecordException, SQLException {
		this.rpm = rpm;
		this.updateToRPM();
	}
	public boolean isRpm(){
		return this.rpm;
	}

	/** Update scores to RPM. Do nothing if RPM transformation is not applicable.*/ 
	protected void updateToRPM(){
		// 
	}
	
	/** This int is just a setting but is NOT translated to a filter! */
	protected int get_f_flag() {
		return f_flag;
	}

	/** This int is just a setting but is NOT translated to a filter! */
	protected void set_f_flag(int f_flag) {
		this.f_flag = f_flag;
	}

	/** You should use a converter to get int from list of filters. 
	 * */
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
	
	public String getPrintableConsensusSequence() throws IOException, InvalidGenomicCoordsException, InvalidColourException{
		return "";
	}

	protected void update() throws MalformedURLException, IOException, InvalidGenomicCoordsException, InvalidRecordException, ClassNotFoundException, SQLException {
		// TODO Auto-generated method stub
	}

	public String getSeqRegex() {
		return null;
	}

	public void setSeqRegex(String seqRegex) throws ClassNotFoundException, IOException, InvalidGenomicCoordsException, InvalidRecordException, SQLException {
		//
	}

	public boolean isHideTrack() {
		return hideTrack;
	}

	protected void setHideTrack(boolean hideTrack) {
		this.hideTrack = hideTrack;
	}

	public boolean isHideTitle() {
		return hideTitle;
	}

	public void setHideTitle(boolean hideTitle) {
		this.hideTitle = hideTitle;
	}

	public void addBookmark(String nameForBookmark) throws IOException, ClassNotFoundException, InvalidRecordException, SQLException, InvalidGenomicCoordsException {
		// TODO Auto-generated method stub
		
	}

	public String getWorkFilename() {
		return workFilename;
	}

	public void setWorkFilename(String workFilename) {
		this.workFilename = workFilename;
	}

	public String printFeaturesToFile() throws IOException, InvalidGenomicCoordsException, InvalidColourException {
		// TODO Auto-generated method stub
		return "";
	}

	/**Return the export file name. The variable %r is expanded to coordinates. 
	 * */
	public String getExportFile() {
		if(exportFile != null){
			String x= exportFile.replaceAll("%r", this.getGc().getChrom() + "_" + this.getGc().getFrom() + "_" + this.getGc().getTo()); 
			return x;
		}
		return exportFile;
	}

	public void setExportFile(String exportFile) {
		this.exportFile = exportFile;
	}

//	public boolean isAppendToExportFile() {
//		return appendToExportFile;
//	}
//
//	public void setAppendToExportFile(boolean appendToExportFile) {
//		this.appendToExportFile = appendToExportFile;
//	}

	protected TrackFormat getTrackFormat() {
		return trackFormat;
	}

	protected void setTrackFormat(TrackFormat trackFormat) {
		this.trackFormat = trackFormat;
	}

	public List<String> getChromosomeNames() {
		throw new RuntimeException("TO BE IMPLEMENTED");	
	}

	protected void setPrintRawLineCount(int count) {

	}
	
}

