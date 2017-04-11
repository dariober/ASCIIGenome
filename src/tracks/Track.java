package tracks;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.UnsupportedEncodingException;
import java.net.MalformedURLException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.apache.commons.io.FileUtils;
import org.apache.commons.validator.routines.UrlValidator;

import com.google.common.base.Joiner;
import com.google.common.base.Splitter;

import coloring.Config;
import coloring.ConfigKey;
import coloring.Xterm256;
import exceptions.InvalidColourException;
import exceptions.InvalidCommandLineException;
import exceptions.InvalidGenomicCoordsException;
import exceptions.InvalidRecordException;
import htsjdk.samtools.filter.AlignedFilter;
import htsjdk.samtools.filter.SamRecordFilter;
import samTextViewer.GenomicCoords;
import samTextViewer.Main;
import samTextViewer.Utils;

public abstract class Track {

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
	protected List<Double> screenScores= new ArrayList<Double>();
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
	protected boolean rpm= false;
	protected static final int f_FLAG= 0; private int f_flag= f_FLAG;
	protected static final int F_FLAG= 4; private int F_flag= F_FLAG;
	protected static final int MAPQ= 0; private int mapq= MAPQ;
	protected List<SamRecordFilter> samRecordFilter= new ArrayList<SamRecordFilter>(); 
	private boolean hideTrack= false; 
	private boolean hideTitle= false;
	private TrackFormat trackFormat;
	protected String awk= ""; 
	private int printRawLineCount= -1; // Number of lines to print. Same as `head -n 10`
	
	/** A file to export track data
	 * */
	private String exportFile= null;
	private String cutScriptForPrinting= "";
//	private boolean appendToExportFile= false;
	
	/** Format the title string to add colour or return title as it is if
	 * no format is set.
	 * @throws InvalidColourException 
	 * */
	protected String formatTitle(String title) throws InvalidColourException{

		if(this.isNoFormat()){
			return title;
		} else {
			int colourCode= Config.get256Color(ConfigKey.title_colour);
			if(this.titleColour != null){
				colourCode= Xterm256.colorNameToXterm256(this.titleColour);
			}
			return "\033[48;5;" + Config.get256Color(ConfigKey.background) + ";38;5;" + colourCode + "m" + title;
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

	public void addBookmark(GenomicCoords gc, String nameForBookmark) throws IOException, ClassNotFoundException, InvalidRecordException, SQLException, InvalidGenomicCoordsException {
		throw new UnsupportedOperationException();		
	}

	public String getWorkFilename() {
		return workFilename;
	}

	public void setWorkFilename(String workFilename) {
		this.workFilename = workFilename;
	}

//	public abstract String printFeaturesToFile() throws IOException, InvalidGenomicCoordsException, InvalidColourException;

	/** Print raw lines after having processed them through `cut`, `clip etc.`.
	 * This method also add formatting so it returns a single string. 
	 * It should be used only for printing and not for any computation.
	 * If an output file has been set via this.setExportFile(exportFile), 
	 * write to file and return empty string. 
	 * @throws IOException 
	 * @throws InvalidGenomicCoordsException 
	 * @throws InvalidColourException 
	 * @throws InvalidCommandLineException 
	 * */
	public String printLines() throws InvalidGenomicCoordsException, IOException, InvalidColourException, InvalidCommandLineException{

		List<String> rawList= this.getRecordsAsStrings();
		
		if(this.getExportFile() != null && ! this.getExportFile().isEmpty()){
			// If an output file has been set, send output there and return. 
			BufferedWriter wr= null;
			try{
				wr = new BufferedWriter(new FileWriter(this.getExportFile(), true));
				for(String line : rawList){
					wr.write(line + "\n");
				}
				wr.close();
			} catch(IOException e){
				System.err.println("Cannot write to " + this.getExportFile());
				throw e;
			}
			return "";
		}
		
		int windowSize= this.getGc().getUserWindowSize();
		if(this.getPrintMode().equals(PrintRawLine.FULL)){
			windowSize= Integer.MAX_VALUE;
		} else if(this.getPrintMode().equals(PrintRawLine.CLIP)){
			// Keep windowSize as it is
		} else {
			return "";
		} 
		
		int count= this.getPrintRawLineCount();
		
		List<String> featureList= new ArrayList<String>();
		String omitString= "";
		for(String line : rawList){
			String parsedLine= this.cutLine(line);
			featureList.add(parsedLine);
			count--;
			if(count == 0){
				int omitted= rawList.size() - this.getPrintRawLineCount();
				if(omitted > 0){
					omitString= "[" + omitted + "/"  + rawList.size() + " features omitted]";
				}
				break;
			}
		}
		List<String> tabList= Utils.tabulateList(featureList);
		StringBuilder sb= new StringBuilder();
		if( ! omitString.isEmpty()){
			sb.append(omitString + "\n");
		}
		for(String x : tabList){
			if(x.length() > windowSize){
				x= x.substring(0, windowSize);
			}			
			sb.append(x + "\n");
		}
		
		if(this.isNoFormat()){
			return sb.toString();
		}

		String formatted=  "\033[38;5;" + Config.get256Color(ConfigKey.foreground) + 
		";48;5;" + Config.get256Color(ConfigKey.background) + "m" + sb.toString();
		
		return formatted; 
		
	}
	
	/** Interprets the content of this.cutScriptForPrinting to cut the input line
	 * in way similar to Unix cut. Args:
	 * -d <regex> Delimiter to split line, default '\t'
	 * -f <idx> Indexes to extract columns, no spaces allowed. E.g. 1,3,5-7
	 * @throws InvalidCommandLineException
	 * The parsing and interpretation of this.cutScriptForPrinting is done for a single line
	 * meaning that for many lines in input this is inefficient. However it should be ok for now,
	 * but consider passing List<String> as input so the interpretation is done once only. 
	 * */
	private String cutLine(String line) throws InvalidCommandLineException {
		if(this.cutScriptForPrinting == null || this.cutScriptForPrinting.isEmpty()){
			return line;
		}
		List<String> args= Utils.tokenize(this.cutScriptForPrinting, " ");
		String fields= Utils.getArgForParam(args, "-f");
		if(fields == null){
			//No fields given: Nothing to cut just return line as it is
			return line; 
		}
		String delim= Utils.getArgForParam(args, "-d");
		if(delim == null){
			delim= "\t"; 
		}
		// Split the line at delim
		List<String> lst= Arrays.asList(line.split(delim));
		List<Integer> idxs = (expandStringOfFieldsToIndexes(fields));
		List<String> outLst= new ArrayList<String>();
		for(int i : idxs){
			i--; // To make it zero-based
			if(i < 0){
				throw new InvalidCommandLineException();
			}
			if(i >= lst.size()){
				i= lst.size() - 1;
			}
			outLst.add(lst.get(i));
		}
		// We always re-join on tab, regardless of delimiter!
		return Joiner.on("\t").join(outLst);
	}

	/** Expand the string of fileds to return the individual indexes.
	 * E.g. "1-3,5-7,10-" -> [1,2,3,5,6,7,10, Integer.MAX_VALUE].
	 * Integer.MAX_VALUE signals that the fields from the last index to the end should be returned.    
	 * @throws InvalidCommandLineException 
	 * */
	private List<Integer> expandStringOfFieldsToIndexes(String fields) throws InvalidCommandLineException{
		List<Integer> idxs= new ArrayList<Integer>();
		
		// A bunch of checks for the validity of the input
		String checkOnlyDigits= fields.replaceAll(",", "").replaceAll("-", "");
		if( ( ! checkOnlyDigits.matches("[0-9]+")) || 
			  checkOnlyDigits.startsWith("0") ||
			  fields.contains("--") ||
			  fields.startsWith("-")){
			throw new InvalidCommandLineException();
		}
		if(fields.endsWith("-")){
			fields= fields.replaceAll("-$", "");
		}
		List<String> lst= Splitter.on(",").omitEmptyStrings().splitToList(fields);
		for(String x : lst){
			if(x.contains("-")){
				// This part expands the string "5-8" to [5,6,7,8]. "8-5" expanded to [8,7,6,5]
				String[] fromTo = x.split("-");
				int from= Integer.parseInt(fromTo[0]);
				int to= Integer.parseInt(fromTo[1]);
				if(from < to){
					for(int i= from; i <= to; i++){
						idxs.add(i);
					}
				} else {
					for(int i= from; i >= to; i--){
						idxs.add(i);
					}
				}
			} else {
				idxs.add(Integer.parseInt(x));
			}
		}
		return idxs;
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
		if(count < 0){
			count= Integer.MAX_VALUE;
		}
		this.printRawLineCount= count;
	}

	protected int getPrintRawLineCount() {
		return this.printRawLineCount;
	}

	/** Returns the records under the current genomic coordinates. As far as possible, 
	 * records are returned exactly as they appear in the raw input, e.g. raw vcf lines, raw sam lines 
	 * etc.
	 * */
	protected abstract List<String> getRecordsAsStrings();
	
	protected void setCutScriptForPrinting(String cutScriptForPrinting){
		this.cutScriptForPrinting= cutScriptForPrinting;
	}
	
	/**Return a single string where title and track have been concatenated.
	 * Concatenation is done in such way that "title" is not followed by newline if
	 * the topmost line of "track" has enough white spaces to accommodate the title.
	 * E.g.
	 * ```
	 * data.bed#1     >>>>
	 *       >>>>        >>>>>
	 * ```
	 * instead of
	 * ```
	 *  data.bed#1     
	 *                >>>>
	 *       >>>>        >>>>>
	 * ```
	 * @throws IOException 
	 * @throws InvalidGenomicCoordsException 
	 * @throws InvalidColourException 
	 * */
	public String concatTitleAndTrack() throws InvalidColourException, InvalidGenomicCoordsException, IOException{
		// * Strip ascii escapes
		// * Get length of leading whitespaces on topmost line of tracks
		// * if len(leadine whitespaces) > len(title):
		// * Remove the len(title) leading whitespaces from profile
		// * Return title + profile
		String track= this.printToScreen();
		String title= this.getTitle();
		int titleLen= Utils.stripAnsiCodes(title).trim().length();

		String sProfile= Utils.stripAnsiCodes(track);
		if(sProfile.trim().isEmpty()){ // No features in this profile
			return title.replaceAll("\n", "") + track; 	
		}
		int leadingSpaces= sProfile.indexOf(sProfile.trim());
		if(leadingSpaces > titleLen){
			while(titleLen > 0){
				track= track.replaceFirst(" ", "");
				titleLen--;
			}
			title= title.replaceAll("\n", "");
		}
		return title + track; 
	}
	
}

