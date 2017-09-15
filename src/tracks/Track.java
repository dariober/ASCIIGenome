package tracks;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.io.UnsupportedEncodingException;
import java.net.HttpURLConnection;
import java.net.MalformedURLException;
import java.net.URL;
import java.nio.file.Files;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.regex.Pattern;

import org.apache.commons.io.IOUtils;
import org.apache.commons.validator.routines.UrlValidator;

import com.google.common.base.Joiner;

import coloring.Config;
import coloring.ConfigKey;
import coloring.Xterm256;
import exceptions.InvalidColourException;
import exceptions.InvalidCommandLineException;
import exceptions.InvalidGenomicCoordsException;
import exceptions.InvalidRecordException;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.filter.AggregateFilter;
import htsjdk.samtools.filter.AlignedFilter;
import htsjdk.samtools.filter.SamRecordFilter;
import samTextViewer.GenomicCoords;
import samTextViewer.Main;
import samTextViewer.Utils;

public abstract class Track {

	public static String awkFunc= "";
	
	static {
		  try {
			  InputStream in= Main.class.getResourceAsStream("/functions.awk");
			  awkFunc= IOUtils.toString(in);
			  in.close();
		  }
		  catch (Exception ex) {
		    /* Handle exception. */
		  }
		}
	
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
	protected final int MAX_REGION_SIZE= 1000001;   
	
	protected String titleColour= null;
	protected boolean bisulf= false;

	private String gtfAttributeForName= null;
	/** Should features on with same coords be squashed into a single one? */
	private PrintRawLine printMode= PrintRawLine.OFF;
	private FeatureDisplayMode featureDisplayMode= FeatureDisplayMode.EXPANDED;
	private int gap= 1;
	protected boolean readsAsPairs= false;
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
	private GenotypeMatrix genotypeMatrix= new GenotypeMatrix();
	/** A file to export track data
	 * */
	private String exportFile= null;
	private String systemCommandForPrint;
	private boolean printNormalizedVcf= false;
	private long lastModified;
	final static public String HIDE_REGEX= "^$"; protected String hideRegex= HIDE_REGEX;
	final static public String SHOW_REGEX= ".*"; protected String showRegex= SHOW_REGEX;

	
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
				new Xterm256();
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
	public abstract String getTitle() throws InvalidColourException, InvalidGenomicCoordsException, IOException;
	
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

	public abstract void setAwk(String awk) throws ClassNotFoundException, IOException, InvalidGenomicCoordsException, InvalidRecordException, SQLException;
	
	public abstract String getAwk();
	
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

	public abstract void update() throws MalformedURLException, IOException, InvalidGenomicCoordsException, InvalidRecordException, ClassNotFoundException, SQLException;

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

		List<String> rawList= this.execSystemCommand(this.getRecordsAsStrings(), this.getSystemCommandForPrint());
		
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
			featureList.add(line);
			count--;
			if(count == 0){
				int omitted= rawList.size() - this.getPrintRawLineCount();
				if(omitted > 0){
					omitString= "[" + omitted + "/"  + rawList.size() + " features omitted]";
				}
				break;
			}
		}
		List<String> tabList= Utils.tabulateList(featureList, this.getGc().getUserWindowSize());
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

	public String getSystemCommandForPrint() {
		return this.systemCommandForPrint;
	}

	/**Stream the list of string recordsAsStrings through the system command(s) given in 
	 * sysCmd. The system command(s) must read from stdin and write to stdout.
	 * @throws IOException 
	 * @throws InvalidCommandLineException 
	 * @throws InterruptedException 
	 * */
	private List<String> execSystemCommand(List<String> recordsAsStrings, String sysCmd) throws IOException {
		if(sysCmd == null || sysCmd.isEmpty()){
			return recordsAsStrings;
		}
		File tmp= Files.createTempFile("asciigenome.", ".print.tmp").toFile();
		tmp.deleteOnExit();
		BufferedWriter writer = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(tmp.getAbsolutePath())));
        for(String line : recordsAsStrings){
        	writer.append(line.replaceAll("\n$", ""));
        	writer.append("\n");
        }
		writer.close();
		
		ArrayList<String> cmd= new ArrayList<String>();
		cmd.add("bash");
		cmd.add("-c");
		cmd.add("cat " + tmp.getAbsolutePath() + " | " + sysCmd);
		// this.setSystemCommandForPrint(null); // Reset after having consumed sys cmd. 

		ProcessBuilder pb = new ProcessBuilder().command(cmd);
		pb.redirectErrorStream(true);
		Process p= pb.start();
		
		BufferedReader reader= new BufferedReader(new InputStreamReader(p.getInputStream()));
    
		List<String> outRecords= new ArrayList<String>();
		String line = "";
		while ((line = reader.readLine())!= null) {
			outRecords.add(line);
		}
		reader.close();

		try {
			p.waitFor();
		} catch (InterruptedException e) {
			e.printStackTrace();
		}
        tmp.delete();
		return outRecords;
	}

	public void setSystemCommandForPrint(String systemCommandForPrint){
		this.systemCommandForPrint= systemCommandForPrint; 
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

	protected boolean getPrintNormalizedVcf(){
		return this.printNormalizedVcf;
	}
	
	protected void setPrintNormalizedVcf(boolean printNormalizedVcf){
		this.printNormalizedVcf= printNormalizedVcf;
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
	
	/**Returns a list of boolean indicating whether the reads in samReader pass the 
	 * sam and awk filters.
	 * @throws IOException 
	 * */
	public List<Boolean> filterReads(SamReader samReader) throws IOException{

		Iterator<SAMRecord> filterSam= samReader.query(this.getGc().getChrom(), this.getGc().getFrom(), this.getGc().getTo(), false);
		
		AggregateFilter aggregateFilter= new AggregateFilter(this.getSamRecordFilter());

		// This array will contain true/false to indicate whether a record passes the 
		// sam filters AND the awk filter (if given).
		// boolean[] results= new boolean[(int) this.nRecsInWindow];
		List<Boolean> results= new ArrayList<Boolean>();
		
		StringBuilder awkDataInput= new StringBuilder();
		while(filterSam.hasNext()){ 
			// Record whether a read passes the sam filters. If necessary, we also 
			// store the raw reads for awk.
			SAMRecord rec= filterSam.next();
			boolean passed;
			if(!rec.getReadUnmappedFlag() && 
			        !aggregateFilter.filterOut(rec) &&
			        rec.getAlignmentEnd() >= rec.getAlignmentStart()){
				passed= true;
			} else {
				passed= false;
			}
			if(passed){
				// grep
				String raw= rec.getSAMString().replaceAll("\n$", "");
				boolean showIt= true;
				if(this.showRegex != null && ! this.showRegex.equals(".*")){
					showIt= Pattern.compile(this.showRegex).matcher(raw).find();
				}
				boolean hideIt= false;
				if(!this.hideRegex.isEmpty()){
					hideIt= Pattern.compile(this.hideRegex).matcher(raw).find();	
				}
				if(!showIt || hideIt){
					passed= false;
				}
			}
			results.add(passed);
			if(passed && this.getAwk() != null && ! this.getAwk().isEmpty()){
				// We pass to awk only records that have been kept so far.
				awkDataInput.append(rec.getSAMString());	
			}
		}
		
		// Apply the awk filter, if given
		if(this.getAwk() != null && ! this.getAwk().isEmpty()){
			String[] rawLines= awkDataInput.toString().replaceAll("\n$", "").split("\n");
			boolean[] awkResults= Utils.passAwkFilter(rawLines, this.getAwk());
			// Compare the results array with awk filtered. Flip as appropriate the results array
			int awkIdx= 0;
			int i= 0;
			for(boolean isPassed : results){
				if(isPassed){
					if( ! awkResults[awkIdx]){
						results.set(i, false);
					}
					awkIdx++;
				}
				i++;
			}
		}
		return results;
	}

	protected void setColorForRegex(List<Argument> xcolorForRegex) {
		
	}

	public GenotypeMatrix getGenotypeMatrix() {
		return genotypeMatrix;
	}

	/** Iterate through the features in this track and set background colour.
	 * colorForRegex: Key= Regex to capture features; Value= Colour to use for the captures features.
	 * @throws InvalidColourException 
	 * */
	protected void changeFeatureColor(List<Argument> list) throws InvalidColourException {
		
	}

	protected void setLastModified() throws IOException {
		UrlValidator urlValidator = new UrlValidator();
		if(urlValidator.isValid(this.getFilename())){
			URL url = new URL(this.getFilename());
			HttpURLConnection httpCon = (HttpURLConnection) url.openConnection();
		    this.lastModified= httpCon.getLastModified();
		} else {
			this.lastModified= new File(this.getFilename()).lastModified();
		}
	}
	protected long getLastModified(){
		return this.lastModified;
	}

	public boolean getReadsAsPairs(){
		return this.readsAsPairs;
	}

	public void setReadsAsPairs(boolean readsAsPairs) throws InvalidGenomicCoordsException, IOException {
		this.readsAsPairs= readsAsPairs;
	}
	
}

