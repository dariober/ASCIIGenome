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

import coloring.Config;
import coloring.ConfigKey;
import coloring.Xterm256;
import exceptions.InvalidColourException;
import exceptions.InvalidCommandLineException;
import exceptions.InvalidGenomicCoordsException;
import exceptions.InvalidRecordException;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.filter.AggregateFilter;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
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
	private boolean hideTrack= false; 
	private boolean hideTitle= false;
	private TrackFormat trackFormat;
	 
	private int printRawLineCount= -1; // Number of lines to print. Same as `head -n 10`
	private GenotypeMatrix genotypeMatrix= new GenotypeMatrix();
	/** A file to export track data
	 * */
	private String exportFile= null;
	private String systemCommandForPrint;
	private boolean printNormalizedVcf= false;
	private long lastModified;
	
	private FeatureFilter featureFilter= new FeatureFilter(); 
	
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


	void setSamRecordFilter(List<SamRecordFilter> samRecordFilter) throws MalformedURLException, ClassNotFoundException, IOException, InvalidGenomicCoordsException, InvalidRecordException, SQLException {
		this.getFeatureFilter().setSamRecordFilter(samRecordFilter);
		this.update();
	}
	
	public boolean isBisulf() { return this.bisulf; }
	public void setBisulf(boolean bisulf) { this.bisulf= bisulf; }

	public void setAwk(String awk) throws ClassNotFoundException, IOException, InvalidGenomicCoordsException, InvalidRecordException, SQLException {
		this.getFeatureFilter().setAwk(awk);
		this.update();
	};
	
	public String getAwk(){
		return this.getFeatureFilter().getAwk();
	};
	
	protected FeatureFilter getFeatureFilter(){
		return this.featureFilter;
	}
	
	/** Setter for both showRegex and hideRegex, if only one is set, use
	 * Tracks.SHOW_REGEX or Tracks.HIDE_REGEX for the other. This is to prevent
	 * calling update() twice when only one is needed.*/
	public void setShowHideRegex(String showRegex, String hideRegex) throws InvalidGenomicCoordsException, IOException, ClassNotFoundException, InvalidRecordException, SQLException{
		this.getFeatureFilter().setShowHideRegex(showRegex, hideRegex);
		this.update();
	}

	public String getHideRegex() { 
		return this.getFeatureFilter().getHideRegex();
	}
	public String getShowRegex() {
		return this.getFeatureFilter().getShowRegex();
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
		return this.getFeatureFilter().get_f_flag();
	}

	/** This int is just a setting but is NOT translated to a filter! */
	protected void set_f_flag(int f_flag) {
		this.getFeatureFilter().set_f_flag(f_flag);
	}

	/** You should use a converter to get int from list of filters. 
	 * */
	/** This int is just a setting but is NOT translated to a filter! */
	protected int get_F_flag() {
		return this.getFeatureFilter().get_F_flag();
	}

	/** This int is just a setting but is NOT translated to a filter! */
	protected void set_F_flag(int F_flag) {
		this.getFeatureFilter().set_F_flag(F_flag);
	}

	/** This int is just a setting but is NOT translated to a filter! */
	public int getMapq() {
		return this.getFeatureFilter().getMapq();
	}

	/** This int is just a setting but is NOT translated to a filter! */
	protected void setMapq(int mapq) {
		this.getFeatureFilter().setMapq(mapq);
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
	
	public List<Boolean> filterReads(SamReader samReader, String chrom, int from, int to) throws IOException{

		Iterator<SAMRecord> filterSam= samReader.query(chrom, from, to, false);
		
		AggregateFilter aggregateFilter= new AggregateFilter(this.getFeatureFilter().getSamRecordFilter());

		// This array will contain true/false to indicate whether a record passes the 
		// sam filters AND the awk filter (if given).
		// boolean[] results= new boolean[(int) this.nRecsInWindow];
		List<Boolean> results= new ArrayList<Boolean>();

		byte[] faSeq= null;
		if(this.getFeatureFilter().getVariantChrom().equals(this.getGc().getChrom()) &&
		   this.getFeatureFilter().getVariantFrom() <= this.getGc().getTo() &&
		   this.getGc().getFrom() <= this.getFeatureFilter().getVariantTo()){
			// Get sequence to test whether the read is variant in the requested interval.
			IndexedFastaSequenceFile faSeqFile = new IndexedFastaSequenceFile(new File(this.getGc().getFastaFile()));
			faSeq= faSeqFile.getSubsequenceAt(
					this.getFeatureFilter().getVariantChrom(), 
					this.getFeatureFilter().getVariantFrom(), 
					this.getFeatureFilter().getVariantTo()).getBases();
			faSeqFile.close();
		}
		List<String> awkDataInput= new ArrayList<String>();
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
			
			// Filter for variant
			if(faSeq != null){
				int varFrom= this.getFeatureFilter().getVariantFrom();
				int varTo= this.getFeatureFilter().getVariantTo();
				// If above we retrieved the sequence it means some positions need to be scanned for variant reads
				passed= false; // Change to true as soon as we find a SNP in this read in the given interval
				int readPos= 0;
				int refPos= rec.getAlignmentStart();
				for(CigarElement cigar : rec.getCigar().getCigarElements()){
					if(cigar.getOperator().equals(CigarOperator.SOFT_CLIP)){
						readPos += cigar.getLength();
					}
					else if(cigar.getOperator().equals(CigarOperator.MATCH_OR_MISMATCH)){
						for(int i= 0; i < cigar.getLength(); i++){
							if(refPos >= varFrom && refPos <= varTo && rec.getReadLength() > 0){
								byte readBase= rec.getReadBases()[readPos];
								byte refBase= faSeq[refPos-varFrom];
								if(readBase != refBase){
									passed= true;
									break;
								}
							}
							readPos++;
							refPos++;
						}
					}
					else if(cigar.getOperator().equals(CigarOperator.DELETION)){ // Consumes ref base, not read base
						// REF  ACTGTTTTACTG
						// READ   TG----AC
						//          ^^^^
						for(int i= 0; i < cigar.getLength(); i++){
							if(refPos >= varFrom && refPos <= varTo){
								passed= true;
								break;
							}
							refPos++;							
						}
					}
					else if(cigar.getOperator().equals(CigarOperator.INSERTION)){ // Consumes read, not ref 
						//  REF ACTG----ACTG
						// READ   TGttttAC
						//         ^    
						for(int i= 0; i < cigar.getLength(); i++){
							if(refPos >= varFrom && refPos <= varTo){
								passed= true;
								break;
							}
							readPos++;							
						}
					} 
					else if(cigar.getOperator().equals(CigarOperator.HARD_CLIP)){ 
						//
					} 
					else if(cigar.getOperator().equals(CigarOperator.SKIPPED_REGION)){ // Same deletion but it's not a mismatch 
						refPos += cigar.getLength();
					} 
					else if(cigar.getOperator().equals(CigarOperator.PADDING)){ 
						// Not sure what to do with this...
					} 
					if(passed){
						break;
					}
				}
			}
			
			String raw= null;

			if(passed && (! this.getFeatureFilter().getShowRegex().equals(FeatureFilter.SHOW_REGEX) || 
					      ! this.getFeatureFilter().getHideRegex().equals(FeatureFilter.HIDE_REGEX))){
				// grep
				raw= rec.getSAMString().trim();
				boolean showIt= true;
				if(! this.getFeatureFilter().getShowRegex().equals(FeatureFilter.SHOW_REGEX)){
					showIt= Pattern.compile(this.getFeatureFilter().getShowRegex()).matcher(raw).find();
				}
				boolean hideIt= false;
				if(! this.getFeatureFilter().getHideRegex().equals(FeatureFilter.HIDE_REGEX)){
					hideIt= Pattern.compile(this.getFeatureFilter().getHideRegex()).matcher(raw).find();	
				}
				if(!showIt || hideIt){
					passed= false;
				}
			}
			results.add(passed);
			if(passed && this.getAwk() != null && ! this.getAwk().isEmpty()){
				// We pass to awk only records that have been kept so far.
				if(raw == null){
					raw= rec.getSAMString().trim();
				}
				awkDataInput.add(raw);
			}
		}

		// Apply the awk filter, if given
		if(this.getAwk() != null && ! this.getAwk().isEmpty()){
			String[] rawLines= new String[awkDataInput.size()];
			rawLines= awkDataInput.toArray(rawLines);
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
	
	/**Returns a list of boolean indicating whether the reads in samReader pass the 
	 * sam and awk filters.
	 * @throws IOException 
	 * */
//	public List<Boolean> filterReads(SamReader samReader) throws IOException{
//		return this.filterReads(samReader, this.getGc().getChrom(), this.getGc().getFrom(), this.getGc().getTo());
//	}

//	private SAMRecord getNextSamRecordOnChrom(String chrom, int from) throws IOException, InvalidGenomicCoordsException{
//		
//		int qend= (from - 1) < 0 ? 0 : (from - 1);
//
//		SamReader samReader= Utils.getSamReader(this.getWorkFilename());
//		SAMRecordIterator iter = samReader.queryAlignmentStart(chrom, qend);
//		while(true){
//			SAMRecord line= iter.next();
//			if(line == null){
//				return null;
//			}
//			// Create a reader at the position of this read and test whether any read here passes
//			// filters.
//			SamReader cur= Utils.getSamReader(this.getWorkFilename());
//			List<Boolean> passed = this.filterReads(cur, line.getReferenceName(), line.getAlignmentStart(), line.getAlignmentStart());
//			
//		}
//		
//		// TabixBigBedIterator iter= this.getReader().query(chrom, qend, Integer.MAX_VALUE);
//		while(true){
//			SAMRecord line= iter.next();
//			if(line == null){
//				return null;
//			}
//			IntervalFeature x= new IntervalFeature(line, this.getTrackFormat(), this.getVCFCodec());
//			if(x.getFrom() > from && this.featureIsVisible(x.getRaw())){
//				return x;
//			}
//		}
//	}
	
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
	
	/**Set filter to extract reads containing variant at the given interval.
	 * from, to: 1-based coordinates (first base of chr1 is `chr1:1-1`).
	 * @throws SQLException 
	 * @throws InvalidRecordException 
	 * @throws InvalidGenomicCoordsException 
	 * @throws IOException 
	 * @throws ClassNotFoundException 
	 * @throws MalformedURLException 
	 * */
	public void setVariantReadInInterval(String chrom, int from, int to) throws MalformedURLException, ClassNotFoundException, IOException, InvalidGenomicCoordsException, InvalidRecordException, SQLException{
		if(from > to){
			System.err.println("Invalid coordinates for filter from > to: " + from + ", " + to);
			throw new InvalidGenomicCoordsException();
		}
		this.getFeatureFilter().setVariantReadInInterval(chrom, from, to);
		this.update();
	}
	
}

