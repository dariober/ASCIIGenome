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
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.commons.lang3.StringUtils;
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
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.filter.AggregateFilter;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import samTextViewer.GenomicCoords;
import samTextViewer.Main;
import samTextViewer.Utils;

public abstract class Track {

	public static String awkFunc= "";
	static {
		  try {
			  InputStream in= Main.class.getResourceAsStream("/functions.awk");
			  BufferedReader br= new BufferedReader(new InputStreamReader(in));
			  String line = null;
			  StringBuilder x= new StringBuilder();
			  while((line = br.readLine()) != null) {
			      if(line.trim().startsWith("#")){
			    	  continue;
			      }
				  x.append(line + '\n');
			  }
			  awkFunc= x.toString();
			  br.close();
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
	protected List<Float> screenScores= new ArrayList<Float>();
	private GenomicCoords gc;
	private boolean noFormat= false; 
	private float yLimitMin= Float.NaN; // Same as R ylim()
	private float yLimitMax= Float.NaN;
	/** Max size of genomic region before the track shuts down to prevent excessive slow down */
	protected final int MAX_REGION_SIZE= 1000001;   
	
	protected String titleColour= null;
	protected boolean bisulf= false;

	private String gtfAttributeForName= null;
	/** Should features on with same coords be squashed into a single one? */
	private PrintRawLine printMode= PrintRawLine.OFF;
	/** Number of decimal places to show when printing raw lines */
	private int printNumDecimals= 3;
	private boolean explainSamFlag= false;
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
	private String printFormattedVep= null;
	private boolean printNormalizedVcf= false;
	private long lastModified;
	
	private FeatureFilter featureFilter= new FeatureFilter(); 
	private VCFHeader vcfHeader;
	private Pattern highlightPattern;
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
		return Utils.reformatFileName(trackTag, false);
	}
	
	public void setTrackTag(String trackTag) { 
		this.trackTag = trackTag; 
	}
	
	protected List<Float> getScreenScores() {
		return screenScores;
	}
	protected void setScreenScores(List<Float> screenScores) {
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

	public Float getYLimitMin() { 
		return yLimitMin;
	}
	
	public void setYLimitMin(float ymin) { 
		this.yLimitMin = ymin;
	}

	public Float getYLimitMax() { 
		return yLimitMax; 
	}
	public void setYLimitMax(float ymax) { 
		this.yLimitMax = ymax; 
	}


	void setSamRecordFilter(List<SamRecordFilter> samRecordFilter) throws MalformedURLException, ClassNotFoundException, IOException, InvalidGenomicCoordsException, InvalidRecordException, SQLException {
		this.getFeatureFilter().setSamRecordFilter(samRecordFilter);
		this.update();
	}
	
	public boolean isBisulf() { return this.bisulf; }
	public void setBisulf(boolean bisulf) { this.bisulf= bisulf; }

	public void setAwk(String awk) throws ClassNotFoundException, IOException, InvalidGenomicCoordsException, InvalidRecordException, SQLException {
		
		if( ! awk.trim().isEmpty()){
			List<String> arglst= Utils.tokenize(awk, " ");
			
			// Do we need to set tab as field sep?
			if(arglst.size() == 1 || ! arglst.contains("-F")){ // It would be more stringent to check for the script.
				awk= "-F '\\t' " + awk; 
			}
		}
		this.getFeatureFilter().setAwk(awk);
		this.update();
	}
	
	public String getAwk(){
		return this.getFeatureFilter().getAwk();
	}
	
	protected FeatureFilter getFeatureFilter(){
		return this.featureFilter;
	}
	
	/** Setter for both showRegex and hideRegex, if only one is set, use
	 * Tracks.SHOW_REGEX or Tracks.HIDE_REGEX for the other. This is to prevent
	 * calling update() twice when only one is needed.*/
	public void setShowHideRegex(Pattern showRegex, Pattern hideRegex) throws InvalidGenomicCoordsException, IOException, ClassNotFoundException, InvalidRecordException, SQLException{
		this.getFeatureFilter().setShowHideRegex(showRegex, hideRegex);
		this.update();
	}

	public Pattern getHideRegex() { 
		return this.getFeatureFilter().getHideRegex();
	}
	public Pattern getShowRegex() {
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
		for(int i= 0; i < rawList.size(); i++){
			rawList.set(i, rawList.get(i).trim());
		}
		
		if(this.getExportFile() != null && ! this.getExportFile().isEmpty()){
			// If an output file has been set, send output there and return. 
			BufferedWriter wr= null;
			try{
				wr = new BufferedWriter(new FileWriter(this.getExportFile(), true));
				for(String line : rawList){
					wr.write(line.replace("\n", "").replace("\r", "") + "\n");
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
			
			if(this.getTrackFormat().equals(TrackFormat.BAM) && this.explainSamFlag){
				line= this.explainSamFlag(line);
			}
			
			if(this.getTrackFormat().equals(TrackFormat.VCF) && this.getPrintFormattedVep() != null){
				String vepTag= this.getPrintFormattedVep().isEmpty() ? "CSQ" : this.getPrintFormattedVep(); 
				line= this.printFormattedVep(line, vepTag);
			}
			
			if(this.getPrintMode().equals(PrintRawLine.CLIP) && (this.getSystemCommandForPrint() == null || this.getSystemCommandForPrint().isEmpty())){
				
				line= this.shortenPrintLine(line);
				
			}
			line= Utils.roundNumbers(line, this.getPrintNumDecimals(), this.getTrackFormat());
			
			if( this.getHighlightPattern() != null && ! this.getHighlightPattern().pattern().isEmpty()){
				line= this.highlightPattern(line, this.getHighlightPattern());
				if(this.getTrackFormat().equals(TrackFormat.VCF)){
					line= this.highlightVcfFormat(line, this.getHighlightPattern());
				}
			} 
			
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
		
		// Trim to fit size ignoring formatting
		for(String x : tabList){
			StringBuilder line= new StringBuilder();
			boolean isAnsiEscape= false;
			int nAnsiEscape= 0;
			int nchars= 0;
			for(int i = 0; i < x.length(); i++){
			    char c= x.charAt(i);
			    if(c == '\033'){
			    	isAnsiEscape= true;
			    	nAnsiEscape++;
			    }
			    if( ! isAnsiEscape ){
			    	nchars++;
			    }
			    if(isAnsiEscape && c == 'm'){
			    	isAnsiEscape= false;
			    }
			    line.append(c);
			    if(nchars >= windowSize){
			    	break;
			    }
			}
			if(nAnsiEscape % 2 != 0){
				line.append("\033[27m"); // If the closing escape have been trimmed off, put it back
			}
			sb.append(line.toString() + "\n");
		}
		
		if(this.isNoFormat()){
			return sb.toString();
		} else {
			String formatted=  "\033[38;5;" + Config.get256Color(ConfigKey.foreground) + 
			";48;5;" + Config.get256Color(ConfigKey.background) + "m" + sb.toString();
			
			return formatted; 
		}
	}

	/** Parse line assuming it is a SAM record and add explanation to sam flag;
	 * */
	private String explainSamFlag(String line) {
		
		String[] xline= line.split("\t");
		
		int flag= Integer.valueOf(xline[1]);
		StringBuilder ex= new StringBuilder();
		ex.append(flag);
		ex.append("|");
		if((flag & 16) == 16) {ex.append("-|");} else {ex.append("+|");}
		if((flag & 1) == 1) ex.append("pair|");
		if((flag & 2) == 2) ex.append("prop-p|");
		if((flag & 4) == 4) ex.append("unmapped|");
		if((flag & 8) == 8) ex.append("mate-unmap|");
		if((flag & 32) == 32) ex.append("mate-rev|");
		if((flag & 64) == 64) ex.append("1st|");
		if((flag & 128) == 128) ex.append("2nd|");
		if((flag & 256) == 256) ex.append("not-primary|");
		if((flag & 512) == 512) ex.append("QC-fail|");
		if((flag & 1024) == 1024) ex.append("dupl|");
		if((flag & 2048) == 2048) ex.append("suppl|");
		
		xline[1]= ex.toString().replaceAll("\\|$", "");
		
		return Joiner.on("\t").join(xline);
	}

	/** Print VCF record with VEP annotation in a readable format;
	 * */
	private String printFormattedVep(String line, String vepTag) {
		
		List<String> vepArgs = this.parseVepTag(vepTag);
		String infoTag= vepArgs.get(0);
		List<String> requiredHeader= new ArrayList<String>();
		for(int i= 1; i < vepArgs.size(); i++){
			requiredHeader.add(vepArgs.get(i).toLowerCase());
		}
		
		VCFInfoHeaderLine vcfHdr = this.getVcfHeader().getInfoHeaderLine(infoTag);
		if(vcfHdr == null){
			return line;
		}
		String description= vcfHdr.getDescription();
		Map<Integer, String> hdrIndexToName= new HashMap<Integer, String>();
		List<String> header= Splitter.on("|").splitToList(description.replaceAll(".*Format: ", ""));
		for(int i= 0; i < header.size(); i++){
			hdrIndexToName.put(i, header.get(i).toLowerCase());
		}
		int pad= 0;
		for(String x : hdrIndexToName.values()){
			if(requiredHeader.size() == 0 || requiredHeader.contains(x)){
				if(x.length() > pad){
					pad= x.length(); 
				}
			}
		}
				
		VCFCodec vcfCodec= new VCFCodec();
		vcfCodec.setVCFHeader(this.getVcfHeader(), Utils.getVCFHeaderVersion(this.getVcfHeader()));
		VariantContext variantContext= vcfCodec.decode(line);
		
		List<String> annotation= variantContext.getAttributeAsStringList(infoTag, "");
		
		String plugable= "...";
		if( ! requiredHeader.contains("null")){ 
			List<String> fmtCsq= new ArrayList<String>();
			for(String feature : annotation){
				List<String> lst= Splitter.on("|").splitToList(feature);
				List<String> fmtFeature= new ArrayList<String>();
				for(int i= 0; i < lst.size(); i++){
					if(lst.get(i).isEmpty()){
						continue;
					}
					if(requiredHeader.size() > 0 && ! requiredHeader.contains(hdrIndexToName.get(i))){
						continue;
					}
					String hdr= header.get(i);
					hdr= StringUtils.rightPad(hdr, pad, " ");
					fmtFeature.add("" + hdr + ": " + lst.get(i));
				}
				fmtCsq.add(Joiner.on("\n").join(fmtFeature));
			}
			plugable= "\n" + Joiner.on("\n----\n").join(fmtCsq) + "\n";
		}
		String[] xline= line.split("\t");
		xline[7]= xline[7].replace(Joiner.on(",").join(annotation), plugable); 
		
		return Joiner.on("\t").join(xline).trim();
	}
	
	/**Parse vepTag containing the INFO tag for the VEP annotation (typically CSQ)
	 * and an optional list of header names to keep for printing.*/
	private List<String> parseVepTag(String vepTag){
		List<String> vepArgs= Splitter.on(",").splitToList(vepTag);
		return vepArgs;
	}
	
	/**Abbreviate printable fields such as long sequences in VCF 
	 * alleles and read sequence and quality strings in BAM.
	 * */
	private String shortenPrintLine(String line){
		
		if(this.getTrackFormat().equals(TrackFormat.BAM)){
			// Make SEQ and QUAL shorter
			String[] bamrec= line.split("\t");
			int max= 5;
			if(bamrec[9].length() > 20){
				bamrec[9]= bamrec[9].substring(0, max) + "[+" + (bamrec[9].length()-max)  + "]";
			}
			if(bamrec[10].length() > 20){
				bamrec[10]= bamrec[10].substring(0, max)  + "[+" + (bamrec[10].length()-max)  + "]";
			}
			line= Joiner.on("\t").join(bamrec);
		}
		if(this.getTrackFormat().equals(TrackFormat.VCF)){
			// Make REF and ALT shorter
			int max= 5;
			String[] vcfrec= line.split("\t");
			if(vcfrec[3].length() > 20 && vcfrec[3].toLowerCase().matches("[actgn]+")){
				vcfrec[3]= vcfrec[3].substring(0, max) + "[+" + (vcfrec[3].length()-max)  + "]";
			}
			if(vcfrec[4].length() > 20 && vcfrec[4].toLowerCase().matches("[actgn]+")){
				vcfrec[4]= vcfrec[4].substring(0, max)  + "[+" + (vcfrec[4].length()-max)  + "]";
			}
			line= Joiner.on("\t").join(vcfrec);
		}
		// Trim any field that is really long
		String[] rec= line.split("\t");
		for(int i= 0; i < rec.length; i++){
			if(rec[i].length() > 150 && ! rec[i].contains(" ")){
				rec[i]= rec[i].substring(0, 5) + "..." + rec[i].substring(rec[i].length()-5); 
			}
		}
		line= Joiner.on("\t").join(rec);
		return line;
	}
	
	/**Find the values in the VCF record whose TAG matches pattern and add ANSI highlighting. 
	 * */
	private String highlightVcfFormat(String vcf, Pattern p){
		
		List<String> xvcf = Splitter.on("\t").splitToList(vcf);
		if(xvcf.size() < 10){
			return vcf;
		}
		
		// Get the index of the value(s) to be formatted
		List<String> tags= Splitter.on(":").splitToList(Utils.stripAnsiCodes(xvcf.get(8))); 
		List<Integer> hlValueIdx= new ArrayList<Integer>();
		for(int i= 0; i < tags.size(); i++){
			// If format tags are GT:AF:AC and pattern is AF, you'll get hlValueIdx= [1] 
			if(p.matcher(tags.get(i)).find()){
				hlValueIdx.add(i);
			}
		}
		if(hlValueIdx.size() == 0){
			return vcf; // Nothing replace
		}
		
		List<String> o= xvcf.subList(0, 9); // Constant fields
		List<String> outvcf= new ArrayList<String>(o);
		// For each sample, format the values at index(es) stored in hlValueIdx 
		for(int i= 9; i < xvcf.size(); i++){
			String sample= xvcf.get(i); // String of values for this sample, e.g. 0|1:0.14:FR
			List<String> v= Splitter.on(":").splitToList(sample);
			List<String> values= new ArrayList<String>(v);
			for(int idx : hlValueIdx){
				if(idx < values.size()){
					values.set(idx, "\033[7m" + values.get(idx) + "\033[27m");
				}
			}
			String fmtSample= Joiner.on(":").join(values);
			outvcf.add(fmtSample);
		}
		return Joiner.on("\t").join(outvcf);
	}
	
	/**Put some ANSI highlighting to string x where there is a match to pattern
	 * */
	private String highlightPattern(String x, Pattern p){
		Matcher m= p.matcher(x);
		StringBuilder sb= new StringBuilder();
		int idx= 0;
		while(m.find()){
			String prefix= x.substring(idx, m.start());
			sb.append(prefix);
			sb.append("\033[7m");
			sb.append(m.group());
			sb.append("\033[27m");
			idx= m.end();
		}
		String prefix= x.substring(idx, x.length());
		sb.append(prefix);
		return sb.toString();
	}
	
	private Pattern getHighlightPattern(){
		return this.highlightPattern;
	}
	
	public void setHighlightPattern(Pattern x){
		this.highlightPattern= x;
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
		
		if(this.trackFormat.equals(TrackFormat.VCF)){
			List<String> vcf= Utils.vcfHeaderToStrings(this.getVcfHeader());
			vcf.addAll(recordsAsStrings);
			recordsAsStrings= vcf;
		}
		File tmp= Utils.createTempFile(".asciigenome.", ".print.tmp");
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
			if(this.trackFormat.equals(TrackFormat.VCF) && line.startsWith("#")){
				continue;
			}
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

	/**Return a string to be plugged into the title line showing which filters are active on the track*/
	protected abstract String getTitleForActiveFilters();
	
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
	
	public List<Boolean> filterReads(SamReader samReader, String chrom, int from, int to) throws IOException {

		Iterator<SAMRecord> filterSam= samReader.query(chrom, from, to, false);
		
		AggregateFilter aggregateFilter= new AggregateFilter(this.getFeatureFilter().getSamRecordFilter());

		// This array will contain true/false to indicate whether a record passes the 
		// sam filters AND the awk filter (if given).
		// boolean[] results= new boolean[(int) this.nRecsInWindow];
		List<Boolean> results= new ArrayList<Boolean>();

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

			// Filter for variant reads: Do it only if there is an intersection between variant interval and current genomic window
			if(this.getFeatureFilter().getVariantChrom().equals(Filter.DEFAULT_VARIANT_CHROM.getValue())){
				// Variant read filter is not set. Nothing to do.
			}
			else if(passed){ 
				// Memo: Test filter(s) only if a read is passed==true. 
				passed= this.isSNVRead(rec, this.getFeatureFilter().isVariantOnly());
			}
			
			String raw= null;

			if(passed && (! this.getFeatureFilter().getShowRegex().equals(Filter.DEFAULT_SHOW_REGEX.getValue()) || 
					      ! this.getFeatureFilter().getHideRegex().equals(Filter.DEFAULT_HIDE_REGEX.getValue()))){
				// grep
				raw= rec.getSAMString().trim();
				boolean showIt= true;
				if(! this.getFeatureFilter().getShowRegex().equals(Filter.DEFAULT_SHOW_REGEX.getValue())){
					showIt= this.getFeatureFilter().getShowRegex().matcher(raw).find();
				}
				boolean hideIt= false;
				if(! this.getFeatureFilter().getHideRegex().equals(Filter.DEFAULT_HIDE_REGEX.getValue())){
					hideIt= this.getFeatureFilter().getHideRegex().matcher(raw).find();	
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
		} // End loop through reads

		// Apply the awk filter, if given
		if(this.getAwk() != null && ! this.getAwk().equals(Filter.DEFAULT_AWK.getValue())){
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

	/**Return true if samrecord contains a mismatch or insertion/deletion in the target region.
	 * */
	private boolean isSNVRead(SAMRecord rec, boolean variantOnly) {
		boolean passed= false;
		
		int varFrom= this.getFeatureFilter().getVariantFrom();
		int varTo= this.getFeatureFilter().getVariantTo();

		if(this.getFeatureFilter().getVariantChrom().equals(rec.getReferenceName()) &&
		        varFrom <= rec.getAlignmentEnd() &&
		        rec.getAlignmentStart() <= varTo){
			// Variant read filter is set and this read overlaps it. 
			if( ! variantOnly){
				return true; // No need to check whether read is variant.
			}
			int readPos= 0;
			int refPos= rec.getAlignmentStart();
			for(CigarElement cigar : rec.getCigar().getCigarElements()){
				if(cigar.getOperator().equals(CigarOperator.SOFT_CLIP)){
					readPos += cigar.getLength();
				}
				else if(cigar.getOperator().equals(CigarOperator.M) ||
						cigar.getOperator().equals(CigarOperator.EQ) || 
						cigar.getOperator().equals(CigarOperator.X)){
					for(int i= 0; i < cigar.getLength(); i++){
						if(refPos >= varFrom && refPos <= varTo && rec.getReadLength() > 0){
							byte readBase= rec.getReadBases()[readPos];
							byte refBase= this.getFeatureFilter().getFaSeq()[refPos-varFrom];
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
		else {
			// Variant read filter is set and this read does not overlap it.
			passed= false;
		}
		return passed;
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
	
	/**Set filter to extract reads containing variant at the given interval.
	 * from, to: 1-based coordinates (first base of chr1 is `chr1:1-1`).
	 * @throws SQLException 
	 * @throws InvalidRecordException 
	 * @throws InvalidGenomicCoordsException 
	 * @throws IOException 
	 * @throws ClassNotFoundException 
	 * @throws MalformedURLException 
	 * */
	public void setVariantReadInInterval(String chrom, int from, int to, boolean variantOnly) throws MalformedURLException, ClassNotFoundException, IOException, InvalidGenomicCoordsException, InvalidRecordException, SQLException{
		this.getFeatureFilter().setVariantOnly(variantOnly);
		if(chrom.equals(Filter.DEFAULT_VARIANT_CHROM.getValue())){
			this.getFeatureFilter().setVariantReadInInterval(chrom, from, to, null);
			this.update();
			return;
		}
		if(from > to){
			System.err.println("Invalid coordinates for filter from > to: " + from + ", " + to);
			throw new InvalidGenomicCoordsException();
		}
		IndexedFastaSequenceFile faSeqFile = new IndexedFastaSequenceFile(new File(this.getGc().getFastaFile()));
		byte[] faSeq= faSeqFile.getSubsequenceAt(chrom, from, to).getBases();
		faSeqFile.close();
		this.getFeatureFilter().setVariantReadInInterval(chrom, from, to, faSeq);
		this.update();
	}
	
	/**
	 * @throws InvalidGenomicCoordsException Interpret the variantInterval string, like 123+-5, to return the coordinates 
	 * of the interval to filter for variants. Return an array of length 2 with
	 * [start, end], 1-based.
	 * @throws  
	 * */
//	private int[] variantIntervalToCoords(String variantInterval) throws InvalidGenomicCoordsException {
//
//		int from= -1;
//		int to= -1;
//		if(variantInterval != null){
//			
//			// Find whether the region is a single pos +/- a value or 
//			variantInterval= variantInterval.trim().replaceAll(" ", "").replaceAll(",", "").replaceAll("/", "");
//			try{
//				if(variantInterval.contains("+-") || variantInterval.contains("-+")){
//					String[] fromTo= variantInterval.replaceAll("\\+", "").split("-");
//					int offset= this.posFromGenomicOrScreen(fromTo[1]);
//					from= this.posFromGenomicOrScreen(fromTo[0]) - offset;
//					to= this.posFromGenomicOrScreen(fromTo[0]) + offset;
//				} 
//				else if(variantInterval.contains("-")){
//					String[] fromTo= variantInterval.split("-");
//					from= this.posFromGenomicOrScreen(fromTo[0]) - this.posFromGenomicOrScreen(fromTo[1]);
//					to= this.posFromGenomicOrScreen(fromTo[0]);
//				} 
//				else if(variantInterval.contains("+")){
//					String[] fromTo= variantInterval.split("\\+");
//					from= this.posFromGenomicOrScreen(fromTo[0]);
//					to= this.posFromGenomicOrScreen(fromTo[0]) + this.posFromGenomicOrScreen(fromTo[1]);
//				}
//				else if(variantInterval.contains(":")){
//					String[] fromTo= variantInterval.split(":");
//					from= this.posFromGenomicOrScreen(fromTo[0]);
//					to= this.posFromGenomicOrScreen(fromTo[1]);
//				}
//				else {
//					from= this.posFromGenomicOrScreen(variantInterval);
//					to= this.posFromGenomicOrScreen(variantInterval);
//				}
//			} catch(NumberFormatException e) {
//				try {
//					throw new NumberFormatException(Utils.padEndMultiLine("Cannot parse region into integers: " + variantInterval, Utils.getTerminalWidth()));
//				} catch (IOException e1) {
//					e1.printStackTrace();
//				}
//			}
//			if(from < 1) from= 1;
//			if(to < 1) to= 1;
//			if(from > to) from= to;
//		}
//		return new int[]{from, to};
//	}
	
	/**Argument position is parsable to either either an integer or a 
	 * decimal number. If decimal, return the genomic position corresponding to
	 * the screen percentile. E.g. 0.5 returns the genomic position at 50% of the screen.  
	 * @throws InvalidCommandLineException 
	 * @throws InvalidGenomicCoordsException 
	 * @throws IOException 
	 * */
//	private int posFromGenomicOrScreen(String position) throws InvalidGenomicCoordsException {
//		int reg;
//		try{
//			reg= Integer.parseInt(position);
//			return reg;
//		} catch(NumberFormatException e){
//			//
//		}
//		Double pct= Double.valueOf(position);
//		double offset= this.getGc().getGenomicWindowSize() * pct;
//		if(offset != 0){
//			offset--;
//		}
//		reg= this.getGc().getFrom() + (int)Math.round(offset);
//		return reg;
//	}

	
	protected VCFHeader getVcfHeader() {
		return vcfHeader;
	}

	protected void setVcfHeader(VCFHeader vcfHeader) {
		this.vcfHeader = vcfHeader;
	}

	public int getPrintNumDecimals() {
		return printNumDecimals;
	}

	public void setPrintNumDecimals(int printNumDecimals) {
		this.printNumDecimals = printNumDecimals;
	}
	
	public abstract void reload() throws InvalidGenomicCoordsException, IOException, ClassNotFoundException, InvalidRecordException, SQLException;

	/**Close readers associated to this track. 
	 * */
	public abstract void close();
	
	/**
	 * @param explainSamFlag the explainSamFlag to set
	 */
	protected void setExplainSamFlag(boolean explainSamFlag) {
		this.explainSamFlag = explainSamFlag;
	}

	public String getPrintFormattedVep() {
		return printFormattedVep;
	}

	public void setPrintFormattedVep(String printFormattedVep) {
		this.printFormattedVep = printFormattedVep;
	}
	
}

