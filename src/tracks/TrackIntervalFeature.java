package tracks;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Pattern;

import org.apache.commons.lang3.StringUtils;
import org.broad.igv.bbfile.BBFileReader;

import com.google.common.collect.Lists;

import exceptions.InvalidCommandLineException;
import exceptions.InvalidGenomicCoordsException;
import exceptions.InvalidRecordException;
import htsjdk.tribble.readers.TabixReader;
import samTextViewer.GenomicCoords;
import samTextViewer.Utils;
import sortBgzipIndex.MakeTabixIndex;
import ucsc.UcscGenePred;

public class TrackIntervalFeature extends Track {
 
	protected List<IntervalFeature> intervalFeatureList= new ArrayList<IntervalFeature>();  
	/**For GTF/GFF data: Use this attribute to get the feature names 
	 * */
	final static private String HIDE_REGEX= "^$";  private String hideRegex= HIDE_REGEX;
	final static private String SHOW_REGEX= ".*"; 	private String showRegex= SHOW_REGEX;
	protected TabixReader tabixReader; // Leave *protected* for TrackBookmark to work
	private BBFileReader bigBedReader;
	private int printRawLineCount= -1; // Number of lines to print. Same as `head -n 10`
	
	/* C o n s t r u c t o r */

	public TrackIntervalFeature(final String filename, GenomicCoords gc) throws IOException, InvalidGenomicCoordsException, ClassNotFoundException, InvalidRecordException, SQLException{
		
		this.setFilename(filename);

		if(Utils.isUcscGenePredSource(filename)){
			UcscGenePred ucsc = null;
			try {
				ucsc = new UcscGenePred(filename, -1);
				this.setWorkFilename(ucsc.getTabixFile());
				this.setTrackFormat(TrackFormat.GTF);
				//this.type= TrackFormat.GTF;
				this.tabixReader= new TabixReader(new File(this.getWorkFilename()).getAbsolutePath());

			} catch (InvalidCommandLineException e) {
				//
			}
		
		} else if(Utils.getFileTypeFromName(filename).equals(TrackFormat.BIGBED)){
			
			this.bigBedReader = new BBFileReader(filename);  // or url for remote access.
			if(!this.bigBedReader.getBBFileHeader().isBigBed()){
				throw new RuntimeException("File " + filename + " is not bigBed.");
			}
			
			this.setWorkFilename(filename);
			this.setTrackFormat(TrackFormat.BIGBED);
			//this.type= TrackFormat.BIGBED;
			
		} else if( ! Utils.hasTabixIndex(new File(filename).getAbsolutePath())){
			// Tabix index not found for this file. Sort and index input to tmp.

			String suffix= new File(filename).getName();
			if( ! suffix.endsWith(".gz")){
				suffix += ".gz";
			}
			String tmpWorkFile= File.createTempFile("asciigenome.", "." + suffix).getAbsolutePath();
			new File(tmpWorkFile).deleteOnExit();
			new File(new File(tmpWorkFile).getAbsolutePath() + ".tbi").deleteOnExit();
			this.setWorkFilename(tmpWorkFile);

			this.setTrackFormat(Utils.getFileTypeFromName(new File(filename).getName()));
			//this.type= Utils.getFileTypeFromName(new File(filename).getName());
			new MakeTabixIndex(filename, new File( this.getWorkFilename() ), Utils.trackFormatToTabixFormat(this.getTrackFormat()));	

			this.setWorkFilename(tmpWorkFile);
			this.tabixReader= new TabixReader(new File(this.getWorkFilename()).getAbsolutePath());
			
		} else { // This means the input is tabix indexed.
			this.setWorkFilename(filename);
			this.setTrackFormat(Utils.getFileTypeFromName(new File(filename).getName()));
			// this.type= Utils.getFileTypeFromName(new File(filename).getName());
			this.tabixReader= new TabixReader(new File(this.getWorkFilename()).getAbsolutePath());
		}
		this.setGc(gc);
	}
	

	protected TrackIntervalFeature(GenomicCoords gc){
		
	}
	
	/* M e t h o d s */
	@Override
	/** Collect features mapping to the current genomic coordinates and update the list of interval features
	 * for the track. Also update the mapping of features to the terminal screen.
	 * This method should be called only from within methods that change the features being displayed. 
	 * E.g. setGc(), which changes the coordinates, or setHideRegex() & setShowRegex() which change the visible 
	 * features. 
	 * update() should not change anything other than the list of features and the mapping. 
	 * */
	protected void update() throws IOException, InvalidGenomicCoordsException, ClassNotFoundException, InvalidRecordException, SQLException{
		this.intervalFeatureList = this.getFeaturesInInterval(
				this.getGc().getChrom(), this.getGc().getFrom(), this.getGc().getTo());
		
		int windowSize= this.getGc().getUserWindowSize();
		for(IntervalFeature ift : intervalFeatureList){
			ift.mapToScreen(this.getGc().getMapping(windowSize));
		}
	}
	
	public List<IntervalFeature> getFeaturesInInterval(String chrom, int from, int to) throws IOException, InvalidGenomicCoordsException{

		if(from < 1){
			System.err.println("from < 1: " + from + "; resetting to 1."); 
			from= 1;
		}
		
		if(from > to || to < 1){
			System.err.println("Invalid coordinates: from: " + from + "; to: " + to 
					+ "; Resetting to initial 1-" + Integer.MAX_VALUE);
			from= 1;
			to= Integer.MAX_VALUE;
			throw new InvalidGenomicCoordsException();
		}		
		List<IntervalFeature> xFeatures= new ArrayList<IntervalFeature>();

		TabixBigBedIterator qry= this.getReader().query(chrom, from-1, to);
		while(true){
			String q = qry.next();
			if(q == null){
				break;
			}
			IntervalFeature intervalFeature= new IntervalFeature(q, this.getTrackFormat());

			if(intervalFeature.getRaw().contains("\t__ignore_me__")){ // Hack to circumvent issue #38
				continue;
			}
			xFeatures.add(intervalFeature);
		} 
		
		// Remove hidden features
		List<IntervalFeature> xFeaturesFiltered= new ArrayList<IntervalFeature>();
		for(IntervalFeature x : xFeatures){
			if(this.featureIsVisible(x.getRaw())){
				xFeaturesFiltered.add(x);
			}
		}
		return xFeaturesFiltered;
	}


	/** Return true if string is visible, i.e. it
	 * passes the regex filters. Note that regex filters are applied to the raw string.
	 * */
	protected boolean featureIsVisible(String x){
		
		if(x.contains("__ignore_me__")){
			return false;
		}
		
		boolean showIt= Pattern.compile(this.showRegex).matcher(x).find();
		boolean hideIt= false;
		if(!this.hideRegex.isEmpty()){
			hideIt= Pattern.compile(this.hideRegex).matcher(x).find();	
		}
		if(showIt && !hideIt){
			return true;
		} else {
			return false;
		} 
	}

	/**Return the coordinates of the next feature so that the start coincide with the start of the feature and
	 * the end is the start + windowSize.  
	 * */
	public GenomicCoords coordsOfNextFeature(GenomicCoords currentGc, boolean getPrevious) throws InvalidGenomicCoordsException, IOException {
		
		IntervalFeature nextFeature;
		if(getPrevious){
			nextFeature= this.getPreviousFeature(currentGc.getChrom(), currentGc.getFrom());
		} else {
			nextFeature= this.getNextFeature(currentGc.getChrom(), currentGc.getTo());
		}
		if(nextFeature == null){
			return currentGc;
		}
		GenomicCoords nextGc= new GenomicCoords(
				Utils.coordinatesToString(nextFeature.getChrom(), 
						                  nextFeature.getFrom(), 
						                  nextFeature.getFrom() + currentGc.getGenomicWindowSize() -1),  
				currentGc.getSamSeqDict(),
				currentGc.getFastaFile());
		return nextGc;
	}

	private IntervalFeature getNextFeature(String startChrom, int from) throws IOException, InvalidGenomicCoordsException{
		
		IntervalFeature next = getNextFeatureOnChrom(startChrom, from);
		if(next != null){
			return next;
		} // There is no feature left on the starting chrom at the given position.

		// Search the remaining chroms in order starting from position 1:
		List<String> chroms= this.getChromListStartingAt(startChrom);
		if(chroms.contains(startChrom)){ 
			// The start chroms is searched last. The `if` statement controls for the case where
			// the startChrom is not present at all in this track.
			chroms.remove(startChrom); 	
			chroms.add(startChrom);			
		}
		for(String chrom : chroms){
			next = getNextFeatureOnChrom(chrom, 1);
			if(next != null){
				return next;
			}			
		}
		return null;
	}

	/** Get the next feature on chrom after "from" position or null if no 
	 * feature found. This function should be used only by getNextFeature() which is
	 * more general since it searches all the chromosomes in turn.    
	 * @throws IOException 
	 * @throws InvalidGenomicCoordsException */
	private IntervalFeature getNextFeatureOnChrom(String chrom, int from) throws IOException, InvalidGenomicCoordsException{
		
		TabixBigBedIterator iter= this.getReader().query(chrom, from-1, Integer.MAX_VALUE);
		while(true){
			String line= iter.next();
			if(line == null){
				return null;
			} 
			IntervalFeature x= new IntervalFeature(line, this.getTrackFormat());
			if(x.getFrom() > from && this.featureIsVisible(x.getRaw())){
				return x;
			}
		}
	}

	private IntervalFeature getPreviousFeature(String startChrom, int pos) throws IOException, InvalidGenomicCoordsException{
		
		IntervalFeature prev = getPreviousFeatureOnChrom(startChrom, pos);
		if(prev != null){
			return prev;
		} 
		// There is no feature left on the starting chrom at the given position.

		// Search the remaining chroms in order starting from last:
		List<String> chroms= Lists.reverse(this.getChromListStartingAt(startChrom));
		for(String chrom : chroms){
			prev = getPreviousFeatureOnChrom(chrom, Integer.MAX_VALUE);
			if(prev != null){
				return prev;
			}			
		}
		return null;
	}
	
	/** Experimental method 
	 * Find the previous feature relative to the given chrom-pos position.
	 * Read the chunk of file *before* the pos position and return the last feature found. I.e. find the feature immediately to 
	 * to left of pos and not overlapping pos.
	 * */
	private IntervalFeature getPreviousFeatureOnChrom(String chrom, int pos) throws IOException, InvalidGenomicCoordsException{
		
		int chunkFrom= Math.max(0, pos-1000000);
		int chunkTo= pos - 1; // -1 because we don't include the current position
		IntervalFeature last= null;

		while(chunkTo > 0){
			TabixBigBedIterator iter= this.getReader().query(chrom, chunkFrom, chunkTo);
			while(true){
				// Find the last feature in this chunk where end coordinate is less then pos
				String line= iter.next();
				if(line == null){
					break;
				} 
				line= line.trim();
				IntervalFeature candidate= new IntervalFeature(line, this.getTrackFormat());
				if(candidate.getTo() < pos && this.featureIsVisible(line)){ 
					// This is a candidate feature but we don't know yet of it's the last one
					last= new IntervalFeature(line, this.getTrackFormat());
				}
			}
			if(last != null){
				break; // The last feature is not null and valid. Stop looking 
			} else {
				// Move to previous chunk;
				chunkTo= chunkFrom; 
				chunkFrom= Math.max(0, chunkFrom-1000000);
			}
		}
		if(last == null || last.getTo() > pos || ! this.featureIsVisible(last.getRaw())){
			return null; 
		}
		return last;
	}
		
	protected GenomicCoords startEndOfNextFeature(GenomicCoords currentGc, boolean getPrevious) throws InvalidGenomicCoordsException, IOException {
		IntervalFeature nextFeature;
		if(getPrevious){
			nextFeature= getPreviousFeature(currentGc.getChrom(), currentGc.getFrom());
		} else {
			nextFeature= getNextFeature(currentGc.getChrom(), currentGc.getTo());
		}
		if(nextFeature == null){
			return currentGc;
		}
		GenomicCoords nextGc= new GenomicCoords(
				Utils.coordinatesToString(nextFeature.getChrom(), nextFeature.getFrom(), nextFeature.getTo()), currentGc.getSamSeqDict(),
				currentGc.getFastaFile());
		return nextGc;		
	}
	
	public GenomicCoords findNextMatch(GenomicCoords currentGc, String query) throws IOException, InvalidGenomicCoordsException{

		IntervalFeature nextFeature= findNextRegexInGenome(query, currentGc.getChrom(), currentGc.getTo());
		if(nextFeature == null){
			return currentGc;
		}
		GenomicCoords nextGc= new GenomicCoords(
				Utils.coordinatesToString(nextFeature.getChrom(), 
						                  nextFeature.getFrom(), 
						                  nextFeature.getFrom() + currentGc.getGenomicWindowSize() - 1),
				currentGc.getSamSeqDict(),
				currentGc.getFastaFile());
		return nextGc;
	}

	/** Execute findAllChromRegexInGenome() and return the extreme coordinates of the matched features */
	protected GenomicCoords genomicCoordsAllChromMatchInGenome(String query, GenomicCoords currentGc) throws IOException, InvalidGenomicCoordsException{

		List<IntervalFeature> matchedFeatures = this.findAllChromMatchInGenome(query, currentGc);
		
		if(matchedFeatures.size() == 0){
			return currentGc;
		}
		
		// Now get the coords of the first and last feature matched.
		String chrom= matchedFeatures.get(0).getChrom();
		int startFrom= matchedFeatures.get(0).getFrom();
		int endTo= matchedFeatures.get(matchedFeatures.size()-1).getTo();
		GenomicCoords allMatchesGc= new GenomicCoords(Utils.coordinatesToString(chrom, startFrom, endTo),
				currentGc.getSamSeqDict(),
				currentGc.getFastaFile());
		return allMatchesGc;
		
	}

	/** Find all the feature matching regex.
	 * Only the feature on one chromosome are returned and this chromsome is the first one to have a match.
	 * The search starts from the beginning of the current chrom and if nothing is found continues
	 * to the other chroms. 
	 * @throws InvalidGenomicCoordsException */
	private List<IntervalFeature> findAllChromMatchInGenome(String query, GenomicCoords currentGc) throws IOException, InvalidGenomicCoordsException{
		
		// Accumulate features here
		List<IntervalFeature> matchedFeatures= new ArrayList<IntervalFeature>(); 

		// We start search from input chrom
		List<String> chromSearchOrder= null;
		chromSearchOrder = getChromListStartingAt(currentGc.getChrom());
		
		chromSearchOrder.add(currentGc.getChrom());		
		for(String curChrom : chromSearchOrder){
		
			TabixBigBedIterator iter = this.getReader().query(curChrom , 0, Integer.MAX_VALUE);
			// Iterator iter = this.iteratorFromQuery(curChrom , 0, Integer.MAX_VALUE); // this.tabixReader.query(curChrom , 0, Integer.MAX_VALUE);
			while(true){
				String line= iter.next();
				if(line == null) break;
				boolean matched= Pattern.compile(query).matcher(line).find();
				if(matched){
					IntervalFeature x= new IntervalFeature(line, this.getTrackFormat());
					if(this.featureIsVisible(x.getRaw())){
						matchedFeatures.add(x);
					}
				}
			}
			if(matchedFeatures.size() > 0){
				// At least one feature matching regex found on this chrom.
				// Chech we are at the same position as the beginning. if so, continue to other chroms
				if(matchedFeatures.get(0).getChrom().equals(currentGc.getChrom()) && 
				   matchedFeatures.get(0).getFrom() == currentGc.getFrom() &&
				   matchedFeatures.get(matchedFeatures.size()-1).getTo() == currentGc.getTo()){
				   // Discard results and keep searching other chroms.
					matchedFeatures= new ArrayList<IntervalFeature>();
				} else {
					break;
				}
			}
		} // Loop chrom
		return matchedFeatures;
	}

	@Override
	public String printToScreen() throws InvalidGenomicCoordsException {
	
		List<String> printable= new ArrayList<String>();		
		int nLines= 0;
		try {
			for(List<IntervalFeature> listToPrint : this.stackFeatures()){
				
				nLines++;
				if(nLines > this.yMaxLines){
					// Limit the number of lines in output
					break;
				}
				printable.add(this.printToScreenOneLine(listToPrint));
			}
		} catch (Exception e) {
			e.printStackTrace();
		}

		return StringUtils.join(printable, "\n");
	}
	
	/** Return a string of a single line of (typically de-stacked) reads
	 * @throws IOException 
	 * @throws InvalidGenomicCoordsException 
	 * */
	private String printToScreenOneLine(List<IntervalFeature> listToPrint) throws InvalidGenomicCoordsException, IOException {
		
		int windowSize= this.getGc().getUserWindowSize();

		List<String> printable= new ArrayList<String>(); // Each item in this list occupies a character space in the terminal. 
		                                                 // NB: Each item is String not char because it might contain the ansi formatting.
		for(int i= 0; i < this.getGc().getMapping(windowSize).size(); i++){ // First create empty line
			printable.add(" ");
		}
		for(IntervalFeature intervalFeature : listToPrint){
			if(intervalFeature.getScreenFrom() == -1){
				throw new RuntimeException(); // Feature doesn't map to screen, this shouldn't happen
			}
			intervalFeature.setGtfAttributeForName( this.getGtfAttributeForName() );
			String[] text = intervalFeature.makeIdeogramFormatted(this.isNoFormat());
			
			int i= 0;
			for(int j= intervalFeature.getScreenFrom(); j <= intervalFeature.getScreenTo(); j++){
				printable.set(j, text[i]);
				i++;
			}
			
		}
		return StringUtils.join(printable, "");
	}

	protected String getUnformattedTitle(){

		if(this.isHideTitle()){
			return "";
		}
		
		String sq= "";
		if (this.getFeatureDisplayMode().equals(FeatureDisplayMode.COLLAPSED)){
			sq= "; collapsed";
		}
		String gapped= "";
		if(this.getGap() == 0){
			gapped= "; ungapped";
		}
		
		String grep= "";
		if( ! this.getShowRegex().equals(SHOW_REGEX) ){
			grep += " -i " + this.getShowRegex(); 
		}
		if( ! this.getHideRegex().equals(HIDE_REGEX) ){
			grep += " -e " + this.getHideRegex(); 
		}
		if( ! grep.isEmpty()){
			grep= "; grep" + grep; 
		}
		String title=  this.getTrackTag() + ";" 
	                 + " N: " + this.intervalFeatureList.size()
					 + grep
	                 + sq
	                 + gapped;
		return title;
	}
	
	@Override
	public String getTitle(){
		return this.formatTitle(this.getUnformattedTitle()) + "\n";
	}
	
	
//	/** Remove positional duplicates from list of interval features for more compact visualization. 
//	 * Squashing is done according to feature field which should be applicable to GTF/GFF only.*/
//	private List<IntervalFeature> squashFeatures(List<IntervalFeature> intervalList){
//
//		List<IntervalFeature> stack= new ArrayList<IntervalFeature>();
//		List<IntervalFeature> squashed= new ArrayList<IntervalFeature>();
//		for(IntervalFeature interval : intervalList){
//			if(stack.size() == 0 || stack.get(0).equalStranded(interval)){
//				// Accumulate features with same coords.
//				stack.add(interval);
//			} else {
//				squashed.add(stack.get(0));
//				stack.clear();
//				stack.add(interval);
//			}
//		}
//		if(stack.size() > 0){
//			squashed.add(stack.get(0));
//		}
//		return squashed;
//	}
	
	/**		
	 * Put in the same list reads that will go in the same line of text. 
	 * This method separates features touching or overlapping each other, useful for visualization.
	 * Each item of the output list going on its own line.
	 * @param space Space between text feature that touch each other on screen. Use 1 to have at least one space
	 * so that distinct features never look merged with adjacent ones. 0 will not put any space and will give more
	 * compact view.  
	 * See also TrackReads.stackReads();
	 * @throws InvalidGenomicCoordsException 
	 * @throws IOException 
	 */
	private List<List<IntervalFeature>> stackFeatures() throws InvalidGenomicCoordsException, IOException{
		
		List<IntervalFeature> intervals; 
		List<IntervalFeature> flatListOfTx = this.flatListOfPrintableFeatures();
		if(this.getFeatureDisplayMode().equals(FeatureDisplayMode.COLLAPSED)) {
			intervals = Utils.mergeIntervalFeatures(flatListOfTx, false);
		} 
		else if(this.getFeatureDisplayMode().equals(FeatureDisplayMode.ONELINE)) {
			intervals = Utils.mergeIntervalFeatures(flatListOfTx, true);
		} 
		else {
			intervals = flatListOfTx;
		}
		
		// Make a copy of the IntervalFeature list. Items will be popped out as they are 
		// added to individual lines. 
		List<IntervalFeature> flatList= new ArrayList<IntervalFeature>(intervals);
				
		List<List<IntervalFeature>> listOfLines= new ArrayList<List<IntervalFeature>>();
		if(flatList.size() == 0){
			return listOfLines;
		}
		List<IntervalFeature> line= new ArrayList<IntervalFeature>();
		line.add(flatList.get(0)); 
		flatList.remove(0);
		listOfLines.add(line);

		while(true){
			ArrayList<IntervalFeature> trToRemove= new ArrayList<IntervalFeature>();
			// Find a read in input whose start is greater then end of current
			for(int i=0; i < flatList.size(); i++){
				IntervalFeature intervalFeature= flatList.get(i);
				// int gap= 1; // Add a space between book-end features
				if(intervalFeature.getScreenFrom() > line.get(line.size()-1).getScreenTo()+this.getGap()){ // +2 because we want some space between adjacent reads
					listOfLines.get(listOfLines.size()-1).add(intervalFeature); // Append to the last line. 
					trToRemove.add(intervalFeature);
				}
			} // At the end of the loop you have put in line as many reads as you can. 
			for(IntervalFeature intervalFeature : trToRemove){ 
				flatList.remove(flatList.indexOf(intervalFeature));
			}
			// Create a new line, add the first IntervalFeature in list
			if(flatList.size() > 0){
				line= new ArrayList<IntervalFeature>();
				line.add(flatList.get(0));
				listOfLines.add(line);
				flatList.remove(0);
			} else {
				break;
			}
		}
		return listOfLines;
	}

	/**Print raw features under track. 
	 * windowSize size the number of characters before clipping occurs. This is 
	 * typically the window size for plotting. windowSize is used only by CLIP mode.  
	 * @throws IOException 
	 * @throws InvalidGenomicCoordsException 
	 * */
	private String printFeatures() throws InvalidGenomicCoordsException, IOException{

		int windowSize= this.getGc().getUserWindowSize();
		if(this.getPrintMode().equals(PrintRawLine.FULL)){
			windowSize= Integer.MAX_VALUE;
		} else if(this.getPrintMode().equals(PrintRawLine.CLIP)){
			// Keep windowSize as it is
		} else {
			return "";
		} 
		
		List<String> featureList= new ArrayList<String>();
		
		int count= this.getPrintRawLineCount();
		for(IntervalFeature ift : intervalFeatureList){
			featureList.add(ift.getRaw());
			count--;
			if(count == 0){
				int omitted= intervalFeatureList.size() - this.getPrintRawLineCount();
				if(omitted > 0){
					System.err.println("[" + omitted + "/"  + intervalFeatureList.size() + " features omitted]");
				}
				break;
			}
		}
		List<String> tabList= Utils.tabulateList(featureList);
		StringBuilder sb= new StringBuilder();
		for(String x : tabList){
			if(x.length() > windowSize){
				x= x.substring(0, windowSize);
			}			
			sb.append(x + "\n");
		}
		return sb.toString(); // NB: Leave last trailing /n
	}

	@Override
	/** Write the features in interval to file by appending to existing file. 
	 * If the file to write to null or empty, return the data that would be
	 * written as string.
	 * printFeaturesToFile aims at reproducing the behavior of Linux cat: print to file, possibly appending or to stdout. 
	 * */
	public String printFeaturesToFile() throws IOException, InvalidGenomicCoordsException {
		
		if(this.getExportFile() == null || this.getExportFile().isEmpty()){
			return this.printFeatures();
		}
		
		BufferedWriter wr= null;
		try{
			wr = new BufferedWriter(new FileWriter(this.getExportFile(), true));
			for(IntervalFeature ift : intervalFeatureList){
				wr.write(ift.getRaw() + "\n");
			}
			wr.close();
		} catch(IOException e){
			System.err.println("Cannot write to " + this.getExportFile());
			throw e;
		}
		return "";
	}
	
	/** Searching the current chrom starting at "from" to find the *next* feature matching the given string. 
	 * If not found, search the other chroms, if not found restart from the beginning of
	 * the current chrom until the "from" position is reached. 
	 * @throws InvalidGenomicCoordsException */
	protected IntervalFeature findNextRegexInGenome(String query, String chrom, int from) throws IOException, InvalidGenomicCoordsException{
		
		int startingPoint= from-1; // -1 because tabix.query from is 0 based (seems so at least)
		List<String> chromSearchOrder = this.getChromListStartingAt(chrom);
		chromSearchOrder.add(chrom);
		for(String curChrom : chromSearchOrder){
			
			TabixBigBedIterator iter= this.getReader().query(curChrom , startingPoint, Integer.MAX_VALUE);
			// Iterator iter = this.iteratorFromQuery(curChrom , startingPoint, Integer.MAX_VALUE); // this.getTabixReader().query(curChrom , startingPoint, Integer.MAX_VALUE);
			while(true){
				String line= iter.next();
				if(line == null) break;
				boolean matched= Pattern.compile(query).matcher(line).find();
				if(matched){
					IntervalFeature x= new IntervalFeature(line, this.getTrackFormat());
					if(x.getFrom() > startingPoint && this.featureIsVisible(x.getRaw())){
						return x;
					}
				} 
			}
			startingPoint= 0;
		} return null; // Not found anywhere
	}

	/** Return the set chroms sorted but with and first chrom set to startChrom.
	 * 	chroms:         chr1 chr2 chr3 chr4 chr5 chr6 chr7
	 *	startChrom:     chr3
	 *  return:         chr3 chr4 chr5 chr6 chr7 chr1 chr2
	 *  */
	private List<String> getChromListStartingAt(String startChrom){
		
		List<String> chroms= this.getChromosomeNames();
		
		int idx= chroms.indexOf(startChrom);
		if(idx == -1){ // If startChrom is not present at all in the bed/gtf file.
			return chroms;
		}
		List<String> chromsStartingAt= new ArrayList<String>();
		chromsStartingAt.addAll(chroms.subList(idx, chroms.size()));
		chromsStartingAt.addAll(chroms.subList(0, idx));
		
		// Sanity check
		if(chroms.size() != chromsStartingAt.size()){ 
			throw new RuntimeException("Error reordering chroms. Expected " + chroms.size() + " chroms got " + chromsStartingAt.size());
		}
		if(! (chromsStartingAt.containsAll(chroms) && chroms.containsAll(chromsStartingAt))){
			throw new RuntimeException("Error re-ordering chromosomes");
		}
		return chromsStartingAt;
	}

	/** Group the features in this genomic window by GFF attribute (typically a transcripts). 
	 * Features that don't have the attribute make each a length=1 list.
	 * */
	private Map<String, List<IntervalFeature>> groupByGFFAttribute(){

		// * First collect the IDs of the transcripts
				
		// Key is transcript ID e.g. ENST00001234. 
		// Values is all the IntervalFeatures captured by this ID and part of a transcript.
		// I.e. their are in txFeature set, 
		Map<String, List<IntervalFeature>> txIds= new LinkedHashMap<String, List<IntervalFeature>>();
		
		// This key:value is for records which are not part transcripts. E.g. features like "chromosome" or rRNA.
		txIds.put("_na_", new ArrayList<IntervalFeature>()); 

		// Now populate the lists of values by assigning to each key the transcript records:
		for(IntervalFeature x : this.getIntervalFeatureList()){
			
			if(FormatGTF.getTxSuperFeatures().contains(x.getFeature().toLowerCase())){
				// Transcript feature. E.g.
				// 7 ensembl_havana mRNA 5527151 5530709 . - . ID=transcript:ENST00000331789;Parent=gene:ENSG00000075624;Name=ACTB-...
				String txId= x.getGFFValueFromKey("ID");
				if( ! txIds.containsKey(txId)){
					txIds.put(txId, new ArrayList<IntervalFeature>());
				}
				txIds.get(txId).add(x);
			} else if(FormatGTF.getTxSubFeatures().contains(x.getFeature().toLowerCase())){
				// Part of transcript, e.g:
				// 7 ensembl_havana exon 5527151 5527891 . - . Parent=transcript:ENST00000331789;Name=ENSE00001902654;constitutive=0;ensembl_end_pha
				String txId= x.getGFFValueFromKey("Parent");
				if( ! txIds.containsKey(txId)){
					txIds.put(txId, new ArrayList<IntervalFeature>());
				}
				txIds.get(txId).add(x);				
			} else {
				// Not a transcript or part thereof. E.g.
				// 7 . biological_region 5529708 5529709 0.999 - . logic_name=eponine
				txIds.get("_na_").add(x);
			}
		}
		// We don't need to return the full Map, only the list of lists (groups) would suffice.
		// However, we need to separate the group of non-trascripts (_na_ key)
		return txIds;	
	}

	
	/** Group the features in this genomic window by GTF attribute (typically a transcripts). 
	 * */
	private Map<String, List<IntervalFeature>> groupByGTFAttribute(){

		// * First collect the IDs of the transcripts
				
		// Key is transcript ID e.g. ENST00001234. 
		// Values is all the IntervalFeatures captured by this ID and part of a transcript.
		// I.e. their are in txFeature set, 
		Map<String, List<IntervalFeature>> txIds= new LinkedHashMap<String, List<IntervalFeature>>();
		
		// This key:value is for records which are not part transcripts. E.g. features like "chromosome" or rRNA.
		txIds.put("_na_", new ArrayList<IntervalFeature>()); 

		// Now populate the lists of values by assigning to each key the transcript records:
		for(IntervalFeature x : this.getIntervalFeatureList()){
			
			if(FormatGTF.getTxSuperFeatures().contains(x.getFeature().toLowerCase()) || FormatGTF.getTxSubFeatures().contains(x.getFeature().toLowerCase())){
				// Transcript feature. E.g.
				// chr7 hg19_wgEncodeGencodeBasicV19 exon       5566782 5567522 0.000000 - . gene_id "ENST00000331789.5"; transcript_id "ENST00000331789.5";
				String txId= x.getGFFValueFromKey("transcript_id");
				if( ! txIds.containsKey(txId)){
					txIds.put(txId, new ArrayList<IntervalFeature>());
				}
				txIds.get(txId).add(x);				
			} else {
				// Not a transcript or part thereof. E.g.
				// 7 . biological_region 5529708 5529709 0.999 - . logic_name=eponine
				txIds.get("_na_").add(x);
			}
		}
		// We don't need to return the full Map, only the list of lists (groups) would suffice.
		// However, we need to separate the group of non-trascripts (_na_ key)
		return txIds;	
	}
	
	/** Collapse the list of features in a single IntervalFeature representing the transcript.
	 * The elements of txFeatures are expected to represent the entire transcript, nothing more
	 * (e.g. "chromosome"). The transcript may not be biologically complete as part of it
	 * may be outside the current genomic coords. 
	 *
	 * mapToScreen: Mapping of genomic coordinates to screen coordinates. This could be obtained inside this
	 * method but better to pass it from outside as it can take time to get the terminal window size several times.  
	 * 
	 * @throws InvalidGenomicCoordsException 
	 * */
	private IntervalFeature collapseGFFTranscript(List<IntervalFeature> txFeatures, List<Double> mapToScreen) throws InvalidGenomicCoordsException{
		
		if(txFeatures.size() == 0){
			System.err.println("Unexpected transcript: Length zero!");
			throw new RuntimeException();
		}
		
		// Collect the genomic and screen coordinates of this transcript
		int gFrom= Integer.MAX_VALUE;
		int gTo= 0;
		int screenFrom= Integer.MAX_VALUE;
		int screenTo= 0;
		for(IntervalFeature x : txFeatures){

			if(x.getFrom() < gFrom){
				gFrom= x.getFrom(); 
			}
			if(x.getTo() > gTo){
				gTo= x.getTo(); 
			}
			if(x.getScreenFrom() < screenFrom){
				screenFrom= x.getScreenFrom(); 
			}
			if(x.getScreenTo() > screenTo){
				screenTo= x.getScreenTo();
			}
			x.setGtfAttributeForName( this.getGtfAttributeForName() );
		}
				
		IntervalFeature transcript= new IntervalFeature(txFeatures.get(0).getChrom(), gFrom, gTo, TrackFormat.GFF); 
		transcript.setStrand(txFeatures.get(0).getStrand());
		transcript.mapToScreen(mapToScreen);
		
		// Now we need to prepare the ideogram
		char[] ideogram= new char[screenTo - screenFrom + 1];
		for(int i= 0; i < ideogram.length; i++){
			ideogram[i]= '-'; // Deafult characater to print. Typically this should apply to introns only.
		}
		HashMap<String, Character> charDict = FormatGTF.getFeatureToTextCharDict().get(transcript.getStrand());

		for(String txSubType : FormatGTF.getTxSubFeatures()){

			for(IntervalFeature subFeature : txFeatures){
				
				if(subFeature.getFeature().toLowerCase().equals(txSubType)){

					char c= charDict.get(subFeature.getFeature().toLowerCase());
					
					// Fill up text with this character
					for(int i= subFeature.getScreenFrom(); i <= subFeature.getScreenTo(); i++){
						ideogram[i - screenFrom]= c;
					}
				}
			}			
		}
		transcript.setIdeogram(ideogram);
		
		//Now we get the name for this transcript
		String txName= "."; // Default: No name
		for(String txSuperType : FormatGTF.getTxSuperFeatures()){
			for(IntervalFeature x : txFeatures){
				if(x.getFeature().toLowerCase().equals(txSuperType)){
					txName= x.getName();
				}
				if( txName != null  && ! txName.isEmpty() && ! txName.equals(".")){ break; } // A name found, break
			}
			if( txName != null  && ! txName.isEmpty() && ! txName.equals(".")){ break; }
		}
		transcript.setName(txName);
		
		// ideogram= transcript.addNameToIdeogram(ideogram);
		// transcript.setIdeogram(ideogram);
		return transcript;
	}
	
	/** List where the original records have been grouped into transcripts, if there are 
	 * transcripts. 
	 * TODO: Process here also squash, merge and gap?
	 * @throws IOException 
	 * @throws InvalidGenomicCoordsException 
	 * */
	private List<IntervalFeature> flatListOfPrintableFeatures() throws InvalidGenomicCoordsException, IOException{
		
		List<IntervalFeature> flatList= new ArrayList<IntervalFeature>(); 

		List<Double> mapToScreen = this.getGc().getMapping(this.getGc().getUserWindowSize());
		
		if(this.getTrackFormat().equals(TrackFormat.GFF) || this.getTrackFormat().equals(TrackFormat.GTF)){
		
			Map<String, List<IntervalFeature>> tx;
			if(this.getTrackFormat().equals(TrackFormat.GFF)){
				tx = this.groupByGFFAttribute();
			} else if (this.getTrackFormat().equals(TrackFormat.GTF)) {
				tx = this.groupByGTFAttribute();
			} else {
				throw new RuntimeException("This is not a GTF or GFF track!");
			}
			
			for(String txId : tx.keySet()){
				if(txId.equals("_na_")){
					flatList.addAll(tx.get(txId));
				} else {
					flatList.add(this.collapseGFFTranscript(tx.get(txId), mapToScreen));
				}
			}
			
		} else {
			flatList.addAll(this.getIntervalFeatureList());
		}
		Collections.sort(flatList);
		return flatList;
	}
	

	/** Merge feature using the screen coordinates to establish overlap instead of the genomic coordinates
	 * */
//	private List<IntervalFeature> mergeIntervalFeaturesOnScreen(List<IntervalFeature> intervalList) throws InvalidGenomicCoordsException{
//		List<IntervalFeature> mergedList= new ArrayList<IntervalFeature>();		 
//		if(intervalList.size() == 0){
//			return mergedList;
//		}
//		
//		String chrom = null;
//		int from= -1;
//		int to= -1;
//		int screenFrom= -1;
//		int screenTo= -1;
//		int numMrgIntv= 1; // Number of intervals in the merged one. 
//		Set<Character> strand= new HashSet<Character>(); // Put here all the different strands found in the merged features.		
//		
//		for(int i= 0; i < (intervalList.size()+1); i++){
//			// We do an additional loop to add to the mergedList the last interval.
//			// The last loop has interval == null so below you need to account for it
//			IntervalFeature interval= null;
//			if(i < intervalList.size()){
//				interval= intervalList.get(i); 
//			}
//			
//			if(from < 0){ // Init interval
//				chrom= interval.getChrom(); 
//				from= interval.getFrom();
//				to= interval.getTo();
//				screenFrom= interval.getScreenFrom();
//				screenTo= interval.getScreenTo();
//				continue;
//			}
//			// Sanity check: The list to be merged is on the same chrom and sorted by start pos.
//			if(i < intervalList.size() && (!chrom.equals(interval.getChrom()) || from > interval.getFrom() || from > to)){
//				System.err.println(chrom + " " + from + " " + to);
//				throw new RuntimeException();
//			} 
//			if(i < intervalList.size() && (screenFrom <= interval.getScreenTo() && screenTo >= (interval.getScreenFrom()-1) )){ 
//				// Overlap: Extend <to> coordinate. See also http://stackoverflow.com/questions/325933/determine-whether-two-date-ranges-overlap
//				to= interval.getTo();
//				screenTo= interval.getScreenTo();
//				strand.add(interval.getStrand());
//				numMrgIntv++;
//			} else {
//				// No overlap add merged interval to list and reset new merged interval
//				IntervalFeature x= new IntervalFeature(chrom + "\t" + (from-1) + "\t" + to, TrackFormat.BED);
//				x.setScreenFrom(screenFrom);
//				x.setScreenTo(screenTo);
//				if(strand.size() == 1){
//					x.setStrand(strand.iterator().next());
//				} 
//				strand.clear();
//				
//				if(x.equals(intervalList.get(i-1)) && numMrgIntv == 1){
//					mergedList.add(intervalList.get(i-1));
//				} else {
//					mergedList.add(x);
//				}
//				
//				if(i < intervalList.size()){
//					// Do not reset from/to if you are in extra loop.
//					from= interval.getFrom();
//					to= interval.getTo();
//					screenFrom= interval.getScreenFrom();
//					screenTo= interval.getScreenTo();
//					strand.add(interval.getStrand());
//					numMrgIntv= 1;
//				}
//			}
//		}
//		return mergedList;
//	}
	
	// SETTERS AND GETTERS
	// -------------------
	
//	private TrackFormat getType() {
//		return this.type;
//	}

	protected List<IntervalFeature> getIntervalFeatureList() {
		return intervalFeatureList;
	}

	protected void setIntervalFeatureList(List<IntervalFeature> intervalFeatureList) {
		this.intervalFeatureList = intervalFeatureList;
	}

	@Override
	public void setHideRegex(String hideRegex) throws ClassNotFoundException, IOException, InvalidGenomicCoordsException, InvalidRecordException, SQLException {
		this.hideRegex= hideRegex;
		this.update();
	}
	@Override
	public String getHideRegex() {
		return this.hideRegex;
	}

	@Override
	public void setShowRegex(String showRegex) throws ClassNotFoundException, IOException, InvalidGenomicCoordsException, InvalidRecordException, SQLException {
		this.showRegex= showRegex;
		this.update();
	}
	@Override
	public String getShowRegex() {
		return this.showRegex;
	}

	@Override
	public List<String> getChromosomeNames(){
		ArrayList<String> x = new ArrayList<String>(this.getReader().getChromosomes());
		Collections.sort(x);
		return x;
	}
	
	private TabixBigBedReader getReader(){
		
		if(this.bigBedReader != null && this.tabixReader != null){
			System.err.println("You cannot have tabix and bigBed readers bot set!");
			throw new RuntimeException();
		}
	
		if(this.bigBedReader != null){
			return new TabixBigBedReader(this.bigBedReader);
		} else if(this.tabixReader != null){
			return new TabixBigBedReader(this.tabixReader);
		} else {
			System.err.println("Tabix and bigBed reader both null.");
			throw new RuntimeException();
		}
	}

	@Override
	protected void setPrintRawLineCount(int count) {
		if(count <= 0){
			count= Integer.MAX_VALUE;
		}
		this.printRawLineCount= count;
	}

	protected int getPrintRawLineCount() {
		return this.printRawLineCount;
	}
	
	/** This setter is for TrackBookmark to work.*/
	protected void setTabixReader(TabixReader tabixReader) {
		this.tabixReader = tabixReader;
	}
	protected TabixReader getTabixReader() {
		return this.tabixReader;
	}
}
