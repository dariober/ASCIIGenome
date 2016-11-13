package tracks;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Set;
import java.util.regex.Pattern;
import org.apache.commons.lang3.StringUtils;

import exceptions.InvalidCommandLineException;
import exceptions.InvalidGenomicCoordsException;
import exceptions.InvalidRecordException;
import htsjdk.tribble.readers.TabixReader;
import htsjdk.tribble.readers.TabixReader.Iterator;
import samTextViewer.GenomicCoords;
import samTextViewer.Utils;
import sortBgzipIndex.MakeTabixIndex;
import ucsc.UcscGenePred;

public class TrackIntervalFeature extends Track {
 
	protected List<IntervalFeature> intervalFeatureList= new ArrayList<IntervalFeature>();  
	/**For GTF/GFF data: Use this attribute to get the feature names 
	 * */
	private String hideRegex= "^$";
	private String showRegex= ".*";
	private TrackFormat type;
	protected TabixReader tabixReader;
	
	/* C o n s t r u c t o r */

	public TrackIntervalFeature(final String filename, GenomicCoords gc) throws IOException, InvalidGenomicCoordsException, ClassNotFoundException, InvalidRecordException, SQLException{
		
		// String sourceFile= filename; // sourceFile is what is actually used to construct the tabix file. 
		this.setWorkFilename(filename);
		
		if(Utils.isUcscGenePredSource(filename)){
			UcscGenePred ucsc = null;
			try {
				ucsc = new UcscGenePred(filename, -1);
				this.setWorkFilename(ucsc.getTabixFile());
				this.type= TrackFormat.GFF;
			} catch (InvalidCommandLineException e) {
				//
			} 
		} else {
			this.type= Utils.getFileTypeFromName(new File(filename).getName());			
			if( ! Utils.hasTabixIndex(new File(filename).getAbsolutePath())){
				// Tabix index not found for this file. Sort and index input to tmp.
	
				String suffix= new File(filename).getName();
				if( ! suffix.endsWith(".gz")){
					suffix += ".gz";
				}
				String tmpWorkFile= File.createTempFile("asciigenome.", "." + suffix).getAbsolutePath();
				new File(tmpWorkFile).deleteOnExit();
				new File(new File(tmpWorkFile).getAbsolutePath() + ".tbi").deleteOnExit();
				this.setWorkFilename(tmpWorkFile);	
				
				new MakeTabixIndex(filename, new File( this.getWorkFilename() ), Utils.trackFormatToTabixFormat(this.type));	
			} 
		}
		this.tabixReader= new TabixReader(new File(this.getWorkFilename()).getAbsolutePath());
		this.setGc(gc);
		this.setFilename(filename);		
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
			System.err.println("Invalid coordinates: from: " + from + ";to: " + to 
					+ "Resetting to initial 1-" + Integer.MAX_VALUE);
			from= 1;
			to= Integer.MAX_VALUE;
			throw new InvalidGenomicCoordsException();
		}		
		List<IntervalFeature> xFeatures= new ArrayList<IntervalFeature>();
		//if(isTabix){
		Iterator qry = this.tabixReader.query(chrom,  from-1, to);
		while(true){
			String q = qry.next();
			if(q == null){
				break;
			}
			IntervalFeature intervalFeature= new IntervalFeature(q, this.type);
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
	public GenomicCoords coordsOfNextFeature(GenomicCoords currentGc) throws InvalidGenomicCoordsException, IOException {
		IntervalFeature nextFeature= this.getNextFeatureOnChrom(currentGc.getChrom(), currentGc.getTo());
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

	/** Get the next feature on chrom after "from" position or null if no 
	 * feature found 
	 * @throws IOException 
	 * @throws InvalidGenomicCoordsException */
	private IntervalFeature getNextFeatureOnChrom(String chrom, int from) throws IOException, InvalidGenomicCoordsException{
		
		Iterator iter = this.tabixReader.query(chrom, from-1, Integer.MAX_VALUE);
		while(true){
			String line= iter.next();
			if(line == null){
				return null;
			} 
			IntervalFeature x= new IntervalFeature(line, this.type);
			if(x.getFrom() > from && this.featureIsVisible(x.getRaw())){
				return x;
			}
		}
	}

	protected GenomicCoords startEndOfNextFeature(GenomicCoords currentGc) throws InvalidGenomicCoordsException, IOException {
		IntervalFeature nextFeature= getNextFeatureOnChrom(currentGc.getChrom(), currentGc.getTo());
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
		chromSearchOrder = getChromListStartingAt(this.tabixReader.getChromosomes(), currentGc.getChrom());
		
		chromSearchOrder.add(currentGc.getChrom());		
		for(String curChrom : chromSearchOrder){
		
			Iterator iter = this.tabixReader.query(curChrom , 0, Integer.MAX_VALUE);
			while(true){
				String line= iter.next();
				if(line == null) break;
				boolean matched= Pattern.compile(query).matcher(line).find();
				if(matched){
					IntervalFeature x= new IntervalFeature(line, this.type);
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
		List<String> printable= new ArrayList<String>();
		for(int i= 0; i < this.getGc().getMapping(windowSize).size(); i++){ // First create empty track
			printable.add(" ");
		}
		for(IntervalFeature intervalFeature : listToPrint){
			if(intervalFeature.getScreenFrom() == -1){
				continue; // Feature doesn't map to screen, this shouldn't happen though
			}
			intervalFeature.setGtfAttributeForName(this.getGtfAttributeForName());
			String nameOnFeature= intervalFeature.getName().trim() + "_";
			int relPos= 0;
			for(int j= intervalFeature.getScreenFrom(); j <= intervalFeature.getScreenTo(); j++){
				
				// Default is to use the feature type as printable text, e.g. 'E', 'T', '>', '<', '|', etc.
				String text= intervalFeature.assignTextToFeature(this.isNoFormat()); 
 
				if((intervalFeature.getScreenTo() - intervalFeature.getScreenFrom() + 1) > 4
						// (j - intervalFeature.getScreenFrom()) > 0 // First char is feature type (E, C, etc.)
						&& j < intervalFeature.getScreenTo() // Last char is feature type (E, C, etc.)
						&& relPos < nameOnFeature.length() 
						&& !nameOnFeature.equals("._")){
					// If these conds are satisfied, use the chars in the name as printable chars. 
					Character x= nameOnFeature.charAt(relPos); 
					if(this.isNoFormat()){
						text= Character.toString(x); 	
					} else {
						text= FormatGTF.format(x, intervalFeature.getStrand());
					}
					relPos += 1;
				}

				printable.set(j, text);
			}
		}
		return StringUtils.join(printable, "");
	}

	protected String getUnformattedTitle(){

		if(this.isHideTitle()){
			return "";
		}
		
		String sq= "";
		if(this.getFeatureDisplayMode().equals(FeatureDisplayMode.SQUASHED)){
			sq= "; squashed";
		} else if (this.getFeatureDisplayMode().equals(FeatureDisplayMode.MERGED)){
			sq= "; merged";
		}
		String gapped= "";
		if(this.getGap() == 0){
			gapped= "; ungapped";
		}
		String title=  this.getTrackTag() + "; " 
	                 + "Incl " + this.getShowRegex()
	                 + " Excl " + this.getHideRegex()
	                 + " N: " + this.intervalFeatureList.size()
	                 + sq
	                 + gapped;
		return title;
	}
	
	@Override
	public String getTitle(){
		return this.formatTitle(this.getUnformattedTitle()) + "\n";
	}
	
	
	/** Remove positional duplicates from list of interval features for more compact visualization. 
	 * Squashing is done according to feature field which should be applicable to GTF/GFF only.*/
	private List<IntervalFeature> squashFeatures(List<IntervalFeature> intervalList){

		List<IntervalFeature> stack= new ArrayList<IntervalFeature>();
		List<IntervalFeature> squashed= new ArrayList<IntervalFeature>();
		for(IntervalFeature interval : intervalList){
			if(stack.size() == 0 || stack.get(0).equalStranded(interval)){
				// Accumulate features with same coords.
				stack.add(interval);
			} else {
				squashed.add(stack.get(0));
				stack.clear();
				stack.add(interval);
			}
		}
		if(stack.size() > 0){
			squashed.add(stack.get(0));
		}
		return squashed;
	}
	
	/**		
	 * Put in the same list reads that will go in the same line of text. 
	 * This method separates features touching or overlapping each other, useful for visualization.
	 * Each item of the output list going on its own line.
	 * @param space Space between text feature that touch each other on screen. Use 1 to have at least one space
	 * so that distinct features never look merged with adjacent ones. 0 will not put any space and will give more
	 * compact view.  
	 * See also TrackReads.stackReads();
	 * @throws InvalidGenomicCoordsException 
	 */
	private List<List<IntervalFeature>> stackFeatures() throws InvalidGenomicCoordsException{
		
		List<IntervalFeature> intervals; 
		if(this.getFeatureDisplayMode().equals(FeatureDisplayMode.SQUASHED)){
			intervals = this.squashFeatures(this.intervalFeatureList);
		} else if(this.getFeatureDisplayMode().equals(FeatureDisplayMode.MERGED)) {
			intervals = Utils.mergeIntervalFeatures(this.intervalFeatureList);
		} else {
			intervals = this.intervalFeatureList;
		}
		
		// Make a copy of the IntervalFeature list. Items will be popped out as they are 
		// added to individual lines. 
		List<IntervalFeature> flatList= new ArrayList<IntervalFeature>();
		for(IntervalFeature x : intervals){ // this.intervalFeatureList
			flatList.add(x);
		}
				
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
		for(IntervalFeature ift : intervalFeatureList){
			featureList.add(ift.getRaw());
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
	/** Write the features in interval to file. if append is true append to existing file. 
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
			wr = new BufferedWriter(new FileWriter(this.getExportFile(), this.isAppendToExportFile()));
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
		List<String> chromSearchOrder = this.getChromListStartingAt(this.getTabixReader().getChromosomes(), chrom);
		chromSearchOrder.add(chrom);
		for(String curChrom : chromSearchOrder){
			
			Iterator iter = this.getTabixReader().query(curChrom , startingPoint, Integer.MAX_VALUE);
			while(true){
				String line= iter.next();
				if(line == null) break;
				boolean matched= Pattern.compile(query).matcher(line).find();
				if(matched){
					IntervalFeature x= new IntervalFeature(line, this.getType());
					if(x.getFrom() > startingPoint && this.featureIsVisible(x.getRaw())){
						return x;
					}
				} 
			}
			startingPoint= 0;
		} return null; // Not found anywhere
	}

	/** Return the set chroms sorted but with and first chrom set to startChrom.
	 *  */
	protected List<String> getChromListStartingAt(Set<String> chroms, String startChrom){
		// Set:            chr1 chr2 chr3 chr4 chr5 chr6 chr7
		// StartChrom:     chr3
		// Ordered chroms: chr3 chr4 chr5 chr6 chr7 chr1 chr2
		List<String> orderedChroms= new ArrayList<String>();
		orderedChroms.addAll(chroms);
		Collections.sort(orderedChroms);
		int idx= orderedChroms.indexOf(startChrom);
		if(idx == -1){ // If startChrom is not present at all in the bed/gtf file.
			return orderedChroms;
		}
		List<String> chromsStartingAt= new ArrayList<String>();
		chromsStartingAt.addAll(orderedChroms.subList(idx, orderedChroms.size()));
		chromsStartingAt.addAll(orderedChroms.subList(0, idx));
		return chromsStartingAt;
	}

	// SETTERS AND GETTERS
	// -------------------
	
	private TrackFormat getType() {
		return this.type;
	}

	private TabixReader getTabixReader() {
		return this.tabixReader;
	}

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

	protected void setType(TrackFormat type) {
		this.type = type;
	}

}
