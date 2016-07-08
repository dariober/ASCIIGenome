package tracks;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.regex.Pattern;
import java.util.regex.PatternSyntaxException;

import org.apache.commons.lang3.text.StrMatcher;
import org.apache.commons.lang3.text.StrTokenizer;

import exceptions.InvalidCommandLineException;
import exceptions.InvalidGenomicCoordsException;
import samTextViewer.GenomicCoords;
import samTextViewer.Utils;

/** Class to hold tracks to be printed. 
 * */
public class TrackSet {
	
	private LinkedHashMap<String, Track> trackSet= new LinkedHashMap<String, Track>();
	
	/*   C o n s t r u c t o r s   */
	
	public TrackSet(){}
	
	/*   M e t h o d s   */

	//public void addOrReplace(Track track){
	//	this.trackSet.put(track.getFileTag(), track);
	//}

	/** From cmdInput extract regex and yMaxLines then iterate through the tracks list to set 
	 * the yMaxLines in the tracks whose filename matches the regex.
	 * The input list is updated in place! 
	*/
	public void setTrackHeightForRegex(String cmdInput) throws InvalidCommandLineException{

		// MEMO of subcommand syntax:
		// 0 trackHeight
		// 1 int    mandatory
		// 2 regex  optional
		
		StrTokenizer str= new StrTokenizer(cmdInput);
		str.setTrimmerMatcher(StrMatcher.spaceMatcher());
		str.setQuoteChar('\'');
		List<String> tokens= str.getTokenList();
		if(tokens.size() < 2){
			System.err.println("Error in trackHeight subcommand. Expected 2 args got: " + cmdInput);
			throw new InvalidCommandLineException();
		}
		String trackNameRegex= ".*"; // Default: Capture everything
		if(tokens.size() == 3){ // If size 3 (trackHeight int regex) user has set a regex. Used that instead of default.
			trackNameRegex= tokens.get(2);
		}
		
		try{
			Pattern.compile(trackNameRegex); // Validate regex
		} catch(PatternSyntaxException e){
	    	System.err.println("Invalid regex in: " + cmdInput);
	    	System.err.println(e.getDescription());
		}
		
		int trackHeight= 0;
		try{
			trackHeight= Integer.parseInt(tokens.get(1));
			trackHeight= trackHeight < 0 ? 0 : trackHeight;
		} catch(NumberFormatException e){
			System.err.println("Number format exception: " + trackHeight);
		}
		for(Track tr : this.trackSet.values()){
			boolean matched= Pattern.compile(trackNameRegex).matcher(tr.getFileTag()).find();
			if(matched){
				tr.setyMaxLines(trackHeight);
			}
		}
	}
	
	public void setBisulfiteModeForRegex(String cmdInput) {

		// MEMO of subcommand syntax:
		// 0 BSseq
		// 1 Regex

		StrTokenizer str= new StrTokenizer(cmdInput);
		str.setTrimmerMatcher(StrMatcher.spaceMatcher());
		str.setQuoteChar('\'');
		List<String> tokens= str.getTokenList();

		// Regex
		String trackNameRegex= "^$"; // Default: Capture nothing
		if(tokens.size() >= 2){
			trackNameRegex= tokens.get(1);
		}
		try{
			Pattern.compile(trackNameRegex); // Validate regex
		} catch(PatternSyntaxException e){
	    	System.err.println("Invalid regex in: " + cmdInput);
	    	System.err.println(e.getDescription());
		}
		
		for(Track tr : this.trackSet.values()){
			boolean matched= Pattern.compile(trackNameRegex).matcher(tr.getFileTag()).find();
			if(matched){
				System.out.println("Setting " + tr);
				if(tr.isBisulf()){ // Invert setting
					tr.setBisulf(false);
				} else {
					tr.setBisulf(true);
				}
			}
		}
	}

	
	public void setTrackColourForRegex(String cmdInput) throws InvalidCommandLineException{

		// MEMO of subcommand syntax:
		// 0 trackColour
		// 1 Colour
		// 2 Regex

		StrTokenizer str= new StrTokenizer(cmdInput);
		str.setTrimmerMatcher(StrMatcher.spaceMatcher());
		str.setQuoteChar('\'');
		List<String> tokens= str.getTokenList();

		// Colour
		String colour= (new Track()).getTitleColour();
		if(tokens.size() >= 2){
			String xcolour= tokens.get(1).toLowerCase();
			if(!Utils.ansiColourCodes().containsKey(xcolour)){
				System.err.println("\nGot invalid colour: " + xcolour + ". Resetting to " + colour);
				System.err.println("Valid colours are: " + Utils.ansiColourCodes().keySet());
			} else {
				colour= xcolour;
			}
		}
		
		// Regex
		String trackNameRegex= ".*"; // Default: Capture everything
		if(tokens.size() >= 3){
			trackNameRegex= tokens.get(2);
		}
		try{
			Pattern.compile(trackNameRegex); // Validate regex
		} catch(PatternSyntaxException e){
	    	System.err.println("Invalid regex in: " + cmdInput);
	    	System.err.println(e.getDescription());
		}
		
		for(Track tr : this.trackSet.values()){
			boolean matched= Pattern.compile(trackNameRegex).matcher(tr.getFileTag()).find();
			if(matched){
				tr.setTitleColour(colour);
			}
		}
	}
	
	/** From cmdInput extract regex and ylimits then iterate through the tracks list to set 
	 * the ylimits in the tracks whose filename matches the regex.
	 * The input list is updated in place! 
	*/
	public void setTrackYlimitsForRegex(String cmdInput) throws InvalidCommandLineException{

		StrTokenizer str= new StrTokenizer(cmdInput);
		str.setTrimmerMatcher(StrMatcher.spaceMatcher());
		str.setQuoteChar('\'');
		List<String> tokens= str.getTokenList();
		if(tokens.size() < 3){
			System.err.println("Error in ylim subcommand. Expected at least 3 args got: " + cmdInput);
			throw new InvalidCommandLineException();
		}
		String trackNameRegex= ".*"; // Default: Capture everything
		if(tokens.size() == 4){
			trackNameRegex= tokens.get(3);
		}
		
		try{
			Pattern.compile(trackNameRegex); // Validate regex
		} catch(PatternSyntaxException e){
	    	System.err.println("Invalid regex in: " + cmdInput);
	    	System.err.println(e.getDescription());
		}
		
		double ymin= Double.NaN;
		try{
			ymin= Double.parseDouble(tokens.get(1));
		} catch(NumberFormatException e){
			ymin= Double.NaN;
		}
		double ymax= Double.NaN;
		try{
			ymax= Double.parseDouble(tokens.get(2));
		} catch(NumberFormatException e){
			ymax= Double.NaN;
		}

		if(ymin > ymax){ // Swap
			Double newMax= ymin;
			ymin= ymax;
			ymax= newMax;			
		}
		//if(ymin >= ymax){
		//	System.err.println("Warning ymin >= ymax. Resetting to default.");
		//	ymin= Double.NaN;
		//	ymax= Double.NaN;							
		//}
		for(Track tr : this.trackSet.values()){
			boolean matched= Pattern.compile(trackNameRegex).matcher(tr.getFileTag()).find();
			if(matched){
				tr.setYLimitMin(ymin);
				tr.setYLimitMax(ymax);
			}
		}
	}

	/** Set visibility for IntervalFeature tracks. 
	*/
	public void setVisibilityForTrackIntervalFeature(String cmdInput) throws InvalidCommandLineException{

		StrTokenizer str= new StrTokenizer(cmdInput);
		str.setTrimmerMatcher(StrMatcher.spaceMatcher());
		str.setQuoteChar('\'');
		List<String> tokens= str.getTokenList();

		// Defaults:
		String showRegex= ".*";  // Show all
		String hideRegex= "^$";    // Hide nothing
		String trackNameRegex= ".*"; // Apply to all tracks
		if(tokens.size() > 1){
			showRegex= tokens.get(1);
		}
		if(tokens.size() > 2){
			hideRegex= tokens.get(2);
		}
		if(tokens.size() > 3){
			trackNameRegex= tokens.get(3);
		}
		try{
			// Validate regex
			Pattern.compile(hideRegex); 
			Pattern.compile(showRegex);
			Pattern.compile(trackNameRegex); 
		} catch(PatternSyntaxException e){
	    	throw new PatternSyntaxException(e.getDescription(), cmdInput, -1);
		}

		System.err.println("Show: '" + showRegex + "'; hide: '" + hideRegex + "'; for tracks captured by '" + trackNameRegex + "':");
		for(Track tr : this.trackSet.values()){
			boolean matched= Pattern.compile(trackNameRegex).matcher(tr.getFileTag()).find();
			if(matched){
				System.err.println(tr.getFileTag());
				tr.setShowRegex(showRegex);
				tr.setHideRegex(hideRegex);
			}
		}
	}
	
	/** Go to the next feature on trackId given the current GenomicCoordinates. 
	 * 
	 * If slop is > 0, the output coordinates are centered on the feature and extended
	 * slop times the size of the feature left and right. With slop= 0 the coordinates 
	 * are exactly spanning the feature. With slop < 0 the output coordinates have 
	 * the feature right at the start.
	 * */
	public GenomicCoords goToNextFeatureOnFile(String trackId, GenomicCoords currentGc, double slop) throws InvalidGenomicCoordsException, IOException, InvalidCommandLineException{

		Track tr= matchIntervalFeatureTrack(trackId.trim());
		if(tr == null){
			return currentGc;
		}
		TrackIntervalFeature tif= (TrackIntervalFeature) tr;
		if(slop < 0){
			return tif.getIntervalFeatureSet().coordsOfNextFeature(currentGc);
		} else {
			GenomicCoords featureGc= tif.getIntervalFeatureSet().startEndOfNextFeature(currentGc);
			if(featureGc.equalCoords(currentGc)){ // No "next feature" found.
				return currentGc;
			} else {
				featureGc.centerAndExtendGenomicCoords(featureGc, featureGc.getGenomicWindowSize(), slop);
				return featureGc;
			}
		}
	}
	
	/** Convenient method to get interval feature tracks by name containing trackTag or matching trackTag by regex.
	 * If no matches are found or the trackSet is empty, return null. If multiple matches are found, 
	 * return the first one with warning.
	 * */
	private Track matchIntervalFeatureTrack(String trackTag){
		
		LinkedHashMap<String, Track> ifTracks = this.getIntervalFeatureTracks().getTrackSet();		
		Track tr= null;
		
		if(ifTracks.size() == 0){
			System.err.println("\nWarning interval feature track is empty.");
			return tr;
		}
		
		if(trackTag.isEmpty() && ifTracks.size() == 1){
			tr= ifTracks.values().iterator().next();
		} else if (trackTag.isEmpty() && ifTracks.size() > 1) {
			tr= ifTracks.values().iterator().next();
			System.err.println("\nWarning: trackId not given default to first track found: " + tr.getFileTag());
		} else {
			List<Track> matched= matchTracks(trackTag);
			if(matched.size() == 0){
				System.err.println("\nWarning '" + trackTag + "' not found in track set:");
				System.err.println(ifTracks.keySet() + "\n");
				return tr;
			} else {
				tr= matched.get(0);
				if(matched.size() > 1){
					System.err.println("\nWarning '" + trackTag + "' matches: " + matched + ". First track is returned.");
				}
			}
		}
		return tr;
	}


	/** Return the tracks whose trackId contains trackTag. If asRegex is true, matching is done by regex.
	 * */
	private List<Track> matchTracks(String trackTag){
		
		List<Track> matchedTracks= new ArrayList<Track>();
		
		Iterator<String> iter = this.trackSet.keySet().iterator();
		while(iter.hasNext()){
			String x= iter.next();
			boolean matched= Pattern.compile(trackTag).matcher(x).find();
			if(matched){
				matchedTracks.add(this.trackSet.get(x));
			}
		}
		return matchedTracks;
	}
	
	public GenomicCoords findNextMatchOnTrack(String query, String trackId, GenomicCoords currentGc, boolean all) throws InvalidGenomicCoordsException, IOException{

		TrackIntervalFeature tif= (TrackIntervalFeature) matchIntervalFeatureTrack(trackId.trim());
		if(tif == null){
			return currentGc;
		}

		System.err.println("Matching on " + tif.getFileTag());
		
		if(all){
			return tif.getIntervalFeatureSet().genomicCoordsAllChromMatchInGenome(query, currentGc);
		} else {
			return tif.getIntervalFeatureSet().findNextMatch(currentGc, query);
		}
	}

	private TrackSet getIntervalFeatureTracks(){
		TrackSet ifSet= new TrackSet();
		for(Track tr : this.trackSet.values()){
			if(Utils.getFileTypeFromName(tr.getFilename()).equals(TrackFormat.BED) 
			   || Utils.getFileTypeFromName(tr.getFilename()).equals(TrackFormat.GFF)
			   || Utils.getFileTypeFromName(tr.getFilename()).equals(TrackFormat.VCF)){
				ifSet.trackSet.put(tr.getFileTag(), tr);
			}
		}
		return ifSet;
	}
	
	public void selectDataColumnForBedgraph(int bdgDataColIdx, String trackIdRegex){
		
		for(Track tr : this.trackSet.values()){
			
			boolean matched= Pattern.compile(trackIdRegex).matcher(tr.getFileTag()).find();
			
			if(Utils.getFileTypeFromName(tr.getFilename()).equals(TrackFormat.BEDGRAPH) &&
					matched) {
		
				TrackWiggles bdg= (TrackWiggles) tr;
				bdg.setBdgDataColIdx(bdgDataColIdx);
			
			}
		}		
	}
	
	/** Reorder tracks with the one in newOrder. Tracks not in newOrder are appended
	 * with order unchanged.
	 * */
	public void orderTracks(List<String> newOrder) {

		// Create a new LinkedHashMap with the new order
		LinkedHashMap<String, Track> newTrackSet= new LinkedHashMap<String, Track>(); 
		for(String query : newOrder){
			List<Track> trList = this.matchTracks(query);
			for(Track xtrack : trList){
				if(!newTrackSet.containsKey(xtrack.getFileTag())){ // This will remove dups
					newTrackSet.put(xtrack.getFileTag(), xtrack);
				}
			}
		}
		
		// Append tracks not in newOrder
		for(String x : this.trackSet.keySet()){
			if(!newTrackSet.containsKey(x)){
				newTrackSet.put(x, this.trackSet.get(x));
			}
		}
		
		// A sanity check we didn't leave anything behind
		if(this.trackSet.size() != newTrackSet.size()){
			throw new RuntimeException("\nReordered track has " + newTrackSet.size() + " tracks. Expected " + this.trackSet.size());
		}
		for(String x : this.trackSet.keySet()){
			if(!newTrackSet.containsKey(x)){
				throw new RuntimeException("\nReordered track does not contain " + x);
			}
		}
		
		// Replace old with new hashmap
		this.trackSet= newTrackSet;
	}
	
	/*   S e t t e r s   and   G e t t e r s  */
	public LinkedHashMap<String, Track> getTrackSet() {
		return trackSet;
	}

}
