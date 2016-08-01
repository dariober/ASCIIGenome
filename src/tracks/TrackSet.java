package tracks;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.regex.Pattern;
import java.util.regex.PatternSyntaxException;

import exceptions.InvalidCommandLineException;
import exceptions.InvalidGenomicCoordsException;
import filter.FlagToFilter;
import htsjdk.samtools.filter.MappingQualityFilter;
import htsjdk.samtools.filter.SamRecordFilter;
import samTextViewer.GenomicCoords;
import samTextViewer.Utils;

/** Class to hold tracks to be printed. 
 * */
public class TrackSet {
	
	private LinkedHashMap<String, Track> trackSet= new LinkedHashMap<String, Track>();
	private List<Pattern> regexForTrackHeight= new ArrayList<Pattern>();
	private int trackHeightForRegex= -1;
	
	public static final String BOOKMARK_TAG= "bookmark";
	
	/*   C o n s t r u c t o r s   */
	
	public TrackSet(){}
	
	/*   M e t h o d s   */

	//public Track getTrackFromTag(){
	//	Track tr;
	//	for(){
	//		this.trackSet.get(key)
	//	}
	//	return tr;
	//}

	public String showTrackInfo(){
		List<String> trackInfo= new ArrayList<String>();
		Iterator<Entry<String, Track>> trx= this.trackSet.entrySet().iterator();
		
		while(trx.hasNext()){
			Entry<String, Track> x = trx.next();
			String hd= x.getValue().getyMaxLines() <= 0 ? "*" : "";
			trackInfo.add(x.getKey() + "\t" 
					+ x.getValue().getFilename() + "\t" 
					+ Utils.getFileTypeFromName(x.getValue().getFilename()) + "\t"
					+ hd);
		}
		StringBuilder sb= new StringBuilder();
		for(String str : Utils.tabulateList(trackInfo)){
			sb.append(str + "\n");
		}
		return sb.toString().trim();
	}


	/** From cmdInput extract regex and yMaxLines then iterate through the tracks list to set 
	 * the yMaxLines in the tracks whose filename matches the regex.
	 * The input list is updated in place! 
	*/
	public void setTrackHeightForRegex(List<String> tokens) throws InvalidCommandLineException{

		// MEMO of subcommand syntax:
		// 0 trackHeight
		// 1 int    mandatory
		// 2 regex  optional
		
		if(tokens.size() < 2){
			System.err.println("Error in trackHeight subcommand. Expected 2 args got: " + tokens);
			throw new InvalidCommandLineException();
		}

		// Get height
		try{
			this.trackHeightForRegex= Integer.parseInt(tokens.get(1));
			this.trackHeightForRegex= this.trackHeightForRegex < 0 ? 0 : this.trackHeightForRegex;
		} catch(NumberFormatException e){
			System.err.println("Number format exception: " + this.trackHeightForRegex);
			throw new InvalidCommandLineException();
		}

		
        // Regex
        List<String> trackNameRegex= new ArrayList<String>();
        if(tokens.size() >= 3){
            trackNameRegex= tokens.subList(2, tokens.size());
        } else {
            trackNameRegex.add(".*"); // Default: Capture everything
        }
        
        // Update
        this.regexForTrackHeight.clear();
        for(String x : trackNameRegex){
        	try{
        		this.regexForTrackHeight.add(Pattern.compile(x));
        	} catch(PatternSyntaxException e){
        		System.err.println("Command: " + tokens);
        		System.err.println("Invalid regex in: " + x);
		    	System.err.println(e.getDescription());
        		throw new InvalidCommandLineException();
        	}
        }
        
        // And set as required:
        List<Track> tracksToReset = this.matchTracks(trackNameRegex, true);
        for(Track tr : tracksToReset){
        	tr.setyMaxLines(this.trackHeightForRegex);
        }
	}

	public void setFeatureDisplayModeForRegex(List<String> tokens) throws InvalidCommandLineException {
		// MEMO of subcommand syntax:
		// 0 squash/merge
		// 1 Regex
		
        // Regex
        List<String> trackNameRegex= new ArrayList<String>();
        if(tokens.size() >= 2){
            trackNameRegex= tokens.subList(1, tokens.size());
        } else {
            trackNameRegex.add(".*"); // Default: Capture everything
        }
                
		FeatureDisplayMode switchTo;
		if(tokens.get(0).equals("squash")){
			switchTo = FeatureDisplayMode.SQUASHED;
		} else if(tokens.get(0).equals("merge")){
			switchTo = FeatureDisplayMode.MERGED;
		} else {
			System.err.println("Unexepected command: " + tokens);
			throw new InvalidCommandLineException();
		}
		
        // And set as required:
		List<Track> tracksToReset = this.matchTracks(trackNameRegex, true);
        for(Track tr : tracksToReset){
        	if(tr.getFeatureDisplayMode().equals(switchTo)){ // Invert setting
				tr.setFeatureDisplayMode(FeatureDisplayMode.EXPANDED);
			} else {
				tr.setFeatureDisplayMode(switchTo);
			}
        }
	}
	
	public void setPrintModeForRegex(List<String> tokens) throws InvalidCommandLineException {
		// MEMO of subcommand syntax:
		// 0 print/printFull
		// 1 Regex
		
        // Regex
        List<String> trackNameRegex= new ArrayList<String>();
        if(tokens.size() >= 2){
            trackNameRegex= tokens.subList(1, tokens.size());
        } else {
            trackNameRegex.add(".*"); // Default: Capture everything
        }
                
		PrintRawLine switchTo; // If printing is OFF do we switch to CLIP or FULL?
		if(tokens.get(0).equals("print")){
			switchTo = PrintRawLine.CLIP;
		} else if(tokens.get(0).equals("printFull")){
			switchTo = PrintRawLine.FULL;
		} else {
			System.err.println("Unexepected command: " + tokens);
			throw new InvalidCommandLineException();
		}

        // And set as required:
		List<Track> tracksToReset = this.matchTracks(trackNameRegex, true);
        for(Track tr : tracksToReset){
        	if(tr.getPrintMode().equals(switchTo)){ // Invert setting
        		tr.setPrintMode(PrintRawLine.OFF);
			} else {
				tr.setPrintMode(switchTo);
			}
        }
	}
	
	public void setBisulfiteModeForRegex(List<String> tokens) throws InvalidCommandLineException {

		// MEMO of subcommand syntax:
		// 0 BSseq
		// 1 Regex
		
        // Regex
        List<String> trackNameRegex= new ArrayList<String>();
        if(tokens.size() >= 2){
            trackNameRegex= tokens.subList(1, tokens.size());
        } else {
            trackNameRegex.add(".*"); // Default: Capture everything
        }
        // And set as required:
        List<Track> tracksToReset = this.matchTracks(trackNameRegex, true);
        for(Track tr : tracksToReset){
			if(tr.isBisulf()){ // Invert setting
				tr.setBisulf(false);
			} else {
				tr.setBisulf(true);
			}
        }
	}

	/*
	public void setPrintPileupForRegex(List<String> tokens) throws InvalidCommandLineException {

		// MEMO of subcommand syntax:
		// 0 pileup
		// 1 Regex
		
        // Regex
        List<String> trackNameRegex= new ArrayList<String>();
        if(tokens.size() >= 2){
            trackNameRegex= tokens.subList(1, tokens.size());
        } else {
            trackNameRegex.add(".*"); // Default: Capture everything
        }
        // And set as required:
        List<Track> tracksToReset = this.matchTracks(trackNameRegex, true);
        for(Track tr : tracksToReset){
			if(tr.isPrintPileup()){ // Invert setting
				tr.setPrintPileup(false);
			} else {
				tr.setPrintPileup(true);
			}
        }
	} */
	
	
	public void setRpmForRegex(List<String> tokens) throws InvalidCommandLineException {

		// MEMO of subcommand syntax:
		// 0 rpm
		// 1 Regex
		
        // Regex
        List<String> trackNameRegex= new ArrayList<String>();
        if(tokens.size() >= 2){
            trackNameRegex= tokens.subList(1, tokens.size());
        } else {
            trackNameRegex.add(".*"); // Default: Capture everything
        }
        // And set as required:
        List<Track> tracksToReset = this.matchTracks(trackNameRegex, true);
        for(Track tr : tracksToReset){
			if(tr.isRpm()){ // Invert setting
				tr.setRpm(false);
			} else {
				tr.setRpm(true);
			}
        }
	}
	public void setTrackColourForRegex(List<String> tokens) throws InvalidCommandLineException{

		// MEMO of subcommand syntax:
		// 0 trackColour
		// 1 Colour
		// 2 Regex

		// Colour
		String colour= (new Track()).getTitleColour();
		if(tokens.size() >= 2){
			String xcolour= tokens.get(1).toLowerCase();
			if(!Utils.ansiColorCodes().containsKey(xcolour)){
				System.err.println("\nInvalid colour: " + xcolour);
				System.err.println("Valid colours are: " + Utils.ansiColorCodes().keySet());
				throw new InvalidCommandLineException();
			} else {
				colour= xcolour;
			}
		}
		
		// Regex
		List<String> trackNameRegex= new ArrayList<String>();
		if(tokens.size() >= 3){
			trackNameRegex= tokens.subList(2, tokens.size());
		} else {
			trackNameRegex.add(".*"); // Default: Capture everything
		}
		// And set as required:
		List<Track> tracksToReset = this.matchTracks(trackNameRegex, true);
		for(Track tr : tracksToReset){
			tr.setTitleColour(colour);
		}
	}
	
	public void setAttributeForGFFName(List<String> tokens) throws InvalidCommandLineException{

		// MEMO of subcommand syntax:
		// 0 gffNameAttr
		// 1 attrName
		// 2 Regex

		String gtfAttributeForName= null; // Null will follow default 
		if(tokens.size() >= 2){
			gtfAttributeForName= tokens.get(1);
			if(gtfAttributeForName.equals("NULL")){
				gtfAttributeForName= null;
			}
		}

        // Regex
        List<String> trackNameRegex= new ArrayList<String>();
        if(tokens.size() >= 3){
            trackNameRegex= tokens.subList(2, tokens.size());
        } else {
            trackNameRegex.add(".*"); // Default: Capture everything
        }
        // And set as required:
        List<Track> tracksToReset = this.matchTracks(trackNameRegex, true);
        for(Track tr : tracksToReset){
        	tr.setGtfAttributeForName(gtfAttributeForName);
        }		
	}
	
	/** From cmdInput extract regex and ylimits then iterate through the tracks list to set 
	 * the ylimits in the tracks whose filename matches the regex.
	 * The input list is updated in place! 
	*/
	public void setTrackYlimitsForRegex(List<String> tokens) throws InvalidCommandLineException{

		// MEMO of subcommand syntax:
		// 0 cmdName
		// 1 min
		// 2 max
		// 3+ regex (opt)

		if(tokens.size() < 3){
			System.err.println("Error in ylim subcommand. Expected at least 2 args got: " + tokens);
			throw new InvalidCommandLineException();
		}
		// Parse min and max
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

        // Regex
        List<String> trackNameRegex= new ArrayList<String>();
        if(tokens.size() >= 4){
            trackNameRegex= tokens.subList(2, tokens.size());
        } else {
            trackNameRegex.add(".*"); // Default: Capture everything
        }
        // And set as required:
        List<Track> tracksToReset = this.matchTracks(trackNameRegex, true);
        for(Track tr : tracksToReset){
    		tr.setYLimitMin(ymin);
			tr.setYLimitMax(ymax);
        }
	}

	/** Set filter for IntervalFeature tracks. 
	*/
	public void setFilterForTrackIntervalFeature(List<String> tokens) throws InvalidCommandLineException{

		// 0 cmdName
		// 1 showRe
		// 2 hideRe
		// 3 trackRe+
		
		// SHOW REGEX
		String showRegex= ".*";  // Show all
		if(tokens.size() > 1){
			showRegex= tokens.get(1);
		}
		try{
			Pattern.compile(showRegex);
		} catch(PatternSyntaxException e){
	    	System.err.println("Invalid regex in: " + tokens);
	    	System.err.println("showRegex: " + showRegex);
	    	throw new InvalidCommandLineException();
		}
		
		// HIDE REGEX
		String hideRegex= "^$";    // Hide nothing
		if(tokens.size() > 2){
			hideRegex= tokens.get(2);
		}
		try{
			Pattern.compile(hideRegex); 
		} catch(PatternSyntaxException e){
	    	System.err.println("Invalid regex in: " + tokens);
	    	System.err.println("hideRegex: " + hideRegex);
	    	throw new InvalidCommandLineException();
	    }
		
		// TRACK REGEXES
        // Regex
        List<String> trackNameRegex= new ArrayList<String>();
        if(tokens.size() >= 4){
            trackNameRegex= tokens.subList(2, tokens.size());
        } else {
            trackNameRegex.add(".*"); // Default: Capture everything
        }
        // And set as required:
        List<Track> tracksToReset = this.matchTracks(trackNameRegex, true);
        for(Track tr : tracksToReset){
			tr.setShowRegex(showRegex);
			tr.setHideRegex(hideRegex);
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
	 * @throws InvalidCommandLineException 
	 * */
	private Track matchIntervalFeatureTrack(String trackTag) throws InvalidCommandLineException{
		
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
			List<String> x= new ArrayList<String>();
			x.add(trackTag);
			List<Track> matched= matchTracks(x, true);
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

	private List<Track> matchTracks(List<String> patterns, boolean asRegex) throws InvalidCommandLineException{

		// Validate regexes
		if(asRegex){
			for(String x : patterns){
				try{
					Pattern.compile(x); 
				} catch(PatternSyntaxException e){
			    	System.err.println("Invalid regex: " + x);
			    	throw new InvalidCommandLineException();
				}		
			}
		}
		
		List<Track> matchedTracks= new ArrayList<Track>();
		
		Iterator<String> iter = this.trackSet.keySet().iterator();
		while(iter.hasNext()){
			String trackId= iter.next();
			for(String pattern : patterns){
				boolean matched= Pattern.compile(pattern).matcher(trackId).find();
				if(matched && !matchedTracks.contains(trackId)){
					matchedTracks.add(this.trackSet.get(trackId));
				}
			}
		}
		return matchedTracks;		
	}
	
	public GenomicCoords findNextMatchOnTrack(String query, String trackId, GenomicCoords currentGc, boolean all) throws InvalidGenomicCoordsException, IOException, InvalidCommandLineException{

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
	
	public void setDataColForRegex(List<String> tokens) throws InvalidCommandLineException{

		// MEMO of subcommand syntax:
		// 0 gffNameAttr
		// 1 attrName
		// 2 Regex

		int dataCol; // Null will follow default 
		if(tokens.size() >= 2){
			try{
				dataCol= Integer.parseInt(tokens.get(1));
			} catch(NumberFormatException e){
				System.err.println("Number format exception: " + tokens.get(1));
				throw new InvalidCommandLineException();
			}
		} else {
			dataCol= 4;
		}

        // Regex
        List<String> trackNameRegex= new ArrayList<String>();
        if(tokens.size() >= 3){
            trackNameRegex= tokens.subList(2, tokens.size());
        } else {
            trackNameRegex.add(".*"); // Default: Capture everything
        }
        // And set as required:
        List<Track> tracksToReset = this.matchTracks(trackNameRegex, true);
        for(Track tr : tracksToReset){
        	if(Utils.getFileTypeFromName(tr.getFilename()).equals(TrackFormat.BEDGRAPH)){
				TrackWiggles bdg= (TrackWiggles) tr;
				bdg.setBdgDataColIdx(dataCol);	
			}
        }		
	}
	
	/** Reorder tracks with the one in newOrder. Tracks not in newOrder are appended
	 * with order unchanged.
	 * @throws InvalidCommandLineException 
	 * */
	public void orderTracks(List<String> newOrder) throws InvalidCommandLineException {

		// Create a new LinkedHashMap with the new order
		LinkedHashMap<String, Track> newTrackSet= new LinkedHashMap<String, Track>(); 
		for(String query : newOrder){
			List<String> x= new ArrayList<String>();
			x.add(query);
			List<Track> trList = this.matchTracks(x, true);
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

	public List<Pattern> getRegexForTrackHeight() {
		return regexForTrackHeight;
	}
	public int getTrackHeightForRegex() {
		return trackHeightForRegex;
	}

	public void addBookmark_IN_PREP(GenomicCoords gc, String name) throws IOException, InvalidGenomicCoordsException {
		
		String raw= (gc.getChrom() + "\t" + (gc.getFrom()-1) + "\t" + gc.getTo() + "\t" + name).trim();
		
		IntervalFeature interval= new IntervalFeature(raw, TrackFormat.BED);
		// This is the map of intervalFeatures that will be used to construct the bookamrk track.
		TrackIntervalFeature bookmarkTrack;
		if(!this.trackSet.containsKey(TrackSet.BOOKMARK_TAG)){
		
			List<IntervalFeature> xl= new ArrayList<IntervalFeature>();
			xl.add(interval);
			Map<String, List<IntervalFeature>> map= new HashMap<String, List<IntervalFeature>>();
			map.put(gc.getChrom(), xl); 
			bookmarkTrack = new TrackIntervalFeature(new IntervalFeatureSet(map, TrackFormat.BED), gc);
			bookmarkTrack.setFileTag(TrackSet.BOOKMARK_TAG);
		
		} else {
			
			// Get the map of the IntervalFeature bookmarks
			Map<String, List<IntervalFeature>> map= new HashMap<String, List<IntervalFeature>>();
			map = ((TrackIntervalFeature)this.trackSet.get(TrackSet.BOOKMARK_TAG)).
					getIntervalFeatureSet().
					getIntervalMap();

			// Get the list of bookmarks on this chrom
			List<IntervalFeature> xlist= new ArrayList<IntervalFeature>();
			if(map.containsKey(interval.getChrom())){
				xlist = map.get(interval.getChrom());
			}
			// Add the new bookmark
			xlist.add(interval);

			// Update map with augmented list 
			map.put(interval.getChrom(), xlist);

			// Recreate the IntervalFeatureSet. You HAVE TO recreate it in order to have the list sorted.
			IntervalFeatureSet intervalFeatureSet= new IntervalFeatureSet(map, TrackFormat.BED);

			// In the existing bookmark track, replace the old intervalset with the newly created one
			bookmarkTrack = (TrackIntervalFeature) this.trackSet.get(TrackSet.BOOKMARK_TAG);
			bookmarkTrack.setIntervalFeatureSet(intervalFeatureSet);
		}
		this.trackSet.put(TrackSet.BOOKMARK_TAG, bookmarkTrack);
	}

	/** Method to set any of the three alignment filters: -F,-f, mapq */
	public void setFilterFlagForRegex(List<String> tokens) throws InvalidCommandLineException {
		// MEMO of subcommand syntax:
		// 0 -F
		// 1 INT
		// 2... Regexes

		if(tokens.size() < 2){
			System.err.println("Expected at least two arguments. Got: " + tokens);
			throw new InvalidCommandLineException();
		}
		
		int flag= 0; // Null will follow default 
		if(tokens.size() >= 2){
			try{
				flag= Integer.parseInt(tokens.get(1));
				if(tokens.get(0).equals("-F") && (4 & flag) == 0){
					flag += 4;
				}
			} catch(NumberFormatException e){
				System.err.println("Number format exception: " + tokens.get(1));
				throw new InvalidCommandLineException();
			}
		} 

        // Regex
        List<String> trackNameRegex= new ArrayList<String>();
        if(tokens.size() >= 3){
            trackNameRegex= tokens.subList(2, tokens.size());
        } else {
            trackNameRegex.add(".*"); // Default: Capture everything
        }
        // And set as required:
        List<Track> tracksToReset = this.matchTracks(trackNameRegex, true);
        for(Track tr : tracksToReset){
            
        	List<SamRecordFilter> filters= null;
        	
        	// NMB: When you set one filetr, e.g. mapq, you have to restore the others!
            if(tokens.get(0).equals("-F")){
        		tr.set_F_flag(flag);
        		filters= FlagToFilter.flagToFilterList(tr.get_f_flag(), flag);
        		filters.add(new MappingQualityFilter(tr.getMapq()));
        		
        	} else if(tokens.get(0).equals("-f")){
        		tr.set_f_flag(flag);
        		filters= FlagToFilter.flagToFilterList(flag, tr.get_F_flag());
        		filters.add(new MappingQualityFilter(tr.getMapq()));
        		
        	} else if(tokens.get(0).equals("mapq")){
        		tr.setMapq(flag);
        		filters= FlagToFilter.flagToFilterList(tr.get_f_flag(), tr.get_F_flag());
        		filters.add(new MappingQualityFilter(flag));
        	
        	} else {
				System.err.println("Unexpected command: " + tokens);
				throw new InvalidCommandLineException();
        	}
        	tr.setSamRecordFilter(filters);
        }		
	}

	public void setFeatureGapForRegex(List<String> tokens) throws InvalidCommandLineException {

		// MEMO of subcommand syntax:
		// 0 gap
		// 1 Regex
		
        List<String> trackNameRegex= new ArrayList<String>();
        if(tokens.size() >= 2){
            trackNameRegex= tokens.subList(1, tokens.size());
        } else {
            trackNameRegex.add(".*"); // Default: Capture everything
        }

        // And set as required:
        List<Track> tracksToReset = this.matchTracks(trackNameRegex, true);
        for(Track tr : tracksToReset){
			if(tr.getGap() == 0){ // Invert setting
				tr.setGap(1);
			} else {
				tr.setGap(0);
			}
        }
	}

	/*
	public void set_f_flagForRegex(ArrayList<String> tokens) throws InvalidCommandLineException {
		// MEMO of subcommand syntax:
		// 0 -F
		// 1 INT
		// 2... Regexes

		if(tokens.size() < 2){
			System.err.println("Expected at least two arguments. Got: " + tokens);
			throw new InvalidCommandLineException();
		}
		
		int flag= 0; // Null will follow default 
		try{
			flag= Integer.parseInt(tokens.get(1));
		} catch(NumberFormatException e){
			System.err.println("Number format exception: " + tokens.get(1));
			throw new InvalidCommandLineException();
		}
		
        // Regex
        List<String> trackNameRegex= new ArrayList<String>();
        if(tokens.size() >= 3){
            trackNameRegex= tokens.subList(2, tokens.size());
        } else {
            trackNameRegex.add(".*"); // Default: Capture everything
        }
        // And set as required:
        List<Track> tracksToReset = this.matchTracks(trackNameRegex, true);
        for(Track tr : tracksToReset){
        	tr.set_f_flag(flag);
        }			
	}*/
	
}
