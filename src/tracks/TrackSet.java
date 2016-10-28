package tracks;

import java.io.File;
import java.io.IOException;
import java.net.MalformedURLException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.regex.Pattern;
import java.util.regex.PatternSyntaxException;

import exceptions.BamIndexNotFoundException;
import exceptions.InvalidCommandLineException;
import exceptions.InvalidGenomicCoordsException;
import exceptions.InvalidRecordException;
import filter.FlagToFilter;
import htsjdk.samtools.filter.MappingQualityFilter;
import htsjdk.samtools.filter.SamRecordFilter;
import samTextViewer.GenomicCoords;
import samTextViewer.Utils;

/** Class to hold tracks to be printed. 
 * */
public class TrackSet {
	
	// private LinkedHashMap<String, Track> trackSet_DEPRECATED= new LinkedHashMap<String, Track>();
	private List<Track> trackList= new ArrayList<Track>();
	private List<Pattern> regexForTrackHeight= new ArrayList<Pattern>();
	private int trackHeightForRegex= -1;
	public static final String BOOKMARK_TAG= "bookmark";
	
	/*   C o n s t r u c t o r s   */
	
	public TrackSet(){}
	
	public TrackSet(List<String> inputFileList, GenomicCoords gc) throws IOException, InvalidGenomicCoordsException, InvalidRecordException, ClassNotFoundException, SQLException{
		
		for(String sourceName : inputFileList){

			if(Utils.getFileTypeFromName(sourceName).equals(TrackFormat.BAM)){
				//
				// BAM FILE
				//
				if(!Utils.bamHasIndex(sourceName)){
					System.err.println("\nNo index found for '" + sourceName + "'. Index can be generated with ");
					System.err.println("samtools index '" + sourceName + "'\n");
					System.exit(1);
				}
				
				/* Coverage track */
				TrackCoverage trackCoverage= new TrackCoverage(sourceName, gc, false);

				trackCoverage.setId(this.getMaxTrackId() + 1);
				trackCoverage.setTrackTag(new File(sourceName).getName() + "#" + (this.getMaxTrackId() + 1));
				this.trackList.add(trackCoverage);
				
				/* Read track */
				TrackReads trackReads= new TrackReads(sourceName, gc);
				trackReads.setId(this.getMaxTrackId() + 1);
				trackReads.setTrackTag(new File(sourceName).getName() + "@" + (this.getMaxTrackId() + 1));
				this.trackList.add(trackReads);
			}
			
			else if(    Utils.getFileTypeFromName(sourceName).equals(TrackFormat.BED) 
		        || Utils.getFileTypeFromName(sourceName).equals(TrackFormat.GFF)
			    || Utils.getFileTypeFromName(sourceName).equals(TrackFormat.VCF)){
				//
				// Annotatation
				//
				TrackIntervalFeature tif= new TrackIntervalFeature(sourceName, gc);
				//tif.setTrackTag(new File(sourceName).getName() + "#" + (this.getMaxTrackId()+1));
				//tif.setId(this.getMaxTrackId() + 1);
				this.add(tif, new File(sourceName).getName());
			} 

			else if(Utils.getFileTypeFromName(sourceName).equals(TrackFormat.BIGWIG) 
					|| Utils.getFileTypeFromName(sourceName).equals(TrackFormat.TDF) 
					|| Utils.getFileTypeFromName(sourceName).equals(TrackFormat.BEDGRAPH)){
				//
				// Wiggles
				//
				TrackWiggles tw= new TrackWiggles(sourceName, gc, 4);
				//tw.setTrackTag(new File(sourceName).getName() + "#" + (this.getMaxTrackId()+1));
				//tw.setId(this.getMaxTrackId()+1);
				this.add(tw, new File(sourceName).getName());
				
			} else {
				System.err.println("Unable to classify " + sourceName + "; skipping"); 								
			}			
		}		
		// TrackWiggles gcProfile= gc.getGCProfile();
	}
	
	/*   M e t h o d s   */

	/**
	 * Add this track with given baseTag. The suffix "#id" will be appended to the baseTag string. 
	 * NB: Adding a track resets the track's ID
	 * */
	public void add(Track track, String baseTag) {
		int idForTrack= this.getMaxTrackId() + 1;
		String trackTag= baseTag + "#" + idForTrack;
		track.setTrackTag(trackTag);
		track.setId(idForTrack);
		this.trackList.add(track);
	}

	/** Add track from given file or URL.
	 * @throws BamIndexNotFoundException 
	 * @throws IOException 
	 * @throws InvalidGenomicCoordsException 
	 * @throws InvalidRecordException 
	 * @throws SQLException 
	 * @throws ClassNotFoundException 
	 * */
	public void add(String sourceName, GenomicCoords gc) throws IOException, BamIndexNotFoundException, InvalidGenomicCoordsException, InvalidRecordException, ClassNotFoundException, SQLException{

		if(Utils.getFileTypeFromName(sourceName).equals(TrackFormat.BAM)){
			this.addBamTrackFromSourceName(sourceName, gc);
		
		} else if(Utils.getFileTypeFromName(sourceName).equals(TrackFormat.BED) 
		          || Utils.getFileTypeFromName(sourceName).equals(TrackFormat.GFF)
			      || Utils.getFileTypeFromName(sourceName).equals(TrackFormat.VCF)){
			this.addIntervalFeatureTrackFromSourceName(sourceName, gc);
		
		} else if(Utils.getFileTypeFromName(sourceName).equals(TrackFormat.BIGWIG) 
				|| Utils.getFileTypeFromName(sourceName).equals(TrackFormat.TDF) 
				|| Utils.getFileTypeFromName(sourceName).equals(TrackFormat.BEDGRAPH)){
			this.addWiggleTrackFromSourceName(sourceName, gc);
		}
	}
	
	private void addWiggleTrackFromSourceName(String sourceName, GenomicCoords gc) throws IOException, InvalidRecordException, InvalidGenomicCoordsException{
		
		int idForTrack= this.getMaxTrackId() + 1;
		String trackId= new File(sourceName).getName() + "#" + idForTrack;
		
		TrackWiggles tw= new TrackWiggles(sourceName, gc, 4);
		tw.setId(idForTrack);
		tw.setTrackTag(trackId);
		this.trackList.add(tw);
	}
	
	private void addIntervalFeatureTrackFromSourceName(String sourceName, GenomicCoords gc) throws IOException, InvalidGenomicCoordsException, ClassNotFoundException, InvalidRecordException, SQLException{
		
		int idForTrack= this.getMaxTrackId() + 1;
		
		String trackId= new File(sourceName).getName() + "#" + (idForTrack);
		TrackIntervalFeature tif= new TrackIntervalFeature(sourceName, gc);
		tif.setTrackTag(trackId);
		tif.setId(idForTrack);
		this.trackList.add(tif);
	}
	
	private void addBamTrackFromSourceName(String sourceName, GenomicCoords gc) throws IOException, BamIndexNotFoundException, InvalidGenomicCoordsException{

		int idForTrack= this.getMaxTrackId() + 1;
		
		if(!Utils.bamHasIndex(sourceName)){
			System.err.println("\nNo index found for '" + sourceName + "'. Index can be generated with ");
			System.err.println("samtools index '" + sourceName + "'\n");
			throw new BamIndexNotFoundException();
		}
		
		/* BAM Coverage track */
		String coverageTrackId= new File(sourceName).getName() + "#" + (idForTrack);
		
		TrackCoverage trackCoverage= new TrackCoverage(sourceName, gc, false);
		trackCoverage.setId(idForTrack);
		trackCoverage.setTrackTag(coverageTrackId);
		this.trackList.add(trackCoverage);
		
		/* Reads */
		idForTrack= this.getMaxTrackId() + 1;
		String trackId= new File(sourceName).getName() + "@" + (idForTrack+1);

		TrackReads trackReads= new TrackReads(sourceName, gc);
		trackReads.setTrackTag(trackId);
		trackReads.setId(idForTrack);
		trackReads.setTrackTag(trackId);

		this.trackList.add(trackReads);
	} 
	
	private int getMaxTrackId(){
		
		if(this.trackList.size() == 0){
			return 0;
		}
		
		int id= Integer.MIN_VALUE;
		for(Track tr : this.trackList){
			if(tr.getId() > id){
				id= tr.getId();
			}
		}
		return id; 
	}
	
	/** Show track info as nicely printable string.  
	 * */
	public String showTrackInfo(){
		List<String> trackInfo= new ArrayList<String>();
		
		for(Track track : this.getTrackList()){
			String hd= track.getyMaxLines() <= 0 ? "*" : "";
			trackInfo.add(track.getTrackTag() + "\t" 
					+ track.getFilename() + "\t" 
					+ Utils.getFileTypeFromName(track.getFilename()) + "\t"
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
        
        // Update
//        this.regexForPrintMode.clear();
//        for(String x : trackNameRegex){
//        	try{
//        		this.regexForPrintMode.add(Pattern.compile(x));
//        	} catch(PatternSyntaxException e){
//        		System.err.println("Command: " + tokens);
//        		System.err.println("Invalid regex in: " + x);
//		    	System.err.println(e.getDescription());
//        		throw new InvalidCommandLineException();
//        	}
//        }
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

	public void setHideTitleForRegex(List<String> tokens) throws InvalidCommandLineException {

		// MEMO of subcommand syntax:
		// 0 hideTitle
		// 1 Regex
		
		if(tokens.size() == 2 && tokens.get(1).equals("/hide_all/")){
			for(Track tr : this.getTrackList()){
				tr.setHideTitle(true);
			}
			return;
		}
		if(tokens.size() == 2 && tokens.get(1).equals("/show_all/")){
			for(Track tr : this.getTrackList()){
				tr.setHideTitle(false);
			}
			return;
		}		
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
			if(tr.isHideTitle()){ // Invert setting
				tr.setHideTitle(false);
			} else {
				tr.setHideTitle(true);
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
			return tif.coordsOfNextFeature(currentGc);
		} else {
			GenomicCoords featureGc= tif.startEndOfNextFeature(currentGc);
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

		List<Track> ifTracks = this.getIntervalFeatureTracks().getTrackList();	
				
		Track tr= null;
		
		if(ifTracks.size() == 0){
			System.err.println("\nWarning interval feature track is empty.");
			return tr;
		}
		
		if(trackTag.isEmpty() && ifTracks.size() == 1){
			tr= ifTracks.get(0);
		} else if (trackTag.isEmpty() && ifTracks.size() > 1) {
			tr= ifTracks.get(0);
			System.err.println("\nWarning: trackId not given default to first track found: " + tr.getTrackTag());
		} else {
			List<String> x= new ArrayList<String>();
			x.add(trackTag);
			List<Track> matched= matchTracks(x, true);
			if(matched.size() == 0){
				System.err.println("\nWarning '" + trackTag + "' not found in track set:");
				System.err.println(ifTracks + "\n");
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
		
		for(Track track : this.getTrackList()){
			String trackId= track.getTrackTag();
			for(String pattern : patterns){
				boolean matched= Pattern.compile(pattern).matcher(trackId).find();
				if(matched && !matchedTracks.contains(trackId)){
					matchedTracks.add(track);
				}
			}
		}
		return matchedTracks;		
	}
	
	/** Simple method to get track from track object. See also this.getTrackFromTag. Return null if track not found. 
	 * */
	protected Track getTrack(Track track){
		
		int idx= this.getTrackList().indexOf(track);
		if(idx == -1){
			return null;
		}
		return this.getTrackList().get(idx);
		
	}

	/** Get track given a track tag. See also this.getTrack. Returns null if trackTag not found.
	 * */
	public Track getTrackFromTag(String trackTag){
		
		for(Track track : this.getTrackList()){
			if(track.getTrackTag().equals(trackTag)){
				return track;
			}
		}
		return null;
	}

	
	public GenomicCoords findNextMatchOnTrack(String query, String trackId, GenomicCoords currentGc, boolean all) throws InvalidGenomicCoordsException, IOException, InvalidCommandLineException{

		TrackIntervalFeature tif= (TrackIntervalFeature) matchIntervalFeatureTrack(trackId.trim());
		if(tif == null){
			return currentGc;
		}

		System.err.println("Matching on " + tif.getTrackTag());
		
		if(all){
			return tif.genomicCoordsAllChromMatchInGenome(query, currentGc);
		} else {
			return tif.findNextMatch(currentGc, query);
		}
	}

	private TrackSet getIntervalFeatureTracks(){
		TrackSet ifSet= new TrackSet();
		for(Track tr : this.getTrackList()){
			if(Utils.getFileTypeFromName(tr.getFilename()).equals(TrackFormat.BED) 
			   || Utils.getFileTypeFromName(tr.getFilename()).equals(TrackFormat.GFF)
			   || Utils.getFileTypeFromName(tr.getFilename()).equals(TrackFormat.VCF)){
				ifSet.add(tr, new File(tr.getFilename()).getName());
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
	
		if(newOrder.size() == 0){
			this.sortTracksByTagName();
			return;
		}
		
		// Create a new list that will have the new order
		List<Track> newTrackList= new ArrayList<Track>();
		
		for(String query : newOrder){
			List<String> x= new ArrayList<String>();
			x.add(query);
			List<Track> trList = this.matchTracks(x, true);
			for(Track xtrack : trList){
				if(!newTrackList.contains(xtrack)){ // This will remove dups
					newTrackList.add(xtrack);
				}
			}
		}
		
		// Append tracks not in newOrder
		for(Track xtrack : this.getTrackList()){
			if(!newTrackList.contains(xtrack)){
				newTrackList.add(xtrack);
			}			
		}
		
		// A sanity check we didn't leave anything behind
		if(this.getTrackList().size() != newTrackList.size()){
			throw new RuntimeException("\nReordered track has " + newTrackList.size() + " tracks. Expected " + this.getTrackList().size());
		}
		for(Track x : this.getTrackList()){
			if(!newTrackList.contains(x)){
				throw new RuntimeException("\nReordered track does not contain " + x.getTrackTag());
			}
		}
		
		// Replace old with new hashmap
		this.trackList= newTrackList;
	}
	
	private void sortTracksByTagName(){
		
		Collections.sort(this.trackList, new Comparator<Track>() {
		    @Override
		    public int compare(Track o1, Track o2) {
		        return o1.getTrackTag().compareTo(o2.getTrackTag());
		    }
		});
	}
	
	/** Simple method to get list of track tags
	 * */
	protected List<String> getTrackTags(){
		
		List<String> trackTags= new ArrayList<String>();
		for(Track tr : this.getTrackList()){
			trackTags.add(tr.getTrackTag());
		}
		return trackTags;
		
	}

	/** Simple method to get list of filenames
	 * */
	public List<String> getFilenameList(){
		
		List<String> filenames= new ArrayList<String>();
		for(Track tr : this.getTrackList()){
			filenames.add(tr.getTrackTag());
		}
		return filenames;
		
	}
	
	/** Drop from TrackSet the track with given hashcode. 
	 * Return true if the track was found and dropped, false otherwise. 
	 * If trackHashCode do nothing and return false.  
	 * */
	public boolean dropTrackWithHashCode(Integer trackHashCode){
		
		if(trackHashCode == null){
			return false;
		}
		
		for(int i= 0; i < this.trackList.size(); i++){
			if(this.trackList.get(i).hashCode() == trackHashCode){
				this.trackList.remove(i);
				return true;
			}
		}
		return false;
	}
	
	/** Drop from TrackSet the track with given tag. 
	 * Return true if the track was found and dropped, false otherwise. 
	 * If trackTag is null do nothing and return false.  
	 * */
	public boolean dropTrackWithTrackTag(String trackTag){
		
		if(trackTag == null){
			return false;
		}
		
		for(int i= 0; i < this.trackList.size(); i++){
			if(this.trackList.get(i).getTrackTag().equals(trackTag)){
				this.trackList.remove(i);
				return true;
			}
		}
		return false;
	}
	
	
	/*   S e t t e r s   and   G e t t e r s  */
	//public LinkedHashMap<String, Track> getTrackSet_DEPRECATED() {
	//	return trackSet_DEPRECATED;
	//}
	
	public List<Track> getTrackList() {
		return trackList;
	}
	
	public List<Pattern> getRegexForTrackHeight() {
		return regexForTrackHeight;
	}
	
	public int getTrackHeightForRegex() {
		return trackHeightForRegex;
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
	
	/** Call the update method for all tracks in trackset. Updating can be time
	 * consuming especially for tracks associated to bam files. But be careful when 
	 * avoiding updating as it can lead to catastrophic out of sync data. 
	 * @throws SQLException 
	 * @throws InvalidRecordException 
	 * @throws InvalidGenomicCoordsException 
	 * @throws IOException 
	 * @throws ClassNotFoundException 
	 * @throws MalformedURLException */
	public void setGenomicCoordsAndUpdateTracks(GenomicCoords gc) throws MalformedURLException, ClassNotFoundException, IOException, InvalidGenomicCoordsException, InvalidRecordException, SQLException{

		for(Track tr : this.getTrackList()){
			tr.setGc(gc);
			tr.update();
		}
	}
	
	@Override
	/** For debugging and convenience only. This method not to be used for seriuous stuff. 
	 * */
	public String toString(){
		String x= "";
		for(Track tr : this.trackList){
			x += tr.toString() + "\n";
		}
		return x;
	}

	/** Iterate through track list and set regex to the track TrackSeqRegex.
	 * */
	public void setSeqRegexForTracks(String seqRegex) {
		
		for(Track tr : this.getTrackList()){
				tr.setSeqRegex(seqRegex);
		}
	}

}
