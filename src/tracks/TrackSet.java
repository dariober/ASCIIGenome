package tracks;

import java.io.File;
import java.io.IOException;
import java.net.MalformedURLException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.concurrent.TimeUnit;
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
// 	public static final String BOOKMARK_TAG= "bookmark";
	
	/*   C o n s t r u c t o r s   */
	
	public TrackSet(){}
	
	public TrackSet(List<String> inputFileList, GenomicCoords gc) throws IOException, InvalidGenomicCoordsException, InvalidRecordException, ClassNotFoundException, SQLException{
		
		for(String sourceName : inputFileList){
			try{
				if(Utils.getFileTypeFromName(sourceName).equals(TrackFormat.BAM)){
					//
					// BAM FILE
					//
					if(!Utils.bamHasIndex(sourceName)){
						System.err.println("\nNo index found for '" + sourceName + "'. Index can be generated with ");
						System.err.println("samtools index '" + sourceName + "'\n");
						throw new BamIndexNotFoundException();
					}
					
					/* Coverage track */
					TrackCoverage trackCoverage= new TrackCoverage(sourceName, gc, false);
	
					trackCoverage.setTrackTag(new File(sourceName).getName() + "#" + this.getNextTrackId());
					this.trackList.add(trackCoverage);
					
					/* Read track */
					TrackReads trackReads= new TrackReads(sourceName, gc);
					trackReads.setTrackTag(new File(sourceName).getName() + "@" + this.getNextTrackId());
					this.trackList.add(trackReads);
				}
				
				else if(    Utils.getFileTypeFromName(sourceName).equals(TrackFormat.BED) 
			        || Utils.getFileTypeFromName(sourceName).equals(TrackFormat.GFF)
				    || Utils.getFileTypeFromName(sourceName).equals(TrackFormat.VCF)){
					//
					// Annotatation
					//
					TrackIntervalFeature tif= new TrackIntervalFeature(sourceName, gc);
					this.add(tif, new File(sourceName).getName());
				} 
	
				else if(Utils.getFileTypeFromName(sourceName).equals(TrackFormat.BIGWIG) 
						|| Utils.getFileTypeFromName(sourceName).equals(TrackFormat.TDF) 
						|| Utils.getFileTypeFromName(sourceName).equals(TrackFormat.BEDGRAPH)){
					//
					// Wiggles
					//
					TrackWiggles tw= new TrackWiggles(sourceName, gc, 4);
					this.add(tw, new File(sourceName).getName());
					
				} else {
					// NB: You never get here because Utils.getFileTypeFromName returns
					// BED for any file that cannot be classified.
					System.err.println("Unable to classify " + sourceName + "; skipping"); 								
				}
			} catch(Exception e){
				System.err.println("Unable to classify " + sourceName + "; skipping");
				try {
					TimeUnit.SECONDS.sleep(3);
				} catch (InterruptedException e1) {
					e1.printStackTrace();
				}
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
		int idForTrack= this.getNextTrackId();
		String trackTag= baseTag + "#" + idForTrack;
		track.setTrackTag(trackTag);
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
	
	private void addWiggleTrackFromSourceName(String sourceName, GenomicCoords gc) throws IOException, InvalidRecordException, InvalidGenomicCoordsException, ClassNotFoundException, SQLException{
		
		int idForTrack= this.getNextTrackId();
		String trackId= new File(sourceName).getName() + "#" + idForTrack;
		
		TrackWiggles tw= new TrackWiggles(sourceName, gc, 4);
		tw.setTrackTag(trackId);
		this.trackList.add(tw);
	}
	
	private void addIntervalFeatureTrackFromSourceName(String sourceName, GenomicCoords gc) throws IOException, InvalidGenomicCoordsException, ClassNotFoundException, InvalidRecordException, SQLException{
		
		int idForTrack= this.getNextTrackId();
		String trackId= new File(sourceName).getName() + "#" + idForTrack;
		TrackIntervalFeature tif= new TrackIntervalFeature(sourceName, gc);
		tif.setTrackTag(trackId);
		this.trackList.add(tif);
	}
	
	private void addBamTrackFromSourceName(String sourceName, GenomicCoords gc) throws IOException, BamIndexNotFoundException, InvalidGenomicCoordsException, ClassNotFoundException, InvalidRecordException, SQLException{

		int idForTrack= this.getNextTrackId();
		
		if(!Utils.bamHasIndex(sourceName)){
			System.err.println("\nNo index found for '" + sourceName + "'. Index can be generated with ");
			System.err.println("samtools index '" + sourceName + "'\n");
			throw new BamIndexNotFoundException();
		}
		
		/* BAM Coverage track */
		String coverageTrackId= new File(sourceName).getName() + "#" + idForTrack;
		
		TrackCoverage trackCoverage= new TrackCoverage(sourceName, gc, false);
		trackCoverage.setTrackTag(coverageTrackId);
		this.trackList.add(trackCoverage);
		
		/* Reads */
		idForTrack= this.getNextTrackId();
		String trackId= new File(sourceName).getName() + "@" + idForTrack;

		TrackReads trackReads= new TrackReads(sourceName, gc);
		trackReads.setTrackTag(trackId);
		trackReads.setTrackTag(trackId);

		this.trackList.add(trackReads);
	} 
	
	private Integer getNextTrackId(){
		
		Integer id= 1;
		for(Track tr : this.getTrackList()){
			Integer x= this.getIdOfTrackName(tr.getTrackTag());
			if(x != null && x >= id){
				id= x + 1;
			}
		}
		return id; 
	}
	
	/** Oarse the track name to extract the ID suffix. Typically in the form
	 * name#1 or name@1 -> 1
	 * */
	private Integer getIdOfTrackName(String trackName){
		String x= trackName.trim().replaceAll(".*(@|#)", "");
		try{
			return Integer.parseInt(x);
		} catch (NumberFormatException e){
			return null;
		}
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
	 * @throws SQLException 
	 * @throws InvalidRecordException 
	 * @throws InvalidGenomicCoordsException 
	 * @throws IOException 
	 * @throws ClassNotFoundException 
	 * @throws MalformedURLException 
	*/
	public void setTrackHeightForRegex(List<String> tokens) throws InvalidCommandLineException, MalformedURLException, ClassNotFoundException, IOException, InvalidGenomicCoordsException, InvalidRecordException, SQLException{

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
	
	public void setPrintModeForRegex(List<String> cmdInput) throws InvalidCommandLineException {
		// MEMO of subcommand syntax:
		// 0 print/printFull
		// 1 Regex
		
		List<String> args= new ArrayList<String>();
		for(String x : cmdInput){
			args.add(x);
		}		
		
		boolean printFull= false;
		if(args.contains("-full")){
			printFull= true;
			args.remove("-full");
		}
		
        // Regex
        List<String> trackNameRegex= new ArrayList<String>();
        if(args.size() >= 2){
            trackNameRegex= args.subList(1, args.size());
        } else {
            trackNameRegex.add(".*"); // Default: Capture everything
        }
                
		PrintRawLine switchTo; // If printing is OFF do we switch to CLIP or FULL?
		if( printFull ){
			switchTo = PrintRawLine.FULL;
		} else {
			switchTo = PrintRawLine.CLIP;			
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
	
	
	public void setRpmForRegex(List<String> tokens) throws InvalidCommandLineException, MalformedURLException, ClassNotFoundException, IOException, InvalidGenomicCoordsException, InvalidRecordException, SQLException {

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
	 * @throws SQLException 
	 * @throws InvalidRecordException 
	 * @throws InvalidGenomicCoordsException 
	 * @throws IOException 
	 * @throws ClassNotFoundException 
	 * @throws MalformedURLException 
	*/
	public void setTrackYlimitsForRegex(List<String> tokens) throws InvalidCommandLineException, MalformedURLException, ClassNotFoundException, IOException, InvalidGenomicCoordsException, InvalidRecordException, SQLException{

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
	 * @throws SQLException 
	 * @throws InvalidRecordException 
	 * @throws InvalidGenomicCoordsException 
	 * @throws IOException 
	 * @throws ClassNotFoundException 
	*/
	public void setFilterForTrackIntervalFeature(List<String> cmdInput) throws InvalidCommandLineException, ClassNotFoundException, IOException, InvalidGenomicCoordsException, InvalidRecordException, SQLException{

		List<String> args= new ArrayList<String>();
		for(String x : cmdInput){
			args.add(x);
		}		
		args.remove(0); // Remove command name

		String showRegex= ".*";  // Show all
		String hideRegex= "^$";  // Hide nothing
		
		// Get args:
		if(args.contains("-i")){
			int idx= args.indexOf("-i") + 1; 
			showRegex= args.get(idx);
			args.remove(idx);
			args.remove("-i");
		}
		if(args.contains("-e")){
			int idx= args.indexOf("-e") + 1; 
			hideRegex= args.get(idx);
			args.remove(idx);
			args.remove("-e");
		}
		// What is left is positional args of regexes
		List<String> trackNameRegex= new ArrayList<String>();
		trackNameRegex.addAll(args);
		if(trackNameRegex.size() == 0){
			trackNameRegex.add(".*"); // Default regex for matching tracks
		}		
		
		// SHOW REGEX
		try{
			Pattern.compile(showRegex);
		} catch(PatternSyntaxException e){
	    	System.err.println("Invalid regex in: " + args);
	    	System.err.println("showRegex: " + showRegex);
	    	throw new InvalidCommandLineException();
		}
		
		// HIDE REGEX
		try{
			Pattern.compile(hideRegex); 
		} catch(PatternSyntaxException e){
	    	System.err.println("Invalid regex in: " + args);
	    	System.err.println("hideRegex: " + hideRegex);
	    	throw new InvalidCommandLineException();
	    }
		
		// TRACK REGEXES
        // Regex
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

		Track tr= this.matchIntervalFeatureTrack(trackId.trim());

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

		List<TrackIntervalFeature> ifTracks = this.getIntervalFeatureTracks();	
		
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

	private List<TrackIntervalFeature> getIntervalFeatureTracks(){
		
		// TrackSet ifSet= new TrackSet();
		List<TrackIntervalFeature> ifSet= new ArrayList<TrackIntervalFeature>();
		
		for(Track tr : this.getTrackList()){
			if(tr instanceof TrackIntervalFeature){
//			if(Utils.getFileTypeFromName(tr.getFilename()).equals(TrackFormat.BED) 
//			   || Utils.getFileTypeFromName(tr.getFilename()).equals(TrackFormat.GFF)
//			   || Utils.getFileTypeFromName(tr.getFilename()).equals(TrackFormat.VCF)){
				ifSet.add((TrackIntervalFeature) tr);
			}
		}
		return ifSet;
	}
	
	public void setDataColForRegex(List<String> tokens) throws InvalidCommandLineException, ClassNotFoundException, IOException, InvalidRecordException, InvalidGenomicCoordsException, SQLException {

		// MEMO of subcommand syntax:
		// 0 gffNameAttr
		// 1 attrName
		// 2 Regex

		int dataCol = 0; // Null will follow default 
		if(tokens.size() >= 2){
			try{
				dataCol= Integer.parseInt(tokens.get(1));
			} catch(NumberFormatException e){
				System.err.println("Number format exception: " + tokens.get(1));
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
	
	/** Stub: Add current position to bookmarks. Book
	 * @throws InvalidGenomicCoordsException 
	 * @throws SQLException 
	 * @throws InvalidRecordException 
	 * @throws IOException 
	 * @throws ClassNotFoundException 
	 * */
	public void addBookmark(GenomicCoords gc, String nameForBookmark) throws ClassNotFoundException, IOException, InvalidRecordException, SQLException, InvalidGenomicCoordsException{
		
		// Check there isn't a bookmark track already:
		boolean addbm = true;
		for(Track tr : this.getTrackList()){
			if(tr instanceof TrackBookmark){
 				tr.addBookmark(nameForBookmark);
 				tr.setGc(gc);
 				addbm= false;
				break;
			}
		}
		if(addbm){
			TrackBookmark tr= new TrackBookmark(gc, nameForBookmark);
			this.add(tr, "Bookmarks");
			tr.setGc(gc);
		}		
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
	
	/** Method to set any of the three alignment filters: -F,-f, mapq 
	 * @throws SQLException 
	 * @throws InvalidRecordException 
	 * @throws InvalidGenomicCoordsException 
	 * @throws IOException 
	 * @throws ClassNotFoundException 
	 * @throws MalformedURLException */
	public void setSamFilterForRegex(List<String> cmdInput) throws InvalidCommandLineException, MalformedURLException, ClassNotFoundException, IOException, InvalidGenomicCoordsException, InvalidRecordException, SQLException {
		// MEMO of subcommand syntax:
		// "-F INT -f INT -q INT re1 re2 ..."
		// * All optional. But -F, -f, -q are given they take one INT as arg.
		List<String> args= new ArrayList<String>();
		for(String x : cmdInput){
			args.add(x);
		}
		args.remove(0); // Remove name of command
		
		// Defaults:
		int f= 0;
		int F= 4; // Always remove unmapped
		int q= 0;
		
		// Get args:
		if(args.contains("-f")){
			int idx= args.indexOf("-f") + 1; 
			f= Integer.parseInt(args.get(idx));
			args.remove(idx);
			args.remove("-f");
		}
		if(args.contains("-F")){
			int idx= args.indexOf("-F") + 1; 
			F= Integer.parseInt(args.get(idx));
			if((4 & F) == 0){
				F += 4;
			}
			args.remove(idx);
			args.remove("-F");
		}
		if(args.contains("-q")){
			int idx= args.indexOf("-q") + 1; 
			q= Integer.parseInt(args.get(idx));
			args.remove(idx);
			args.remove("-q");
		}
		// What is left is positional args of regexes
		List<String> trackNameRegex= new ArrayList<String>();
		trackNameRegex.addAll(args);
		if(trackNameRegex.size() == 0){
			trackNameRegex.add(".*"); // Default regex for matching tracks
		}
		
        // Get tracks to reset:
        List<Track> tracksToReset = this.matchTracks(trackNameRegex, true);
        
        if(tracksToReset.size() == 0){
        	System.err.println("No track matching regex " + trackNameRegex);
        }
        
        for(Track tr : tracksToReset){
            
        	tr.set_f_flag(f);
        	tr.set_F_flag(F);
        	tr.setMapq(q);
        	
        	List<SamRecordFilter> filters= new ArrayList<SamRecordFilter>();
        	filters.addAll(FlagToFilter.flagToFilterList(f, F));
        	filters.add(new MappingQualityFilter(q));
        	tr.setSamRecordFilter(filters);
        }		
	}
	
	public void editNamesForRegex(List<String> cmdInput) throws InvalidCommandLineException {
		// API:
		// editNames <patterns> <replacement> [track_re] ...
		List<String> args= new ArrayList<String>();
		for(String x : cmdInput){
			args.add(x);
		}
		args.remove(0); // Remove name of command
		
		if(args.size() < 2){
			System.err.println("Expected at least two arguments. Got " + cmdInput);
			throw new InvalidCommandLineException();
		}
		
		String pattern= args.get(0); args.remove(0);
		String replacement= args.get(0); args.remove(0);
		if(replacement.equals("\"\"")){
			replacement= "";
		}
		if(args.size() == 0){ // What is left is pos args. Set to .* if nothing left.
			args.add(".*");
		}
		
        // Get tracks to reset:
        List<Track> tracksToReset = this.matchTracks(args, true);
        
        // Dry-run: See if it duplicate or empty names are generated: 
        // These are the track tags that would be generated. Including the edited ones and the untouched ones. 
        List<String> testNames= new ArrayList<String>(); 
        for(Track tr : this.getTrackList()){
        	if(tracksToReset.contains(tr)){
	        	String newTag= tr.getTrackTag().replaceAll(pattern, replacement);
	        	if(newTag.isEmpty()){
	    			System.err.println("Empty tag is not allowed. Empty for " + tr.getTrackTag());
	    			throw new InvalidCommandLineException();
	        	}
	        	testNames.add(newTag);
        	} else {
        		testNames.add(tr.getTrackTag());
        	}
        }
        Set<String> set = new HashSet<String>(testNames);
        if(set.size() < testNames.size() ){
        	System.err.println("Edit not allowed: It generates duplicate tags.");
        	throw new InvalidCommandLineException();
        }
        
        // Now change for real
        for(Track tr : tracksToReset){
        	String newTag= tr.getTrackTag().replaceAll(pattern, replacement);
        	tr.setTrackTag(newTag);
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
			// tr.update();
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
	 * @throws SQLException 
	 * @throws InvalidRecordException 
	 * @throws InvalidGenomicCoordsException 
	 * @throws IOException 
	 * @throws ClassNotFoundException 
	 * */
	public void setSeqRegexForTracks(List<String> cmdInput) throws ClassNotFoundException, IOException, InvalidGenomicCoordsException, InvalidRecordException, SQLException {

		List<String> args= new ArrayList<String>();
		for(String x : cmdInput){
			args.add(x);
		}		
		args.remove(0); // Remove command name
		
		boolean isCaseSensisitive= false;
		if(args.contains("-c")){
			isCaseSensisitive= true;
			args.remove("-c");
		}
		
		String seqRegex= null;
		if(args.size() == 0){
			seqRegex= "";
		} else {
			seqRegex= args.get(0);
			try{
				Pattern.compile(seqRegex);
			} catch(PatternSyntaxException e){
		    	System.err.println("Invalid seqRegex in: " + args);
		    	System.err.println(e.getDescription());
			}
		}

		for(Track tr : this.getTrackList()){
			if(tr instanceof TrackSeqRegex){
				((TrackSeqRegex) tr).setCaseSensitive(isCaseSensisitive);
				tr.setSeqRegex(seqRegex);
			}	
		}
	}

	public void dropTracksWithRegex(List<String> cmdInput) throws InvalidCommandLineException {

		List<String> trackNameRegex= cmdInput.subList(1, cmdInput.size());
       
        List<Track> tracksToDrop= this.matchTracks(trackNameRegex, true);
        for(Track tr : tracksToDrop){
        	boolean removed= this.trackList.remove(tr);
        	if(removed){
        		System.err.println("Dropped: " + tr.getTrackTag());
        	}
        }
        
	}

	
//	/** Attempt to collect source of sequence dictionary 
//	 * */
//	private String getSamSeqDictSource(){
//		String samSeqDictSource= null;
//		// First try to find fasta
//		for(Track tr : this.getTrackList()){
//			if(tr.getGc().getFastaFile() != null){
//				return tr.getGc().getFastaFile();
//			}
//		}
//		// Try directly the get method:
//		for(Track tr : this.getTrackList()){
//			if(tr.getGc().getSamSeqDictSource() != null){
//				return tr.getGc().getSamSeqDictSource();
//			}
//		}
//		return samSeqDictSource;
//	}
	
}
