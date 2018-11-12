package tracks;

import java.io.File;
import java.io.IOException;
import java.net.MalformedURLException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.TimeUnit;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.regex.PatternSyntaxException;

import com.google.common.base.Joiner;
import com.google.common.base.Splitter;
import com.google.common.collect.Lists;

import coloring.Xterm256;
import exceptions.BamIndexNotFoundException;
import exceptions.InvalidColourException;
import exceptions.InvalidCommandLineException;
import exceptions.InvalidGenomicCoordsException;
import exceptions.InvalidRecordException;
import filter.FlagToFilter;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.filter.MappingQualityFilter;
import htsjdk.samtools.filter.SamRecordFilter;
import net.sourceforge.argparse4j.inf.ArgumentParserException;
import samTextViewer.GenomicCoords;
import samTextViewer.Utils;

/** Class to hold tracks to be printed. 
 * */
public class TrackSet {
	
	// private LinkedHashMap<String, Track> trackSet_DEPRECATED= new LinkedHashMap<String, Track>();
	private List<Track> trackList= new ArrayList<Track>();
	private List<Pattern> regexForTrackHeight= new ArrayList<Pattern>();
	private int trackHeightForRegex= -1;
	private String[] yStrLimits= new String[2];
	// private List<String> regexForYLimits= new ArrayList<String>(); 
	private List<Track> tracksForYLimits= new ArrayList<Track>();
	private LinkedHashSet<String> openedFiles= new LinkedHashSet<String>();
	
	/*   C o n s t r u c t o r s   */
	
	public TrackSet(List<String> inputFileList, GenomicCoords gc) throws IOException, InvalidGenomicCoordsException, InvalidRecordException, ClassNotFoundException, SQLException{
		
		for(String sourceName : inputFileList){
			try{
				if(Utils.getFileTypeFromName(sourceName).equals(TrackFormat.BAM)){
					/* Coverage track */
					TrackPileup trackPileup= new TrackPileup(sourceName, gc);
					// trackPileup.setTrackTag(new File(sourceName).getName() + "#" + this.getNextTrackId());
					trackPileup.setTrackTag(sourceName + "#" + this.getNextTrackId());
					this.trackList.add(trackPileup);
					
					/* Read track */
					TrackReads trackReads= new TrackReads(sourceName, gc);
					// trackReads.setTrackTag(new File(sourceName).getName() + "@" + this.getNextTrackId());
					trackReads.setTrackTag(sourceName + "@" + this.getNextTrackId());
					this.trackList.add(trackReads);
				}
				
				else if(Utils.getFileTypeFromName(sourceName).equals(TrackFormat.BED) 
			        || Utils.getFileTypeFromName(sourceName).equals(TrackFormat.GFF)
			        || Utils.getFileTypeFromName(sourceName).equals(TrackFormat.GTF)
				    || Utils.getFileTypeFromName(sourceName).equals(TrackFormat.VCF)
				    || Utils.getFileTypeFromName(sourceName).equals(TrackFormat.BIGBED)){
					//
					// Annotatation
					//
					TrackIntervalFeature tif= new TrackIntervalFeature(sourceName, gc);
					// this.addTrack(tif, new File(sourceName).getName());
					this.addTrack(tif, sourceName);
				} 
	
				else if(Utils.getFileTypeFromName(sourceName).equals(TrackFormat.BIGWIG) 
						|| Utils.getFileTypeFromName(sourceName).equals(TrackFormat.TDF) 
						|| Utils.getFileTypeFromName(sourceName).equals(TrackFormat.BEDGRAPH)){
					//
					// Wiggles
					//
					TrackWiggles tw= new TrackWiggles(sourceName, gc, 4);
					// this.addTrack(tw, new File(sourceName).getName());
					this.addTrack(tw, sourceName);
					
				} else {
					// NB: You never get here because Utils.getFileTypeFromName returns
					// BED for any file that cannot be classified.
					System.err.println("Unable to classify " + sourceName + "; skipping"); 								
				}
			} catch(Exception e){
				System.err.println(e.getMessage());
				e.printStackTrace();
				System.err.println("Cannot add " + sourceName + "; skipping");
				System.exit(1);
			}
		}	
		this.tracksForYLimits.addAll(this.getTrackList());
		this.yStrLimits[0]= "na";
		this.yStrLimits[1]= "na";
		// TrackWiggles gcProfile= gc.getGCProfile();
		for(Track tr : this.getTrackList()){
			this.addToOpenedFiles(tr.getFilename());
		}
		
		Runtime.getRuntime().addShutdownHook(new ShutDownTask(this));
		
	}
	
	private class ShutDownTask extends Thread{
		
		private TrackSet trackSet;

		public ShutDownTask(TrackSet trackSet){
			this.trackSet= trackSet;
		}

		@Override
		public void run(){
			for(Track tr : trackSet.getTrackList()){
				tr.close();
			}
		}
		
	}

	
	/*   M e t h o d s   */

	/**
	 * Add this track with given baseTag. The suffix "#id" will be appended to the baseTag string. 
	 * NB: Adding a track resets the track's ID
	 * */
	public void addTrack(Track track, String baseTag) {
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
	public void addTrackFromSource(String sourceName, GenomicCoords gc, String trackTag) throws IOException, BamIndexNotFoundException, InvalidGenomicCoordsException, InvalidRecordException, ClassNotFoundException, SQLException{

		if(Utils.getFileTypeFromName(sourceName).equals(TrackFormat.BAM)){
			this.addBamTrackFromSourceName(sourceName, gc, trackTag);
		
		} else if( Utils.getFileTypeFromName(sourceName).equals(TrackFormat.VCF) ){
			this.addIntervalFeatureTrackFromVCF(sourceName, gc, trackTag);
			
		} else if(Utils.getFileTypeFromName(sourceName).equals(TrackFormat.BED) 
				  || Utils.getFileTypeFromName(sourceName).equals(TrackFormat.BIGBED)
		          || Utils.getFileTypeFromName(sourceName).equals(TrackFormat.GFF)
		          || Utils.getFileTypeFromName(sourceName).equals(TrackFormat.GTF)){
			this.addIntervalFeatureTrackFromSourceName(sourceName, gc, trackTag);
		
		} else if(Utils.getFileTypeFromName(sourceName).equals(TrackFormat.BIGWIG) 
				|| Utils.getFileTypeFromName(sourceName).equals(TrackFormat.TDF) 
				|| Utils.getFileTypeFromName(sourceName).equals(TrackFormat.BEDGRAPH)){
			this.addWiggleTrackFromSourceName(sourceName, gc, trackTag);
		} else {
			System.err.println("Unexpected file format: " + Utils.getFileTypeFromName(sourceName) + " for " + sourceName);
			throw new RuntimeException();
		}
		
		for(Track tr : this.getTrackList()){
			this.addToOpenedFiles(tr.getFilename());
		}			
	}

	private void addToOpenedFiles(String sourceName){
		if(this.getOpenedFiles().contains(sourceName)){ // Remove and add as last opened
			this.openedFiles.remove(sourceName);
		} 
		this.openedFiles.add(sourceName);
	}
	
	private void addWiggleTrackFromSourceName(String sourceName, GenomicCoords gc, String trackTag) throws IOException, InvalidRecordException, InvalidGenomicCoordsException, ClassNotFoundException, SQLException{
		
		int idForTrack= this.getNextTrackId();
		// String trackId= new File(sourceName).getName() + "#" + idForTrack;
		String trackId= sourceName + "#" + idForTrack;
		
		TrackWiggles tw= new TrackWiggles(sourceName, gc, 4);
		tw.setTrackTag(trackId);
		this.trackList.add(tw);
	}
	
	private void addIntervalFeatureTrackFromSourceName(String sourceName, GenomicCoords gc, String trackTag) throws IOException, InvalidGenomicCoordsException, ClassNotFoundException, InvalidRecordException, SQLException{
		
		int idForTrack= this.getNextTrackId();
		// String trackId= new File(sourceName).getName() + "#" + idForTrack;
		String trackId= sourceName + "#" + idForTrack;
		TrackIntervalFeature tif= new TrackIntervalFeature(sourceName, gc);
		tif.setTrackTag(trackId);
		this.trackList.add(tif);
	}
	
	private void addIntervalFeatureTrackFromVCF(String sourceName, GenomicCoords gc, String trackTag) throws ClassNotFoundException, IOException, InvalidGenomicCoordsException, InvalidRecordException, SQLException {
		
		int idForTrack= this.getNextTrackId();
		// String trackId= new File(sourceName).getName() + "#" + idForTrack;
		String trackId= sourceName + "#" + idForTrack;
		
		// If this VCF file has sequence dictionary, check the coordinates in gc are compatible
		// If they are not, throw an exception which force resetting the coords.
		SAMSequenceDictionary seqDict = Utils.getVCFHeader(sourceName).getSequenceDictionary();
		
		if(seqDict != null && seqDict.getSequence(gc.getChrom()) == null){
			throw new InvalidGenomicCoordsException();
		}
		TrackIntervalFeature tif= new TrackIntervalFeature(sourceName, gc);
		tif.setTrackTag(trackId);
		this.trackList.add(tif);
	}

	
	private void addBamTrackFromSourceName(String sourceName, GenomicCoords gc, String trackTag) throws IOException, BamIndexNotFoundException, InvalidGenomicCoordsException, ClassNotFoundException, InvalidRecordException, SQLException{

		/* BAM Coverage track */
		int idForTrack= this.getNextTrackId();
		// String coverageTrackId= new File(sourceName).getName() + "#" + idForTrack;
		String coverageTrackId= sourceName + "#" + idForTrack;
		
		TrackPileup trackPileup= new TrackPileup(sourceName, gc);
		trackPileup.setTrackTag(coverageTrackId);
		this.trackList.add(trackPileup);
		
		/* Reads */
		idForTrack= this.getNextTrackId();
		// String trackId= new File(sourceName).getName() + "@" + idForTrack;
		String trackId= sourceName + "@" + idForTrack;

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
			trackInfo.add(
					  "------\n"
					+ "Track tag:    " + track.getTrackTag() + "\n" 
					+ "Input source: " + track.getFilename() + "\n"
					+ "Working file: " + track.getWorkFilename() + "\n"
					+ "Track type:   " + Utils.getFileTypeFromName(track.getFilename()) + " " + hd + "\n"
					+ "awk script:   " + (! track.getAwk().trim().isEmpty() ? track.getAwk() : "N/A") + "\n"
					+ "grep:         " + "show: " + track.getShowRegex() + "; hide: " + track.getHideRegex());
		}
		
		StringBuilder sb= new StringBuilder();
		for(String str : trackInfo){
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
		
        boolean invertSelection= Utils.argListContainsFlag(tokens, "-v");
		
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
        List<Track> tracksToReset = this.matchTracks(trackNameRegex, true, invertSelection);
        for(Track tr : tracksToReset){
        	tr.setyMaxLines(this.trackHeightForRegex);
        }
	}

	public void setFeatureDisplayModeForRegex(List<String> tokens) throws InvalidCommandLineException {

		List<String> args= new ArrayList<String>(tokens);
		args.remove(0); // remove command name

		// Display mode
		FeatureDisplayMode mode= null;
		if(args.contains("-collapsed")){
			mode= FeatureDisplayMode.COLLAPSED;
			args.remove("-collapsed");
		}
		if(args.contains("-c")){
			mode= FeatureDisplayMode.COLLAPSED;
			args.remove("-c");
		}
		
		if(args.contains("-expanded")){
			mode= FeatureDisplayMode.EXPANDED;
			args.remove("-expanded");
		}
		if(args.contains("-e")){
			mode= FeatureDisplayMode.EXPANDED;
			args.remove("-e");
		}
		
		if(args.contains("-oneline")){
			mode= FeatureDisplayMode.ONELINE;
			args.remove("-oneline");
		}
		if(args.contains("-o")){
			mode= FeatureDisplayMode.ONELINE;
			args.remove("-o");
		}
		boolean invertSelection= Utils.argListContainsFlag(args, "-v");
		
        // Regex to capture tracks: Everything left after removing command name and args:
        List<String> trackNameRegex= new ArrayList<String>();
        if(args.size() > 0){
            trackNameRegex.addAll(args);
        } else {
            trackNameRegex.add(".*"); // Default: Capture everything
        }
                
        // And set as required:
		List<Track> tracksToReset = this.matchTracks(trackNameRegex, true, invertSelection);
        for(Track tr : tracksToReset){
        	if(mode == null){ // Toggle between expanded and collapsed
        		if(tr.getFeatureDisplayMode().equals(FeatureDisplayMode.EXPANDED)){
        			tr.setFeatureDisplayMode(FeatureDisplayMode.COLLAPSED);
        		} else if(tr.getFeatureDisplayMode().equals(FeatureDisplayMode.COLLAPSED)){
        			tr.setFeatureDisplayMode(FeatureDisplayMode.EXPANDED);
        		} else if(tr.getFeatureDisplayMode().equals(FeatureDisplayMode.ONELINE)){
        			tr.setFeatureDisplayMode(FeatureDisplayMode.EXPANDED);
        		}
        	} else {
        		tr.setFeatureDisplayMode(mode);
        	}
        }
	}
	
	public void setPrintModeAndPrintFeaturesForRegex(List<String> cmdInput) throws InvalidCommandLineException, IOException, InvalidGenomicCoordsException, InvalidColourException, ArgumentParserException {

		// --------------------------------------------------------------------
		// PARSE ARGUMENTS
		List<String> args= new ArrayList<String>(cmdInput);
		args.remove(0); // Remove cmd name.
		
		// Defaults if no args given
		PrintRawLine printMode = null;
        int count= 10; // opts.getInt("nlines"); // Default number of lines to print
        List<String> trackNameRegex= new ArrayList<String>(); // opts.getList("track_regex");
        trackNameRegex.add(".*");
        
        String printToFile= null;
        boolean append= false;
        
		// PARSE ARGS
        // --------------------------------------------------------------------
        boolean invertSelection= Utils.argListContainsFlag(args, "-v");

        boolean esf= Utils.argListContainsFlag(args, "-esf");
        if(esf){
        	printMode= PrintRawLine.CLIP;
        }
        
		String vep= Utils.getArgForParam(args, "-vep", null);
        if(vep != null){
        	printMode= PrintRawLine.FULL;
        }
        
        Pattern highlightPattern= Pattern.compile("");
        if(args.contains("-hl")){
        	if(printMode == null){
        		printMode= PrintRawLine.CLIP;
        	}
    		String p= Utils.getArgForParam(args, "-hl", "");
	        highlightPattern= Pattern.compile(p);
        }
        
        if(args.contains("-clip")){
			printMode= PrintRawLine.CLIP;
			args.remove("-clip");
		}
        
        if(args.contains("-full")){
			printMode= PrintRawLine.FULL;
			args.remove("-full");
		}

		if(args.contains("-off")){
			printMode= PrintRawLine.OFF;
			args.remove("-off");
		}
		
		if(args.contains("-n")){
			printMode= PrintRawLine.NO_ACTION;
			int idx= args.indexOf("-n") + 1; 
			try{
				count= Integer.parseInt(args.get(idx));
				args.remove(idx);
				args.remove("-n");
			} catch (NumberFormatException e){
				System.err.println("Invalid argument to -n. Expected INT. Got '" + args.get(idx) + "'");
				throw new InvalidCommandLineException();
			}
		}
		
		String sys= Utils.getArgForParam(args, "-sys", "");
		if(sys != ""){
			printMode= PrintRawLine.NO_ACTION;
		}
		
    	String round= Utils.getArgForParam(args, "-round", "");
    	Integer printNumDecimals= null;
    	if( ! round.isEmpty()){
    		printMode= PrintRawLine.NO_ACTION;
            try{
            	printNumDecimals= Integer.valueOf(round);
    		} catch (NumberFormatException e){
    			System.err.println("Invalid argument to -round. An integer is expected.");
    			throw new InvalidCommandLineException();
    		}    		
    	}        
		
		// Capture the redirection operator and remove operator and filename.
		int idx= -1; // Position of the redirection operator, -1 if not present.
		if(args.contains(">")){
			idx= args.indexOf(">");
		} else if(args.contains(">>")){
			idx= args.indexOf(">>");
			append= true;
		}
		if(idx >= 0){ // Redirection found, write to file, unless file is not given.
			try{
				printToFile= args.get(idx+1);
				args.remove(idx);
				args.remove(printToFile);
				printToFile= Utils.tildeToHomeDir(printToFile);
			} catch(IndexOutOfBoundsException e){
				System.err.println("No file found to write to.");
				throw new InvalidCommandLineException();
			}
		}

        // Everything left in arg list is positional args of regex for tracks
		if(args.size() >= 1){
            trackNameRegex= args;
		}
		// END PARSING ARGS. 
		// --------------------------------------------------------------------
		
        // Tracks affected by this command:
        List<Track> tracksToReset = this.matchTracks(trackNameRegex, true, invertSelection);

        // If we print to file just do that.
		if(printToFile != null && ! printToFile.isEmpty()){

	        if( ! append){
	        	new File(printToFile).delete();
	        }

			for(Track tr : tracksToReset){
	        	tr.setExportFile(printToFile);
	        	tr.setHighlightPattern(highlightPattern);
	        	tr.printLines();
	        	tr.setPrintMode(PrintRawLine.OFF); // This is not ideal: redirecting to file also set mode to off
	        	tr.setExportFile(null); // Reset to null so we don't keep writing to this file once we go to another position. 
	        }
	        return;
		}

		// Process as required
		for(Track tr : tracksToReset){
			tr.setExplainSamFlag(esf);
			tr.setPrintFormattedVep(vep);
        	tr.setHighlightPattern(highlightPattern);
			tr.setPrintRawLineCount(count);
			tr.setSystemCommandForPrint(sys);
			if(printNumDecimals != null) tr.setPrintNumDecimals(printNumDecimals);
			if(printMode != null && printMode.equals(PrintRawLine.NO_ACTION)){
				if(tr.getPrintMode().equals(PrintRawLine.OFF) && (cmdInput.contains("-n") || cmdInput.contains("-sys") || cmdInput.contains("-round"))){
					// Make -n or -sys switch ON the printing mode
					// This happens if you exec `print -n INT` with the track set to OFF. 
					tr.setPrintMode(PrintRawLine.CLIP); 
				}
				// 
			} else if(printMode != null){
				tr.setPrintMode(printMode);
				
			} else if(tr.getPrintMode().equals(PrintRawLine.OFF)) { // Toggle
				tr.setPrintMode(PrintRawLine.CLIP);
				
			} else {
				tr.setPrintMode(PrintRawLine.OFF);
			}
			//tr.setCutScriptForPrinting(cutScript);
        }        
	}
	
	public void setBisulfiteModeForRegex(List<String> tokens) throws InvalidCommandLineException {

		List<String> args= new ArrayList<String>(tokens);
		args.remove(0);
		
        boolean invertSelection= Utils.argListContainsFlag(args, "-v");
		
		Boolean bisulf= null;
		if(args.contains("-on")){
			bisulf= true;
			args.remove("-on");
		}
		if(args.contains("-off")){
			bisulf= false;
			args.remove("-off");
		}
		
        // Regex
        List<String> trackNameRegex= new ArrayList<String>();
        if(args.size() > 0){
            trackNameRegex= args;
        } else {
            trackNameRegex.add(".*"); // Default: Capture everything
        }

        // And set as required:
        List<Track> tracksToReset = this.matchTracks(trackNameRegex, true, invertSelection);
        for(Track tr : tracksToReset){
        	if(bisulf == null){
    			if(tr.isBisulf()){ // Invert setting
    				tr.setBisulf(false);
    			} else {
    				tr.setBisulf(true);
    			}
        	} else {
        		tr.setBisulf(bisulf);
        	}
        }
	}

	public void setHideTitleForRegex(List<String> tokens) throws InvalidCommandLineException {

		List<String> args= new ArrayList<String>(tokens);
		args.remove(0);
		
        boolean invertSelection= Utils.argListContainsFlag(args, "-v");
		
		Boolean hide= null;
		if(args.contains("-on")){
			hide= true;
			args.remove("-on");
		}
		if(args.contains("-off")){
			hide= false;
			args.remove("-off");
		}
		
        // Regex
        List<String> trackNameRegex= new ArrayList<String>();
        if(args.size() > 0){
            trackNameRegex= args;
        } else {
            trackNameRegex.add(".*"); // Default: Capture everything
        }

        // And set as required:
        List<Track> tracksToReset = this.matchTracks(trackNameRegex, true, invertSelection);
        for(Track tr : tracksToReset){
        	if(hide == null){
    			if(tr.isHideTitle()){ // Invert setting
    				tr.setHideTitle(false);
    			} else {
    				tr.setHideTitle(true);
    			}
        	} else {
        		tr.setHideTitle(hide);
        	}
        }
	}
		
	public void setRpmForRegex(List<String> tokens) throws InvalidCommandLineException, MalformedURLException, ClassNotFoundException, IOException, InvalidGenomicCoordsException, InvalidRecordException, SQLException {

		List<String> args= new ArrayList<String>(tokens);
		args.remove(0);
		
        boolean invertSelection= Utils.argListContainsFlag(args, "-v");
		
		Boolean rpm= null;
		if(args.contains("-on")){
			rpm= true;
			args.remove("-on");
		}
		if(args.contains("-off")){
			rpm= false;
			args.remove("-off");
		}
		
        // Regex
        List<String> trackNameRegex= new ArrayList<String>();
        if(args.size() > 0){
            trackNameRegex= args;
        } else {
            trackNameRegex.add(".*"); // Default: Capture everything
        }

        // And set as required:
        List<Track> tracksToReset = this.matchTracks(trackNameRegex, true, invertSelection);
        for(Track tr : tracksToReset){
        	if(rpm == null){
    			if(tr.isRpm()){ // Invert setting
    				tr.setRpm(false);
    			} else {
    				tr.setRpm(true);
    			}
        	} else {
        		tr.setRpm(rpm);
        	}
        }
	}

	public void setGenotypeMatrix(List<String> cmdTokens) throws InvalidCommandLineException {
		List<String> argList= new ArrayList<String>(cmdTokens);
		argList.remove(0);

		// Collect arguments
		String nRows= Utils.getArgForParam(argList, "-n", null);
		String selectSampleRegex= Utils.getArgForParam(argList, "-s", null);
		List<String> subSampleRegex= Utils.getNArgsForParam(argList, "-r", 2);
		String jsScriptFilter= Utils.getArgForParam(argList, "-f", null);
		
		boolean invertSelection= Utils.argListContainsFlag(argList, "-v");
		
		// Regex to capture tracks: All positional args left
        List<String> trackNameRegex= new ArrayList<String>();
        if(argList.size() > 0){
            trackNameRegex= argList;
        } else {
            trackNameRegex.add(".*"); // Default: Capture everything
        }
        
        // Set as appropriate
        List<Track> tracksToReset = this.matchTracks(trackNameRegex, true, invertSelection);
        for(Track tr : tracksToReset){
        	if(nRows != null){
        		int n= Integer.valueOf(nRows) >= 0 ? Integer.valueOf(nRows) : Integer.MAX_VALUE;
        		tr.getGenotypeMatrix().setnMaxSamples(n);
        	}
        	if(selectSampleRegex != null){
            	tr.getGenotypeMatrix().setSelectSampleRegex(selectSampleRegex);
        	}
        	if(subSampleRegex != null){
        		String rplc= subSampleRegex.get(1);
        		rplc= ! rplc.equals("\"\"") ? rplc : ""; 
            	tr.getGenotypeMatrix().setSubSampleRegex(subSampleRegex.get(0), rplc);
        	}
        	if(jsScriptFilter != null){
            	tr.getGenotypeMatrix().setJsScriptFilter(jsScriptFilter);
        	}
        }

	}	
	
	public void setFeatureColorForRegex(List<String> cmdTokens) throws InvalidCommandLineException, InvalidColourException {

		List<String> argList= new ArrayList<String>(cmdTokens);
		argList.remove(0); // Remove cmd name
		
		// Collect all regex/color pairs from input. We move left to right along the command 
		// arguments and collect -r/-R and set the regex inversion accordingly.
		List<Argument> colorForRegex= new ArrayList<Argument>();
		new Xterm256();
		while(argList.contains("-r") || argList.contains("-R")){
			int r= argList.indexOf("-r") >= 0 ? argList.indexOf("-r") : Integer.MAX_VALUE;
			int R= argList.indexOf("-R") >= 0 ? argList.indexOf("-R") : Integer.MAX_VALUE;
			List<String> pair;
			boolean invert= false;
			if(r < R){
				pair = Utils.getNArgsForParam(argList, "-r", 2);
			} else {
				pair = Utils.getNArgsForParam(argList, "-R", 2);
				invert= true;
			}
			String pattern= pair.get(0);
			try{ // Check valid regex
				Pattern.compile(pattern); 
			} catch(PatternSyntaxException e){
		    	System.err.println("Invalid regex: " + pattern);
		    	throw new InvalidCommandLineException();
			}
			Argument xcolor= new Argument(pair.get(0), pair.get(1), invert);
			Xterm256.colorNameToXterm256(xcolor.getArg()); // Check this is a valid colour 
			colorForRegex.add(xcolor);
		}
		if(colorForRegex.size() == 0){
			colorForRegex= null;
		}
		
		boolean invertSelection= Utils.argListContainsFlag(argList, "-v");
		
		// Regex to capture tracks: All positional args left
        List<String> trackNameRegex= new ArrayList<String>();
        if(argList.size() > 0){
            trackNameRegex= argList;
        } else {
            trackNameRegex.add(".*"); // Default: Capture everything
        }
        
        // Set as appropriate
        List<Track> tracksToReset = this.matchTracks(trackNameRegex, true, invertSelection);
        for(Track tr : tracksToReset){
        	tr.setColorForRegex(colorForRegex);
        }
	}
	
	public void setTrackColourForRegex(List<String> tokens) throws InvalidCommandLineException, InvalidColourException{

		// MEMO of subcommand syntax:
		// 0 trackColour
		// 1 Colour
		// 2 Regex

		boolean invertSelection= Utils.argListContainsFlag(tokens, "-v");
		
		// Colour
		String colour= (new TrackIntervalFeature(null)).getTitleColour();
		if(tokens.size() >= 2){
			String xcolour= tokens.get(1).toLowerCase();

			Xterm256.colorNameToXterm256(xcolour); // This is only to test whether exception is thrown.

			colour= xcolour;
		}
		
		// Regex
		List<String> trackNameRegex= new ArrayList<String>();
		if(tokens.size() >= 3){
			trackNameRegex= tokens.subList(2, tokens.size());
		} else {
			trackNameRegex.add(".*"); // Default: Capture everything
		}
		// And set as required:
		List<Track> tracksToReset = this.matchTracks(trackNameRegex, true, invertSelection);
		for(Track tr : tracksToReset){
			tr.setTitleColour(colour);
		}
	}
	
	public void setAttributeForGFFName(List<String> tokens) throws InvalidCommandLineException{

		// MEMO of subcommand syntax:
		// 0 gffNameAttr
		// 1 attrName
		// 2 Regex

        boolean invertSelection= Utils.argListContainsFlag(tokens, "-v");
		
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
        List<Track> tracksToReset = this.matchTracks(trackNameRegex, true, invertSelection);
        for(Track tr : tracksToReset){
        	tr.setGtfAttributeForName(gtfAttributeForName);
        }		
	}
	
	/** Set ylimits for all tracks using private attributes of ylim and track regex.
	 * @throws InvalidCommandLineException 
	 * */
	public void setAutoYLimits() throws InvalidCommandLineException {

		// Tracks to reset.
		// To be depracated. trackToReset should be replaced by the content of field this.tracksForYlimits;
        List<Track> tracksToReset = this.getTracksForYLimits(); // this.matchTracks(this.getRegexForYLimits(), true, false);
		
        String yStrMin= this.getYStringLimits()[0];
        String yStrMax= this.getYStringLimits()[1];
        
		Float[] yrange= {Float.NaN, Float.NaN};
		if(yStrMin.equals("min") || yStrMax.equals("max")){
			yrange= this.yRangeOfTracks(tracksToReset);
		}

		// Parse min
		double ymin= Double.NaN;
		if(yStrMin.equals("min")){
			ymin= yrange[0];
		} else {
			try{
				ymin= Double.parseDouble(yStrMin);
			} catch(NumberFormatException e){
				ymin= Double.NaN;
			}
		}

		// Parse max
		double ymax= Double.NaN;
		if(yStrMax.equals("max")){
			ymax= yrange[1];
		} else { 
			try{
				ymax= Double.parseDouble(yStrMax);
			} catch(NumberFormatException e){
				ymax= Double.NaN;
			}
		}
		// Swap
		if(ymin > ymax){
			Double newMax= ymin;
			ymin= ymax;
			ymax= newMax;			
		}

		String[] yy = Utils.roundToSignificantDigits(ymin, ymax, 2);
		
        for(Track tr : tracksToReset){
    		tr.setYLimitMin(Float.valueOf(yy[0]));
			tr.setYLimitMax(Float.valueOf(yy[1]));
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

		List<String> args= new ArrayList<String>(tokens);
		
		boolean invertSelection= Utils.argListContainsFlag(args, "-v");
		
		if(args.size() < 3){
			System.err.println("Error in ylim subcommand. Expected at least 2 positional arguments args got: " + args);
			throw new InvalidCommandLineException();
		}

		// Set list of tracks to be reset for ylim
        if(args.size() >= 4){
        	this.setTracksForYLimits(this.matchTracks(args.subList(3, args.size()), true, invertSelection));
        } else {
        	this.setTracksForYLimits(this.matchTracks(Arrays.asList(".*"), true, invertSelection));
        }
		// ----------------------------------
//        if(args.size() >= 4){
//        	this.setRegexForYLimits(args.subList(3, tokens.size()));
//        } else {
//        	this.setRegexForYLimits(Arrays.asList(".*"));
//        }
        // ----------------------------------
        
        // Set ylimits attribute
        String[] yStrLimits= {args.get(1), args.get(2)};
        this.setYStrLimits(yStrLimits);
        
        // Reset ylimits as required
        this.setAutoYLimits();
	}

	/** Get the range of all the screen scores of this list of tracks. I.e. the global min and max.
	 * */
	private Float[] yRangeOfTracks(List<Track> tracks){
		
		List<Float> yall= new ArrayList<Float>();
		for(Track tr : tracks){
			yall.addAll(tr.getScreenScores());
		}		
		return Utils.range(yall);
	}
	
	/** Parse awk command and set awk script for the captured tracks. 
	 * The tokens in cmdInput maybe in the form:
	 * [awk, $3 > 10]
	 * [awk, $3 > 10, track1]
	 * [awk, $3 > 10, track_re1, track_re2]
	 * [awk, -F, sep, -v, n=10, $3 > n, track1, track2]
	 * @throws InvalidCommandLineException 
	 * @throws SQLException 
	 * @throws InvalidRecordException 
	 * @throws InvalidGenomicCoordsException 
	 * @throws IOException 
	 * @throws ClassNotFoundException 
	 * */
	public void setAwkForTrack(List<String> cmdInput) throws InvalidCommandLineException, ClassNotFoundException, IOException, InvalidGenomicCoordsException, InvalidRecordException, SQLException {

		List<String> args= new ArrayList<String>();
		for(String x : cmdInput){
			args.add(x);
		}
		
		args.remove(0); // Remove command name
		
		boolean invertSelection= Utils.argListContainsFlag(args, "-V");
		
		List<String> trackNameRegex= new ArrayList<String>();
		String awk= "";

		if(args.size() == 0){
			// This will turn off everything
			trackNameRegex.add(".*");
		} else if(args.get(0).equals("-off")){
			args.remove(0);
			trackNameRegex.addAll(args); // Everything after -off is track re to turn off. 
		} else {
		
			// We need to find the awk script (i.e. '$3 > 10'), everything after it is track regexes.
			// * Iterate through the args, if arg starts with - skip this arg and the next one.
			// * The first arg (string) not starting with - is the script
			// * Any args after the script are track regex.
			final List<String> awkOpts= Arrays.asList(new String[] {"-F", "-f", "-v", "-t", "-c", "-o", "-z", "-Z", "-d", "-S", "-s", "-x", "-y", "-r", "-ext", "-ni"});
			int idxScript= 0; // Index of the script in the command args. 
			boolean skip= false;
			for(String x : args){
				if(awkOpts.contains(x)){
					idxScript += 2;
					skip= true;
				} else if(skip){
					skip= false;
					continue;
				} else {
					break;
				}
			}
			args.set(idxScript, "'" +  args.get(idxScript) + "'"); // Put back single quotes around the script, exclude cmd line params like -F 
			awk= Joiner.on(" ").join(args.subList(0, idxScript+1));

			// Everything after the script is track regexes
			trackNameRegex.addAll(args.subList(idxScript+1, args.size()));
		}

		if(trackNameRegex.size() == 0){
			trackNameRegex.add(".*"); // Track regex list not given: Set to capture all of them.
		}
		
		// Set script
        List<Track> tracksToReset = this.matchTracks(trackNameRegex, true, invertSelection);
        for(Track tr : tracksToReset){
        	if(awk.contains("getSamTag(") && ! tr.getTrackFormat().equals(TrackFormat.BAM)){
        		System.err.println("\nFunction getSamTag() can be applied to BAM tracks only. Got:\n" + tr.getTrackTag());
        		throw new InvalidCommandLineException();
        	}
        	if((awk.contains("getInfoTag(") || awk.contains("getFmtTag(")) && ! tr.getTrackFormat().equals(TrackFormat.VCF)){
        		System.err.println("\nFunction getInfoTag(), getFmtTag() can be applied to VCF tracks only. Got:\n" + tr.getTrackTag());
        		throw new InvalidCommandLineException();
        	}
        	if((awk.contains("getGtfTag(") || awk.contains("getGtfTag(")) && ! tr.getTrackFormat().equals(TrackFormat.GTF)){
        		System.err.println("\nFunction getGtfTag() can be applied to GTF tracks only. Got:\n" + tr.getTrackTag());
        		throw new InvalidCommandLineException();
        	}
        	if((awk.contains("getGffTag(") || awk.contains("getGffTag(")) && ! tr.getTrackFormat().equals(TrackFormat.GFF)){
        		System.err.println("\nFunction getGffTag() can be applied to GFF tracks only. Got:\n" + tr.getTrackTag());
        		throw new InvalidCommandLineException();
        	}
        	String script= this.replaceAwkHeaders(awk, tr.getTrackFormat());
        	script= this.replaceAwkFuncs(script, tr);
        	tr.setAwk(script);
        }
	}

	/**Replace in awk script the overloaded function name(s) with the actual names and args 
	 * @throws InvalidCommandLineException 
	 * */
	private String replaceAwkFuncs(String awkScript, Track track) throws InvalidCommandLineException{

		final String FUNC= "get"; // Function name to be replaced 
		
		// Search for 'foo FUNC(...) bar'. Use parenthesis to get the group "FUNC(...)".  
		Pattern pattern= Pattern.compile(".*(" + FUNC + "\\(.+?\\)).*");
		
		while(awkScript.matches(".*\\b" + FUNC + "\\(.*")){
			Matcher m = pattern.matcher(awkScript);
			m.matches();
			String func= m.group(1); // this looks like "get(DP, 1, 2)"
			int gapStart= m.start(1);
			int gapEnd= m.start(1) + func.length();

			// Remove func name and brackets, leave only param list. E.g. ["DP", "1", "2"]
			String xargs= func.replaceAll(".*\\(", "").replaceAll("\\).*", "");
			List<String> args= new ArrayList<String>();
			args.addAll(Splitter.on(",").splitToList(xargs));
			
			// Get the first arg i.e., the tag name, possibly in double quotes. 
			String tag= args.get(0).
					replaceAll(".*\\(", "").
					replaceAll("\\).*", "").
					trim().
					replaceAll("^\"", "").
					replaceAll("\"$", "");
			
			// Depending on track type and tag we decide what the replacement function is
			String fname;
			if(track.getTrackFormat().equals(TrackFormat.BAM)){
				fname= "getSamTag";
			} 
			else if(track.getTrackFormat().equals(TrackFormat.VCF)){
				if(tag.trim().startsWith("FMT/") || track.getVcfHeader().getFormatHeaderLine(tag) != null){
					fname= "getFmtTag";
				}
				else if(tag.trim().startsWith("INFO/") || track.getVcfHeader().getInfoHeaderLine(tag) != null){
					fname= "getInfoTag";
				}
				else {
					System.err.println("Tag " + tag + " not found in VCF header of track " + track.getTrackTag() + ".\n"
							+ "Please prepend 'INFO/' or 'FMT/' to the tag to search the INFO or FORMAT fields, respectively.");
					throw new InvalidCommandLineException();
				}
			}
			else if(track.getTrackFormat().equals(TrackFormat.GTF)){
				fname= "getGtfTag";
			}
			else if(track.getTrackFormat().equals(TrackFormat.GFF)){
				fname= "getGffTag";
			}
			else {
				System.err.println("Function " + FUNC + " is not available for track type " + track.getTrackFormat() + ".");
				throw new InvalidCommandLineException();
			}
			args.set(0, '"' + tag + '"');
			func= fname + "(" + Joiner.on(",").join(args) + ")";
			awkScript= awkScript.substring(0, gapStart) + func + awkScript.substring(gapEnd, awkScript.length());
		}
		return awkScript;
	}
	
	private String replaceAwkHeaders(final String awkScript, TrackFormat trackFormat) throws InvalidCommandLineException{
		
		Map<TrackFormat, Map<String, String>> headers= new HashMap<TrackFormat, Map<String, String>>();
		
		Map<String, String> bamHdr= new HashMap<String, String>();
		bamHdr.put("$QNAME", "$1");
		bamHdr.put("$FLAG", "$2");
		bamHdr.put("$RNAME", "$3"); bamHdr.put("$CHROM", "$3"); bamHdr.put("$SEQNAME", "$3");
		bamHdr.put("$POS", "$4"); bamHdr.put("$START", "$4");
		bamHdr.put("$MAPQ", "$5");
		bamHdr.put("$CIGAR", "$6");
		bamHdr.put("$RNEXT", "$7");
		bamHdr.put("$PNEXT", "$8");
		bamHdr.put("$TLEN", "$9");
		bamHdr.put("$SEQ", "$10");
		bamHdr.put("$QUAL", "$11");
		headers.put(TrackFormat.BAM, bamHdr);
		
		Map<String, String> vcfHdr= new HashMap<String, String>();
		vcfHdr.put("$CHROM", "$1"); vcfHdr.put("$RNAME", "$1");
		vcfHdr.put("$POS", "$2"); vcfHdr.put("$START", "$2");
		vcfHdr.put("$ID", "$3");
		vcfHdr.put("$REF", "$4");
		vcfHdr.put("$ALT", "$5");
		vcfHdr.put("$QUAL", "$6");
		vcfHdr.put("$FILTER", "$7");
		vcfHdr.put("$INFO", "$8");
		vcfHdr.put("$FORMAT", "$9");
		headers.put(TrackFormat.VCF, vcfHdr);
		
		Map<String, String> gffHdr= new HashMap<String, String>();
		gffHdr.put("$SEQNAME", "$1"); gffHdr.put("$CHROM", "$1"); gffHdr.put("$RNAME", "$1");
		gffHdr.put("$SOURCE", "$2");
		gffHdr.put("$FEATURE", "$3");
		gffHdr.put("$START", "$4"); gffHdr.put("$POS", "$4");
		gffHdr.put("$END", "$5");
		gffHdr.put("$SCORE", "$6");
		gffHdr.put("$STRAND", "$7");
		gffHdr.put("$FRAME", "$8");
		gffHdr.put("$ATTRIBUTE", "$9");
		headers.put(TrackFormat.GFF, gffHdr);
		headers.put(TrackFormat.GTF, gffHdr);
		
		Map<String, String> bedHdr= new HashMap<String, String>();
		bedHdr.put("$CHROM", "$1"); bedHdr.put("$RNAME", "$1"); bedHdr.put("$SEQNAME", "$1");
		bedHdr.put("$START", "$2"); bedHdr.put("$POS", "$2");
		bedHdr.put("$END", "$3");
		bedHdr.put("$NAME", "$4");
		bedHdr.put("$SCORE", "$5");
		bedHdr.put("$STRAND", "$6");
		bedHdr.put("$THICKSTART", "$7");
		bedHdr.put("$THICKEND", "$8");
		bedHdr.put("$RGB", "$9");
		bedHdr.put("$BLOCKCOUNT", "$10");
		bedHdr.put("$BLOCKSIZES", "$11");
		bedHdr.put("$BLOCKSTARTS", "$12");
		headers.put(TrackFormat.BED, bedHdr);
		
    	Map<String, String> hdrs= headers.get(trackFormat);
    	if(hdrs == null){
    		// This happens for track types which don't have headers. e.g. TrackWiggle
    		hdrs= new HashMap<String, String>();
    		return awkScript;
    	}

    	String awk2= awkScript;
    	for(String hdr : hdrs.keySet()){
    		awk2= awk2.replace(hdr, hdrs.get(hdr));
    	};
    	// Check there are no headers associated to other track formats. 
    	for(TrackFormat fmt : headers.keySet()){
        	Map<String, String> hdrs2 = headers.get(fmt);
        	for(String hdr : hdrs2.keySet()){
        		if(awk2.contains(hdr)){
        			System.err.println("Column header '" + hdr + "' cannot be applied to track format '" + trackFormat + "'");
        			throw new InvalidCommandLineException();
        		}
        	};
    	}
    	return awk2;
	}
	
	public void setFilterVariantReads(List<String> cmdInput) throws InvalidCommandLineException, MalformedURLException, ClassNotFoundException, IOException, InvalidGenomicCoordsException, InvalidRecordException, SQLException {
		
		List<String> args= new ArrayList<String>();
		for(String x : cmdInput){
			args.add(x);
		}
		
		args.remove(0); // Remove command name
		
		// Collect arguments
		boolean variantOnly= ! Utils.argListContainsFlag(args, "-all");
		boolean invertSelection= Utils.argListContainsFlag(args, "-v");

		String region= Utils.getArgForParam(args, "-r", null);
		String chrom= Filter.DEFAULT_VARIANT_CHROM.getValue();
		int from= -1;
		int to= -1;
		if(region != null && this.getTrackList().size() > 0){
		        
		        if(this.getTrackList().get(0).getGc().getFastaFile() == null){
		                System.err.println("A reference sequence is required to execute this command.");
		                throw new InvalidCommandLineException();
		        }
		        
		        chrom= this.getTrackList().get(0).getGc().getChrom();
		        // Find whether the region is a single pos +/- a value or  
		        region= region.trim().replaceAll(" ", "").replaceAll(",", "").replaceAll("/", "");
		        try{
		                if(region.contains("+-") || region.contains("-+")){
		                        String[] fromTo= region.replaceAll("\\+", "").split("-"); 
		                        int offset= Integer.parseInt(fromTo[1]);
		                        from= Integer.parseInt(fromTo[0]) - offset;
		                        to= Integer.parseInt(fromTo[0]) + offset;
		                }               
		                else if(region.contains("-")){
		                        String[] fromTo= region.split("-");
		                        from= Integer.parseInt(fromTo[0]) - Integer.parseInt(fromTo[1]);
		                        to= Integer.parseInt(fromTo[0]);
		                }               
		                else if(region.contains("+")){
		                        String[] fromTo= region.split("\\+");
		                        from= Integer.parseInt(fromTo[0]);
		                        to= Integer.parseInt(fromTo[0]) + Integer.parseInt(fromTo[1]);
		                }               
		                else if(region.contains(":")){
		                        String[] fromTo= region.split(":");
		                        from= Integer.parseInt(fromTo[0]);
		                        to= Integer.parseInt(fromTo[1]);
		                }               
		                else {          
		                        from= Integer.parseInt(region);
		                        to= Integer.parseInt(region);
		                }               
		        } catch(NumberFormatException e) {
		                System.err.println("Cannot parse region into integers: " + region);
		                throw new InvalidCommandLineException();
		        }
		        if(from < 1) from= 1;
		        if(to < 1) to= 1;
		        if(from > to) from= to; 
		}
		// -----------------------------------------------
		// Stub for filtering by records in another track 
//		String track= Utils.getArgForParam(args, "-t", null);
//		if(track != null){
//			List<String>pattern= new ArrayList<String>();
//			pattern.add(track);
//			List<Track> tname= this.matchTracks(pattern, true, false);
//			if(tname.size() == 0){
//				// Handle no track matched
//			}
//			TrackIntervalFeature tr= (TrackIntervalFeature)this.getTrack(tname.get(0));
//			tr.getIntervalFeatureList();
//			// List<String> chroms=
//			// List<Integer> starts=
//			// List<ends> ends=
//		}
		// Modify Track.setVariantReadInInterval() to accept a list of chrom, start, end positions.
		// -----------------------------------------------
		
		List<String> trackNameRegex= new ArrayList<String>();
		
		// Regex for tracks:
		if(args.size() == 0){
			// This will turn off everything
			trackNameRegex.add(".*");
		} else {
			trackNameRegex.addAll(args); // Everything after -off is track re to turn off. 
		}
	
		// Set filter
        List<Track> tracksToReset = this.matchTracks(trackNameRegex, true, invertSelection);
        for(Track tr : tracksToReset){
        	tr.setVariantReadInInterval(chrom, from, to, variantOnly);
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

		Pattern showRegex= Pattern.compile(Filter.DEFAULT_SHOW_REGEX.getValue()); // Default
		Pattern hideRegex= Pattern.compile(Filter.DEFAULT_HIDE_REGEX.getValue());
		
		// Get args:
		boolean invertSelection= Utils.argListContainsFlag(args, "-v");

		int flag= 0;
		if(Utils.argListContainsFlag(args, "-F")){
			flag |= Pattern.LITERAL; 
		}
		if( ! Utils.argListContainsFlag(args, "-c")){
			flag |= Pattern.CASE_INSENSITIVE;  
		}
		
		try{
			if(args.contains("-i")){
				int idx= args.indexOf("-i") + 1;
				showRegex= Pattern.compile(args.get(idx), flag);
				args.remove(idx);
				args.remove("-i");
			}
			if(args.contains("-e")){
				int idx= args.indexOf("-e") + 1; 
				hideRegex= Pattern.compile(args.get(idx), flag);
				args.remove(idx);
				args.remove("-e");
			}
		} catch(PatternSyntaxException e){
	    	System.err.println("Invalid regex");
	    	throw new InvalidCommandLineException();
		}
		// What is left is positional args of regexes
		List<String> trackNameRegex= new ArrayList<String>();
		trackNameRegex.addAll(args);
		if(trackNameRegex.size() == 0){
			trackNameRegex.add(".*"); // Default regex for matching tracks
		}		
		
		// TRACK REGEXES
        // Regex
        // And set as required:
        List<Track> tracksToReset = this.matchTracks(trackNameRegex, true, invertSelection);
        for(Track tr : tracksToReset){
        	tr.setShowHideRegex(showRegex, hideRegex);
        }
	}
	
	/** Go to the next feature on trackId given the current GenomicCoordinates. 
	 * 
	 * If slop is > 0, the output coordinates are centered on the feature and extended
	 * slop times the size of the feature left and right. With slop= 0 the coordinates 
	 * are exactly spanning the feature. With slop < 0 the output coordinates have 
	 * the feature right at the start.
	 * */
	public GenomicCoords goToNextFeatureOnFile(String trackId, GenomicCoords currentGc, double slop, boolean getPrevious) throws InvalidGenomicCoordsException, IOException, InvalidCommandLineException{

		Track tr= this.matchIntervalFeatureTrack(trackId.trim());

		if(tr == null){
			return currentGc;
		}
		
		TrackIntervalFeature tif= (TrackIntervalFeature) tr;
		if(slop < 0){
			return tif.coordsOfNextFeature(currentGc, getPrevious);
		} else {
			GenomicCoords featureGc= tif.startEndOfNextFeature(currentGc, getPrevious);
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
			List<Track> matched= matchTracks(x, true, false);
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

	private List<Track> matchTracks(List<String> patterns, boolean asRegex, boolean invertSelection) throws InvalidCommandLineException{

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
		
		LinkedHashSet<Track> matchedTracks= new LinkedHashSet<Track>(); 
		
		for(Track track : this.getTrackList()){
			String trackId= track.getTrackTag();
			for(String pattern : patterns){
				boolean matched= Pattern.compile(pattern).matcher(trackId).find();
				if(matched){
					matchedTracks.add(track);
				}
			}
		}
		if(invertSelection){
			List<Track>inv= new ArrayList<Track>();
			for(Track x : this.getTrackList()){
				if( ! matchedTracks.contains(x)){
					inv.add(x);
				}
			}
			return inv;
		} else {
			return new ArrayList<Track>(matchedTracks);
		}
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

	public GenomicCoords findNextMatchOnTrack(Pattern pattern, String trackId, GenomicCoords currentGc, boolean all) throws InvalidGenomicCoordsException, IOException, InvalidCommandLineException{

		TrackIntervalFeature tif= (TrackIntervalFeature) matchIntervalFeatureTrack(trackId.trim());
		if(tif == null){
			return currentGc;
		}

		System.err.println("Matching on " + tif.getTrackTag());
		
		if(all){
			return tif.genomicCoordsAllChromMatchInGenome(pattern, currentGc);
		} else {
			return tif.findNextMatch(currentGc, pattern);
		}
	}

	private List<TrackIntervalFeature> getIntervalFeatureTracks(){
		
		// TrackSet ifSet= new TrackSet();
		List<TrackIntervalFeature> ifSet= new ArrayList<TrackIntervalFeature>();
		
		for(Track tr : this.getTrackList()){
			if(tr instanceof TrackIntervalFeature){
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
        boolean invertSelection= Utils.argListContainsFlag(tokens, "-v");

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
        List<Track> tracksToReset = this.matchTracks(trackNameRegex, true, invertSelection);
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
	public void orderTracks(List<String> tokens) throws InvalidCommandLineException {
	
		// List<String> args= new ArrayList<String>(tokens);
		// boolean invertSelection= Utils.argListContainsFlag(args, "-v");
		
		List<String> newOrder= new ArrayList<String>(tokens);
		
		if(newOrder.size() == 0){
			this.sortTracksByTagName();
			return;
		}
		
		// Create a new list that will have the new order
		List<Track> newTrackList= new ArrayList<Track>();
		
		for(String query : newOrder){
			List<String> x= new ArrayList<String>();
			x.add(query);
			List<Track> trList = this.matchTracks(x, true, false);
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
			filenames.add(tr.getFilename());
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
	
	/** Handle bookmarks by processing cmd line args.
	 * @throws InvalidGenomicCoordsException 
	 * @throws SQLException 
	 * @throws InvalidRecordException 
	 * @throws IOException 
	 * @throws ClassNotFoundException 
	 * @throws InvalidCommandLineException 
	 * */
	public String bookmark(GenomicCoords gc, List<String> cmdInput) throws ClassNotFoundException, IOException, InvalidRecordException, SQLException, InvalidGenomicCoordsException, InvalidCommandLineException{
		
		String messages= "";
		
		List<String> args= new ArrayList<String>(cmdInput);
		args.remove(0); // Remove command name

		// Get all arguments. What is left is the positional argument 
		boolean delete= Utils.argListContainsFlag(args, "-d");
		String name= Utils.getArgForParam(args, "-n", null);
		boolean print= Utils.argListContainsFlag(args, "-print");
		String file= Utils.getArgForParam(args, ">", null);

		GenomicCoords bookmarkRegion= null;
		if(args.size() > 0){
			List<String> strRegion= Utils.parseStringCoordsToList(args.get(0));
			if(strRegion.get(1) == null){
				strRegion.set(1, "1");
			}
			if(strRegion.get(2) == null){
				strRegion.set(2, new Integer(Integer.MAX_VALUE).toString());
			}
			bookmarkRegion= new GenomicCoords(strRegion.get(0) + ":" + strRegion.get(1) + "-" + strRegion.get(2), 
					gc.getUserWindowSize(), gc.getSamSeqDict(), gc.getFastaFile());
		} else {
			bookmarkRegion= new GenomicCoords(gc.toStringRegion(), gc.getUserWindowSize(), gc.getSamSeqDict(), gc.getFastaFile());
		}

		if(print){
			for(Track tr : this.getTrackList()){
				if(tr instanceof TrackBookmark){
					List<String> marks= Utils.tabulateList(((TrackBookmark)tr).asList(), -1);
					messages= Joiner.on("\n").join(marks);
					return messages + "\n";
				}
			}
			return messages;
		}

		if(file != null){
			for(Track tr : this.getTrackList()){
				if(tr instanceof TrackBookmark){
					((TrackBookmark)tr).save(file, false);
	 				return messages;
				}
			}
			return messages;
		}
		
		if(delete){
			for(Track tr : this.getTrackList()){
				if(tr instanceof TrackBookmark){
					((TrackBookmark)tr).removeBookmark(bookmarkRegion);
	 				return messages;
				}
			}
			return messages;
		}

		if(name == null){
			name= ".";
		}
		this.addBookmark(bookmarkRegion, name);
		return messages;
	}
	
	/** Add position gc to bookmark track. Track created if it does not exist. 
	 * */
	private void addBookmark(GenomicCoords gc, String nameForBookmark) throws ClassNotFoundException, IOException, InvalidRecordException, SQLException, InvalidGenomicCoordsException{

		// Check there isn't a bookmark track already:
		for(Track tr : this.getTrackList()){
			if(tr instanceof TrackBookmark){ // A Bookmark track exists, add position to it
				tr.addBookmark(gc, nameForBookmark);
 				tr.setGc(gc);
 				return;
			}
		}
		// We need to create the bookmark track.
		TrackBookmark tr= new TrackBookmark(gc, nameForBookmark);
		this.addTrack(tr, "Bookmarks");
		tr.setGc(gc);
		
	}
	
	/*   S e t t e r s   and   G e t t e r s  */

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
		int f= Integer.valueOf(Filter.DEFAULT_f_FLAG.getValue());
		int F= Integer.valueOf(Filter.DEFAULT_F_FLAG.getValue());
		int q= Integer.valueOf(Filter.DEFAULT_MAPQ.getValue());
		
		// Get args:
		boolean invertSelection= Utils.argListContainsFlag(args, "-v");
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
        List<Track> tracksToReset = this.matchTracks(trackNameRegex, true, invertSelection);
        
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
	
	public String editNamesForRegex(List<String> cmdInput) throws InvalidCommandLineException {
		// API:
		// editNames <patterns> <replacement> [track_re] ...
		List<String> args= new ArrayList<String>(cmdInput);
		args.remove(0); // Remove name of command

		boolean invertSelection= Utils.argListContainsFlag(args, "-v");

		boolean test= false;
		if(args.contains("-t")){
			test= true;
			args.remove("-t");
		}
		
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
        List<Track> tracksToReset = this.matchTracks(args, true, invertSelection);
        
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
        	System.err.println("Edit not allowed as it generates duplicate tags.");
        	throw new InvalidCommandLineException();
        }
        
        // Now change for real
        String messages= "";
        if(tracksToReset.size() == 0){
        	messages += "No track matched";
        }
        for(Track tr : tracksToReset){
        	String newTag= tr.getTrackTag().replaceAll(pattern, replacement);
        	messages += "Renaming " + tr.getTrackTag() + " to " + newTag + "\n";
        	if( ! test ){
        		tr.setTrackTag(newTag);
        	}
        }
        return messages;
	}

	public void setReadsAsPairsForRegex(List<String> tokens) throws InvalidCommandLineException, InvalidGenomicCoordsException, IOException {

		List<String> args= new ArrayList<String>(tokens);
		args.remove(0);

		boolean invertSelection= Utils.argListContainsFlag(args, "-v");

		Boolean readsAsPairs= null;
		if(args.contains("-on")){
			readsAsPairs= true;
			args.remove("-on");
		}
		if(args.contains("-off")){
			readsAsPairs= false;
			args.remove("-off");
		}
		
        // Regex
        List<String> trackNameRegex= new ArrayList<String>();
        if(args.size() > 0){
            trackNameRegex= args;
        } else {
            trackNameRegex.add(".*"); // Default: Capture everything
        }

        // And set as required:
        List<Track> tracksToReset = this.matchTracks(trackNameRegex, true, invertSelection);
        for(Track tr : tracksToReset){
        	if(readsAsPairs == null){
    			if(tr.getReadsAsPairs()){ // Invert setting
    				tr.setReadsAsPairs(false);
    			} else {
    				tr.setReadsAsPairs(true);
    			}
        	} else if(readsAsPairs) {
        		tr.setReadsAsPairs(true);
        	} else {
        		tr.setReadsAsPairs(false);
        	}
        }
	}
	
	public void setFeatureGapForRegex(List<String> tokens) throws InvalidCommandLineException {

		List<String> args= new ArrayList<String>(tokens);
		args.remove(0);

		boolean invertSelection= Utils.argListContainsFlag(args, "-v");

		Boolean gap= null;
		if(args.contains("-on")){
			gap= true;
			args.remove("-on");
		}
		if(args.contains("-off")){
			gap= false;
			args.remove("-off");
		}
		
        // Regex
        List<String> trackNameRegex= new ArrayList<String>();
        if(args.size() > 0){
            trackNameRegex= args;
        } else {
            trackNameRegex.add(".*"); // Default: Capture everything
        }

        // And set as required:
        List<Track> tracksToReset = this.matchTracks(trackNameRegex, true, invertSelection);
        for(Track tr : tracksToReset){
        	if(gap == null){
    			if(tr.getGap() == 0){ // Invert setting
    				tr.setGap(1);
    			} else {
    				tr.setGap(0);
    			}
        	} else if(gap) {
        		tr.setGap(1);
        	} else {
        		tr.setGap(0);
        	}
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
	 * @param genomicCoords 
	 * @throws SQLException 
	 * @throws InvalidRecordException 
	 * @throws InvalidGenomicCoordsException 
	 * @throws IOException 
	 * @throws ClassNotFoundException 
	 * @throws InvalidCommandLineException 
	 * */
	public void setRegexForTrackSeqRegex(List<String> cmdInput, GenomicCoords genomicCoords) throws ClassNotFoundException, IOException, InvalidGenomicCoordsException, InvalidRecordException, SQLException, InvalidCommandLineException {

		initTrackSeqRegex(genomicCoords);
		
		List<String> args= new ArrayList<String>(cmdInput);
		args.remove(0); // Remove command name
		
		boolean isCaseSensisitive= false;
		if(args.contains("-c")){
			isCaseSensisitive= true;
			args.remove("-c");
		}

		boolean isIupac= false;
		if(args.contains("-iupac")){
			isIupac= true;
			args.remove("-iupac");
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
				((TrackSeqRegex) tr).setIupac(isIupac);
				tr.setSeqRegex(seqRegex);
			}
		}
	}

	/** Initialize TrackSeqRegex for regex matches in fasta. Do nothing if track already exists. 
	 * Throw exception if cannot set it.
	 * This method only initializes the track. It doesn't set the regex or search for it.*/
	private void initTrackSeqRegex(GenomicCoords genomicCoords) throws InvalidCommandLineException, ClassNotFoundException, IOException, InvalidGenomicCoordsException, InvalidRecordException, SQLException {
		// See if a TrackSeqRegex exists.
		boolean found= false;
		for(Track tr : this.getTrackList()){
			if(tr instanceof TrackSeqRegex){
				found= true;
			}
		}
		if( ! found){
			// Not found: Try to create a new TrackSeqRegex
			if(genomicCoords.getFastaFile() != null){
				TrackSeqRegex re= new TrackSeqRegex(genomicCoords);
				this.addTrack(re, "seqRegex");
				found= true;
			}
		}
		if(! found){
			throw new InvalidCommandLineException();
		}
	}

	public String dropTracksWithRegex(List<String> cmdInput) throws InvalidCommandLineException {

		List<String> args= new ArrayList<String>(cmdInput);
		args.remove(0); //Remove cmd name
		boolean invertSelection= Utils.argListContainsFlag(args, "-v");

		boolean test= false;
		if(args.contains("-t")){
			test= true;
			args.remove("-t");
		}
		List<String> trackNameRegex= new ArrayList<String>(args);
       
		String messages= "";
        List<Track> tracksToDrop= this.matchTracks(trackNameRegex, true, invertSelection);
        if(tracksToDrop.size() == 0){
        	messages += "No track matched by regex: " + trackNameRegex;
        }
        for(Track tr : tracksToDrop){
        	messages += "Dropping: " + tr.getTrackTag() + "\n";
        	if( ! test){
        		this.trackList.remove(tr);
        	}
        }
        return messages;
        
	}

	public GenomicCoords trimCoordsForTrack(List<String> cmdInput) throws InvalidGenomicCoordsException, IOException {

		List<String> args= new ArrayList<String>(cmdInput);
		args.remove(0); //Remove cmd name

		String trackTag= "";
		if(args.size() > 0){
			trackTag= args.get(0);
		}
		
		// Get track to trim: We get the first one matching the given tag.
		for(Track tr : this.getTrackList()){
			if(tr instanceof TrackIntervalFeature && tr.getTrackTag().contains(trackTag)){
				return trimTrack((TrackIntervalFeature) tr);
			}
		}
		return null;
	}

	private GenomicCoords trimTrack(TrackIntervalFeature tr) throws InvalidGenomicCoordsException, IOException{
		GenomicCoords current = tr.getGc();
		TrackIntervalFeature itr= (TrackIntervalFeature)tr;
		List<IntervalFeature> features = itr.getIntervalFeatureList();
		if(features.size() == 0){
			return current;
		}
		
		// Get the leftmost and rightmost coordinates of the feature set
		// These might extend beyond the current window.
		int left = Integer.MAX_VALUE;
		int right= 0;
		for(IntervalFeature f : features){
			if(f.getFrom() < left){
				left= f.getFrom(); 
			}
			if(f.getTo() > right){
				right= f.getTo();
			}
		}
		
		// See if left and right need to be adjusted by window coordinates
		if(left < current.getFrom()){
			left= current.getFrom();
		}
		if(right > current.getTo()){
			right= current.getTo();
		}
		return new GenomicCoords(current.getChrom() + ":" + left + "-" + right, 
				current.getUserWindowSize(), current.getSamSeqDict(), current.getFastaFile());
	}

	public List<Track> getTracksForYLimits() {
		return tracksForYLimits;
	}

	public void setTracksForYLimits(List<Track> tracksForYLimits) {
		this.tracksForYLimits = tracksForYLimits;
	}

	
//	public List<String> getRegexForYLimits() {
//		return regexForYLimits;
//	}
//
//	public void setRegexForYLimits(List<String> regexForYLimits) {
//		this.regexForYLimits = regexForYLimits;
//	}

	public String[] getYStringLimits() {
		return yStrLimits;
	}

	public void setYStrLimits(String[] yStrLimits) {
		this.yStrLimits = yStrLimits;
	}

	public LinkedHashSet<String> getOpenedFiles() {
		return openedFiles;
	}

	public void setOpenedFiles(LinkedHashSet<String> openedFiles) {
		this.openedFiles = openedFiles;
	}

	public String showRecentlyOpened(List<String> cmdInput) throws InvalidCommandLineException {
		
		List<String> args= new ArrayList<String>(cmdInput);
		args.remove(0);
		
		String re= ".*";
		int nmax= Integer.MAX_VALUE;
		if(args.contains("-grep")){
			re= args.get(args.indexOf("-grep") + 1);
		}
		if(args.contains("-n")){
			nmax= Integer.parseInt(args.get(args.indexOf("-n") + 1));
		}
		
		Pattern pattern= Pattern.compile(re); // .matcher(x).find();
		List<String> opened= new ArrayList<String>();
		for(String x : this.getOpenedFiles()){
			if(pattern.matcher(x).find()){
				opened.add(x);
			}
		}
		if(opened.size() > nmax){ // Trim list to 
			opened= opened.subList(opened.size() - nmax, opened.size());
		}
		// Add index. 
		// We add the index of the file from the full list of opened files, not 
		// from the list returned by recentlyOpened.
		List<String> openedFilesRev= Lists.reverse(new ArrayList<String>(this.getOpenedFiles()));
		List<String> toshow= Lists.reverse(opened);
		for(int i= 0; i < opened.size(); i++){
			int idx= openedFilesRev.indexOf(toshow.get(i)) + 1;
			toshow.set(i, idx + "\t" + toshow.get(i));
		}
		toshow= Lists.reverse(toshow);
		return Joiner.on("\n").join(Utils.tabulateList(toshow, -1));
	}

	/** Merge the set of opened files with the given (historic) list. 
	 * */ 
	public void addHistoryFiles(List<String> historyFiles) {
		LinkedHashSet<String> union= new LinkedHashSet<String>();
		for(String x : historyFiles){
			union.add(x);
		}
		LinkedHashSet<String> now= this.getOpenedFiles();
		for(String file : now){ // If a file is in the current track set and in the history file, put it last. I.e. last opened. 
			if(union.contains(file)){
				union.remove(file);
			}
			union.add(file);
		}
		this.setOpenedFiles(union);
	}

	public void reload(List<String> cmdTokens) throws InvalidCommandLineException, ClassNotFoundException, InvalidGenomicCoordsException, IOException, InvalidRecordException, SQLException {
		List<String> args= new ArrayList<String>(cmdTokens);
		args.remove(0); //Remove cmd name
		boolean invertSelection= Utils.argListContainsFlag(args, "-v");

		if(args.size() == 0){
			args.add(".*");
		}
		List<String> trackNameRegex= new ArrayList<String>(args);
       
        List<Track> tracksToReload= this.matchTracks(trackNameRegex, true, invertSelection);
        for(Track tr : tracksToReload){
        	try{
        		tr.reload();
        	} catch(Exception e) {
        		System.err.println(e.toString());
            	this.trackList.remove(tr);
        	}
        }
	}
}
