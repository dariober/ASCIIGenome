package samTextViewer;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Pattern;
import java.util.regex.PatternSyntaxException;

import com.google.common.base.Joiner;
import com.google.common.base.Splitter;

import coloring.Png;
import commandHelp.CommandList;
import exceptions.BamIndexNotFoundException;
import exceptions.InvalidCommandLineException;
import exceptions.InvalidGenomicCoordsException;
import exceptions.InvalidRecordException;
import net.sourceforge.argparse4j.inf.Namespace;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import jline.console.ConsoleReader;
import tracks.Track;
import tracks.TrackIntervalFeature;
import tracks.TrackSet;
import tracks.TrackWiggles;

/**
 * @author berald01
 *
 */
public class Main {
	
	private static String getMemoryStat(){
		float mem= (float) ((float)Runtime.getRuntime().totalMemory() / 1000000d);
		String memStats= "Mem: " +  Math.round(mem * 10)/10 + " MB";
		return memStats;
	}

	
	public static void main(String[] args) throws IOException, InvalidGenomicCoordsException, InvalidCommandLineException, InvalidRecordException, BamIndexNotFoundException, ClassNotFoundException, SQLException {
		
		Runtime.getRuntime().addShutdownHook(new Thread(new Runnable() {
		    public void run() {
		    	System.out.print("\033[0m");
		    	try {
					ConsoleReader x= new ConsoleReader();
					x.clearScreen();
					x.flush();
		    	} catch (IOException e) {
					e.printStackTrace();
				}
		    }
		}));
		
		/* Start parsing arguments * 
		 * *** If you change something here change also in console input ***/
		Namespace opts= ArgParse.argParse(args);
		
		List<String> initFileList= opts.getList("input");
		String region= opts.getString("region");
		final String genome= opts.getString("genome");
		final String fasta= opts.getString("fasta");
		final String exec= opts.getString("exec");
		final boolean noFormat= opts.getBoolean("noFormat");
		final boolean nonInteractive= opts.getBoolean("nonInteractive");
		
		/* Set up console */
		
		int windowSize= 160;
		try{
			windowSize= jline.TerminalFactory.get().getWidth() - 1;			
		} catch(Exception e){
			e.printStackTrace();
		}
		
		Utils.checkFasta(fasta);
		
		/* Test input files exist */
		List<String> inputFileList= new ArrayList<String>();
		try{
			Utils.addSourceName(inputFileList, initFileList);
		} catch(InvalidCommandLineException e){
			//
		}	
		
		if((region == null || region.isEmpty()) && fasta != null){ // Try to initilize from fasta
			IndexedFastaSequenceFile faSeqFile = new IndexedFastaSequenceFile(new File(fasta));
			region= faSeqFile.nextSequence().getName();
			faSeqFile.close();
		}

		GenomicCoordsHistory gch= new GenomicCoordsHistory();
		SAMSequenceDictionary samSeqDict = GenomicCoords.getSamSeqDictFromAnyFile(inputFileList, fasta, genome);
		
		/* Prepare genomic coordinates to fetch. This should probably be a function in itself */
		if(region.isEmpty()){
			System.err.print("Initializing coordinates... ");
			if(!samSeqDict.isEmpty()){
				region= samSeqDict.getSequence(0).getSequenceName();
			} else {
				for(String x : inputFileList){
					try {
						region= Utils.initRegionFromFile(x);
						System.err.println("Done from: " + x);
						break;
					} catch(Exception e){
						System.err.println("\nCould not initilize from file " + x);
					}
				}
			}
		}
		gch.add(new GenomicCoords(region, samSeqDict, windowSize, fasta));

		/* Initialize trackSet */
		TrackSet trackSet= new TrackSet(inputFileList, gch.current());
		
		/* Initialize GC profile */
		if(fasta != null){
			TrackWiggles cgWiggle= gch.current().getGCProfile();
			trackSet.add(cgWiggle, "gcProfile");
		}
		
		// Put here the previous command so that it is re-issued if no imput is given
		// You have to initialize this var outside the while loop that processes input files.
		String currentCmdConcatInput= ""; 
		// String seqRegex= null;
		String snapshotFile= null;
		// boolean snapshotStripAnsi= true;
		boolean execDone= exec.isEmpty() ? true : false; // Do we need to execute commands from --exec? If exec is empty, consider it done.

		if(!noFormat){
			System.out.print("\033[48;5;231m");
		}		
		
		/* =================================== *
		 * Start processing interactive input  *
		 * =================================== */
		ConsoleReader console = CommandList.initConsole();
		
		while(true){ 


			if(!nonInteractive){ // == interactive
				console.clearScreen();
				console.flush();
			}
			
			if(gch.current().getChromIdeogram(20, noFormat) != null && execDone){
				Utils.printer(gch.current().getChromIdeogram(20, noFormat) + "\n", snapshotFile);
			}			
			
			for(Track track : trackSet.getTrackList()){
			
				track.setNoFormat(noFormat);

				track.setGc(gch.current());
				if(track.getyMaxLines() > 0){
					track.update();
				}
			
				if(execDone && track.getyMaxLines() > 0){
					Utils.printer(track.getTitle(), snapshotFile);
					Utils.printer(track.printToScreen() + "\n", snapshotFile);
					Utils.printer(track.getPrintableConsensusSequence(), snapshotFile);
					Utils.printer(track.printFeatures(windowSize), snapshotFile);
				}

			}

			/* Footers and interactive prompt */
			/* ****************************** */
			TrackIntervalFeature seqRegexTrack = gch.findRegex();
			// Track for regex match on fasta sequence See if we need to change height
			for(Pattern p : trackSet.getRegexForTrackHeight()){
				if(p.matcher(seqRegexTrack.getTrackTag()).find()){
					seqRegexTrack.setyMaxLines(trackSet.getTrackHeightForRegex());
					break;
				}
			}
			String seqPattern= seqRegexTrack.printToScreen();
			if(!seqPattern.isEmpty()){
				seqPattern += "\n";
			} 
			if(execDone){
				Utils.printer(seqPattern, snapshotFile);
				
				// Ruler and sequence
				Utils.printer(gch.current().printableRefSeq(noFormat), snapshotFile);
				String ruler= gch.current().printableRuler(10, noFormat);
				Utils.printer(ruler.substring(0, ruler.length() <= windowSize ? ruler.length() : windowSize) + "\n", snapshotFile);

				// Position, memory, etc
				String footer= gch.current().toString() + "; " + Math.rint(gch.current().getBpPerScreenColumn() * 10d)/10d + " bp/char; " + getMemoryStat();
				if(!noFormat){
					Utils.printer("\033[48;5;231;34m" + footer + "\033[48;5;231;30m; \n", snapshotFile);
				} else {
					Utils.printer(footer + "\n", snapshotFile);
				}

				// Optionally convert to png
				if(snapshotFile != null && snapshotFile.endsWith("png")){
					(new Png(new File(snapshotFile))).convert(new File(snapshotFile));
				} 
			}
			
			/* Interactive input */
			/* ================= */
			if(nonInteractive && execDone){
				// Exit now if non-interactive mode is set and no string is passed to be executed. 
				break;
			}
			
			/* =================================== *
			 * I N T E R A C T I V E    I N P U T  *
			 * =================================== */
			
			String cmdConcatInput= null; // String like "zi && -F 16 && mapq 10"
			snapshotFile= null; 

			while(cmdConcatInput == null){
				
				if(execDone){
					console.setPrompt("[h] for help: ");
					cmdConcatInput= console.readLine().trim();
					
					if (cmdConcatInput.isEmpty() || cmdConcatInput == null){
						// Repeat previous command
						cmdConcatInput= currentCmdConcatInput;
					}
				} else {
					cmdConcatInput= exec;
					execDone= true;
				}
				// cmdInputList: List of individual commands in tokens to be issued. 
				// E.g.: [ ["zi"], 
				//         ["-F", "16"], 
				//         ["mapq", "10"] ]
				// Don't check the validity of each cmd now. Execute one by one and if anything goes wrong
				// reset cmdConcatInput=null so that console input is asked again. Of course, what is executed is not 
				// rolled back.
				List<List<String>> cmdInputList= new ArrayList<List<String>>();

				for(String cmd : Splitter.on("&&").trimResults().omitEmptyStrings().split(cmdConcatInput)){
					cmdInputList.add(Utils.tokenize(cmd, " "));
				}
				
				for(List<String> cmdInput : cmdInputList){
					try {
						
						if(cmdInput.get(0).equals("h")){  
							Utils.printer(CommandList.briefHelp(), snapshotFile);
							currentCmdConcatInput= cmdConcatInput;
							cmdConcatInput= null;	
							
						} else if(cmdInput.size() >= 2 && cmdInput.get(1).equals("-h")){ // Help on this command
							Utils.printer("\n" + CommandList.getHelpForCommand(cmdInput.get(0)) + "\n", snapshotFile);
							currentCmdConcatInput= cmdConcatInput;
							cmdConcatInput= null;
						
						} else if(cmdInput.get(0).equals("history")){
							for(GenomicCoords xgc : gch.getHistory()){
								Utils.printer(xgc.toString() + "\n", snapshotFile);
							}
							currentCmdConcatInput= cmdConcatInput;
							cmdConcatInput= null;

						} else if(cmdInput.get(0).equals("showGenome")) {
							System.out.println(Utils.printSamSeqDict(gch.current().getSamSeqDict(), 30));
							currentCmdConcatInput= cmdConcatInput;
							cmdConcatInput= null;
						
						} else if(cmdInput.get(0).equals("infoTracks")) {
							System.out.println(trackSet.showTrackInfo());
							currentCmdConcatInput= cmdConcatInput;
							cmdConcatInput= null;
							
						} else if(cmdInput.get(0).equals("q")){
							System.exit(0);
						
						} else if(cmdInput.get(0).equals("f")
								|| cmdInput.get(0).equals("b")
								|| cmdInput.get(0).equals("ff") 
								|| cmdInput.get(0).equals("bb")
								|| cmdInput.get(0).matches("^\\d+.*")
								|| cmdInput.get(0).matches("^\\-\\d+.*") 
								|| cmdInput.get(0).matches("^\\+\\d+.*")){ // No cmd line args either f/b ops or ints
								// FIXME: You shouldn't join the list of args back to string. You should refactor parseConsoleInput! 
								String newRegion= Utils.parseConsoleInput(Joiner.on(" ").join(cmdInput), gch.current()).trim();
								gch.add(new GenomicCoords(newRegion, samSeqDict, windowSize, fasta));
								
						} else if(cmdInput.get(0).equals("goto") || cmdInput.get(0).startsWith(":")){
							String reg= Joiner.on(" ").join(cmdInput).replaceFirst("goto|:", "").trim();
							gch.add(new GenomicCoords(reg, samSeqDict, windowSize, fasta));
							
						} else if(cmdInput.get(0).equals("dataCol")){
							trackSet.setDataColForRegex(cmdInput);
											
						} else if(cmdInput.get(0).equals("ylim") && cmdInput.size() > 1){
							trackSet.setTrackYlimitsForRegex(cmdInput);
							
						} else if(cmdInput.get(0).equals("trackHeight") && cmdInput.size() > 1){
							trackSet.setTrackHeightForRegex(cmdInput);
							
						} else if((cmdInput.get(0).equals("colorTrack") || cmdInput.get(0).equals("colourTrack")) && cmdInput.size() > 1){
							trackSet.setTrackColourForRegex(cmdInput); 

						} else if(cmdInput.get(0).equals("BSseq")) {
							if( fasta == null ){
								System.err.println("Cannot set BSseq mode without fasta");
								cmdConcatInput= null;
								continue;
							}
						    trackSet.setBisulfiteModeForRegex(cmdInput);
						    
						} else if (cmdInput.get(0).equals("squash") || cmdInput.get(0).equals("merge")){
							trackSet.setFeatureDisplayModeForRegex(cmdInput);

						} else if (cmdInput.get(0).equals("gap")){
							trackSet.setFeatureGapForRegex(cmdInput);
							
						} else if(cmdInput.get(0).equals("gffNameAttr")) {
							trackSet.setAttributeForGFFName(cmdInput);
							
						} else if(cmdInput.get(0).equals("addTracks") && cmdInput.size() > 1){
							cmdInput.remove(0);
							Utils.addSourceName(inputFileList, cmdInput);
		
							if(gch.current().getSamSeqDict().size() == 0){
								samSeqDict = GenomicCoords.getSamSeqDictFromAnyFile(inputFileList, null, null);
								GenomicCoords gc= gch.current();
								gc.setSamSeqDict(samSeqDict);
							}
							
							for(String sourceName : cmdInput){
								trackSet.add(sourceName, gch.current());
							}
							
						} else if(cmdInput.get(0).equals("orderTracks") && cmdInput.size() > 1){
							cmdInput.remove(0);
							trackSet.orderTracks(cmdInput);
							
						} else if (cmdInput.get(0).equals("p")) {
							gch.previous();
							
						} else if (cmdInput.get(0).equals("n")) {
							gch.next();
							
						} else if(cmdInput.get(0).equals("zo")){
							int nz= Utils.parseZoom(Joiner.on(" ").join(cmdInput), 1);
							GenomicCoords gc = (GenomicCoords)gch.current().clone();
							for(int i= 0; i < nz; i++){
								gc.zoomOut();
							}
							gch.add(gc);
							
						} else if(cmdInput.get(0).equals("zi")){
							int nz= Utils.parseZoom(Joiner.on(" ").join(cmdInput), 1);
							GenomicCoords gc = (GenomicCoords)gch.current().clone();
							for(int i= 0; i < nz; i++){
								gc.zoomIn();
							}
							gch.add(gc);
						
						} else if(cmdInput.get(0).equals("print") || cmdInput.get(0).equals("printFull")){
								trackSet.setPrintModeForRegex(cmdInput);

						} else if(cmdInput.get(0).equals("next_start") || cmdInput.get(0).equals("next")){
							GenomicCoords gc= (GenomicCoords)gch.current().clone();
							String trackId= "";
							if(cmdInput.size() > 1){
								trackId= cmdInput.get(1);
							}
							if(cmdInput.get(0).equals("next_start")){
								gch.add(trackSet.goToNextFeatureOnFile(trackId, gc, -1.0));
							} else {
								gch.add(trackSet.goToNextFeatureOnFile(trackId, gc, 5.0));
							}

						} else if(cmdInput.get(0).equals("find_first") || 
								  cmdInput.get(0).equals("find_all")) {  
							if(cmdInput.size() < 2){
								System.err.println("Error in find* subcommand. Expected at least 1 args got: " + cmdInput);
								cmdConcatInput= null;
								continue;
							}
							if(cmdInput.size() == 2){
								cmdInput.add("");
							}
							GenomicCoords gc= (GenomicCoords)gch.current().clone();
							// Determine whether we match first or all
							boolean all= (cmdInput.get(0).equals("find_all")) ? true : false;
							gch.add(trackSet.findNextMatchOnTrack(cmdInput.get(1), cmdInput.get(2), gc, all));
						
						} else if (cmdInput.get(0).equals("seqRegex")){
							if( fasta == null ){
								System.err.println("Cannot find regex in sequence without fasta reference!");
								cmdConcatInput= null;
								continue;
							}
							String seqRegex= null;
							if(cmdInput.size() == 1){
								seqRegex= "";
								gch.setSeqRegex(seqRegex);
							} else {
								seqRegex= cmdInput.get(1);
								try{
									Pattern.compile(seqRegex);
									gch.setSeqRegex(seqRegex);
								} catch(PatternSyntaxException e){
							    	System.err.println("Invalid seqRegex in: " + cmdInput);
							    	System.err.println(e.getDescription());
									cmdConcatInput= null;
									continue;
								}						
							}

						} else if(cmdInput.get(0).equals("filter")){
							trackSet.setFilterForTrackIntervalFeature(cmdInput);

						} else if(cmdInput.get(0).equals("rpm")) {
							trackSet.setRpmForRegex(cmdInput);
							
						} else if(cmdInput.get(0).equals("-F")) {     //
							trackSet.setFilterFlagForRegex(cmdInput); //	
						}											  //
																	  //
						else if(cmdInput.get(0).equals("-f")) {       //
							trackSet.setFilterFlagForRegex(cmdInput); //  -F -f mapq are processed by the same method!
						}                                             //
																	  //
						else if(cmdInput.get(0).equals("mapq")) {     //
							trackSet.setFilterFlagForRegex(cmdInput); //	
						}                                             //
						
					    else if(cmdInput.get(0).equals("save")) {
							snapshotFile= Utils.parseCmdinputToGetSnapshotFile(Joiner.on(" ").join(cmdInput), gch.current());
						} else if(cmdInput.get(0).equals("bookmark")){
							// String name= cmdInput.replaceAll("^bookmark", "").trim();
							//GenomicCoords gc = (GenomicCoords)gch.current().clone();
							//trackSet.addBookmark_IN_PREP(gc, name);
							
						} else {
							System.err.println("Unrecognized argument: " + cmdInput);
							throw new InvalidCommandLineException();
						}
						
					} catch(Exception e){ // You shouldn't catch anything! Be more specific.
						System.err.println("\nError processing input: " + cmdInput + "\n");
						cmdConcatInput= null; 
					} // END PARSING ONE COMMAND
					
					if(cmdConcatInput == null){
						// If something goes wrong or help is invoked, stop executing commands and restart asking for input
						// Unless we are in non-interactive mode
						if(nonInteractive){
							System.exit(1);
						} 
						break;
					}
				} // END PARSING ALL COMMANDS
				if(cmdConcatInput != null){
					currentCmdConcatInput= cmdConcatInput;
				}
			} // END PARSING CONSOLE INPUT 
			
			int newSize= jline.TerminalFactory.get().getWidth() - 1;
			if(newSize != gch.current().getUserWindowSize()){
				// Replace the current genomicCoords obj with a new one having the same coordinates but different windowSize.
				// NB: The current genomic obj might not be the last one in the history list.
				windowSize= newSize;
				String newRegion= gch.current().getChrom() + ":" + gch.current().getFrom() + "-" + gch.current().getTo(); 
				gch.getHistory().add(gch.getHistory().indexOf(gch.current()), new GenomicCoords(newRegion, samSeqDict, windowSize, fasta));
			}
		} // End while loop keep going until quit or if no interactive input set
	}
	
}


/*
for(int i= 0; i < inputFileList.size(); i++){ 
	// Iterate through each input file
	String inputFileName= inputFileList.get(i);
	
	if(Utils.getFileTypeFromName(inputFileName).equals(TrackFormat.BAM)){
	
		if(!Utils.bamHasIndex(inputFileName)){
			System.err.println("\nNo index found for '" + inputFileName + "'. Index can be generated with ");
			System.err.println("samtools index '" + inputFileName + "'\n");
			System.exit(1);
		}
		
		// BAM Coverage track
		String coverageTrackId= new File(inputFileName).getName() + "#" + (idForTrack+1);
		idForTrack++;
		if(!trackSet.getTrackSet_DEPRECATED().containsKey(coverageTrackId)){
			TrackCoverage trackCoverage= new TrackCoverage(inputFileName, gch.current(), false);
			trackCoverage.setTrackTag(coverageTrackId);
			trackSet.getTrackSet_DEPRECATED().put(trackCoverage.getTrackTag(), trackCoverage);
		}
		TrackCoverage trackCoverage= (TrackCoverage) trackSet.getTrackSet_DEPRECATED().get(coverageTrackId);
		trackCoverage.setGc(gch.current());
		if(trackCoverage.getyMaxLines() > 0){
			trackCoverage.update();
		}
		trackCoverage.printToScreen();					

		// Reads
		String trackId= new File(inputFileName).getName() + "@" + (idForTrack+1);
		idForTrack++;
		if(!trackSet.getTrackSet_DEPRECATED().containsKey(trackId)){
			TrackReads trackReads= new TrackReads(inputFileName, gch.current());
			trackReads.setTrackTag(trackId);
			trackSet.getTrackSet_DEPRECATED().put(trackReads.getTrackTag(), trackReads);
			trackReads.setFilename(inputFileName);
			trackReads.setTrackTag(trackId);
		}
		TrackReads trackReads= (TrackReads) trackSet.getTrackSet_DEPRECATED().get(trackId);
		trackReads.setGc(gch.current());
		if(trackReads.getyMaxLines() > 0){
			trackReads.update();
		}
	} // End processing bam file
	
	// Annotatation
	if(    Utils.getFileTypeFromName(inputFileName).equals(TrackFormat.BED) 
        || Utils.getFileTypeFromName(inputFileName).equals(TrackFormat.GFF)
	    || Utils.getFileTypeFromName(inputFileName).equals(TrackFormat.VCF)){
		String trackId= new File(inputFileName).getName() + "#" + (idForTrack+1);
		idForTrack++;
		if(!trackSet.getTrackSet_DEPRECATED().containsKey(trackId)){
			TrackIntervalFeature tif= new TrackIntervalFeature(inputFileName, gch.current());
			tif.setTrackTag(trackId);
			trackSet.getTrackSet_DEPRECATED().put(tif.getTrackTag(), tif);
		}
		TrackIntervalFeature tif= (TrackIntervalFeature) trackSet.getTrackSet_DEPRECATED().get(trackId);
		tif.setGc(gch.current());
		try {
			if(tif.getyMaxLines() > 0){
				tif.update();
			}
		} catch(InvalidGenomicCoordsException e){
			e.printStackTrace();
		}
	} 
	// Wiggles
	if(Utils.getFileTypeFromName(inputFileName).equals(TrackFormat.BIGWIG) 
			|| Utils.getFileTypeFromName(inputFileName).equals(TrackFormat.TDF) 
			|| Utils.getFileTypeFromName(inputFileName).equals(TrackFormat.BEDGRAPH)){

		String trackId= new File(inputFileName).getName() + "#" + (idForTrack+1);
		idForTrack++;
		if(!trackSet.getTrackSet_DEPRECATED().containsKey(trackId)){
			TrackWiggles tw= new TrackWiggles(inputFileName, gch.current(), 4);
			tw.setTrackTag(trackId);
			trackSet.getTrackSet_DEPRECATED().put(tw.getTrackTag(), tw);
		}
		TrackWiggles tw= (TrackWiggles) trackSet.getTrackSet_DEPRECATED().get(trackId);
		tw.setGc(gch.current());
		if(tw.getyMaxLines() > 0){
			tw.update();
		}
		tw.printToScreen();
	} 
} // End loop through files 
*/