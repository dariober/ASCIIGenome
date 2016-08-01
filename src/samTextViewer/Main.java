package samTextViewer;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Pattern;
import java.util.regex.PatternSyntaxException;

import com.google.common.base.Joiner;
import com.google.common.base.Splitter;

import commandHelp.CommandList;
import exceptions.InvalidCommandLineException;
import exceptions.InvalidGenomicCoordsException;
import net.sourceforge.argparse4j.inf.Namespace;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import jline.console.ConsoleReader;
import tracks.Track;
import tracks.TrackCoverage;
import tracks.TrackFormat;
import tracks.TrackIntervalFeature;
import tracks.TrackReads;
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

	public static void main(String[] args) throws IOException, InvalidGenomicCoordsException, InvalidCommandLineException {
		
		/* Start parsing arguments * 
		 * *** If you change something here change also in console input ***/
		Namespace opts= ArgParse.argParse(args);
		
		List<String> inputFileList= opts.getList("input");
		String region= opts.getString("region");
		final String genome= opts.getString("genome");
		final String fasta= opts.getString("fasta");
		final String exec= opts.getString("exec");
		final int maxReadsStack= opts.getInt("maxReadsStack");
		final boolean noFormat= opts.getBoolean("noFormat");
		final boolean nonInteractive= opts.getBoolean("nonInteractive");
		// boolean withReadName= false; 
				
		int windowSize= 160;
		try{
			windowSize= jline.TerminalFactory.get().getWidth() - 1; 
		} catch(Exception e){
			e.printStackTrace();
		}
		
		Utils.checkFasta(fasta);
		
		/* Test input files exist */
		List<String> dropMe= new ArrayList<String>();
		for(String x : inputFileList){
			if(!new File(x).exists() && !Utils.urlFileExists(x)){
				dropMe.add(x);
			} 
		}
		for(String x : dropMe){
			System.err.println("\nWarning: Dropping file " + x + " as it does not exist.\n");
			inputFileList.remove(x);
		}
		
		if(inputFileList.size() == 0 && fasta == null){
			System.err.println("\nNo files in input: Nothing to be done!\n");
			System.exit(1);
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

		TrackSet trackSet= new TrackSet();

		/* Initialize GC profile */
		if(fasta != null){
			TrackWiggles cgWiggle= gch.current().getGCProfile();
			trackSet.getTrackSet().put(cgWiggle.getFileTag(), cgWiggle);
		}
		
		TrackIntervalFeature seqRegexTrack= gch.current().findRegex("^$");
		
		// Put here the previous command so that it is re-issued if no imput is given
		// You have to initialize this var outside the while loop that processes input files.
		String currentCmdConcatInput= ""; 
		String seqRegex= null;
		int idForTrack= 0;
		String snapshotFile= null;
		boolean snapshotStripAnsi= true;
		ConsoleReader console = CommandList.initConsole();
		
		boolean execDone= exec.isEmpty() ? true : false; // Do we need to execute commands from --exec? If exec is empty, consider it done.

		/* =================================== *
		 * Start processing interactive input  *
		 * =================================== */
		
		while(true){ 

			for(int i= 0; i < inputFileList.size(); i++){ 
				/* Iterate through each input file */
				String inputFileName= inputFileList.get(i);
				
				if(Utils.getFileTypeFromName(inputFileName).equals(TrackFormat.BAM)){
				
					/* Coverage and methylation track */
					String coverageTrackId= new File(inputFileName).getName() + "#" + (idForTrack+1);
					idForTrack++;
					if(!trackSet.getTrackSet().containsKey(coverageTrackId)){
						TrackCoverage trackCoverage= new TrackCoverage(inputFileName, gch.current(), false);
						trackCoverage.setFileTag(coverageTrackId);
						trackSet.getTrackSet().put(trackCoverage.getFileTag(), trackCoverage);
					}
					TrackCoverage trackCoverage= (TrackCoverage) trackSet.getTrackSet().get(coverageTrackId);
					trackCoverage.setGc(gch.current());
					if(trackCoverage.getyMaxLines() > 0){
						trackCoverage.update();
					}
					trackCoverage.printToScreen();					
					/* Methylation profile disable until a better representation is prepared
				
					if(bs && trackCoverage.getScreenLocusInfoList().size() > 0){
						if(maxMethylLines < 0){
							maxMethylLines= 0;
						}
						coverageTrackId= new File(inputFileName).getName() + "#" + (idForTrack+1);
						idForTrack++;
						if(!trackSet.getTrackSet().containsKey(coverageTrackId)){
							TrackMethylation trackMethylation= new TrackMethylation(inputFileName, trackCoverage.getScreenLocusInfoList());
							trackMethylation.setFileTag(coverageTrackId);
							trackSet.getTrackSet().put(trackMethylation.getFileTag(), trackMethylation);
						}
						TrackMethylation trackMethylation= (TrackMethylation) trackSet.getTrackSet().get(coverageTrackId);
						trackMethylation.setScreenLocusInfoList(trackCoverage.getScreenLocusInfoList());
						trackMethylation.setyMaxLines(maxMethylLines);
					}
					*/
										
					/* Reads */
					String trackId= new File(inputFileName).getName() + "@" + (idForTrack+1);
					idForTrack++;
					if(!trackSet.getTrackSet().containsKey(trackId)){
						TrackReads trackReads= new TrackReads(inputFileName, gch.current(), maxReadsStack);
						trackReads.setFileTag(trackId);
						trackSet.getTrackSet().put(trackReads.getFileTag(), trackReads);
						trackReads.setFilename(inputFileName);
						trackReads.setFileTag(trackId);
					}
					TrackReads trackReads= (TrackReads) trackSet.getTrackSet().get(trackId);
					trackReads.setGc(gch.current());
					// trackReads.setWithReadName(withReadName);
					if(trackReads.getyMaxLines() > 0){
						trackReads.update();
					}
				} // End processing bam file
				
				/* Annotatation */
				if(    Utils.getFileTypeFromName(inputFileName).equals(TrackFormat.BED) 
			        || Utils.getFileTypeFromName(inputFileName).equals(TrackFormat.GFF)
				    || Utils.getFileTypeFromName(inputFileName).equals(TrackFormat.VCF)){
					String trackId= new File(inputFileName).getName() + "#" + (idForTrack+1);
					idForTrack++;
					if(!trackSet.getTrackSet().containsKey(trackId)){
						TrackIntervalFeature tif= new TrackIntervalFeature(inputFileName, gch.current());
						tif.setFileTag(trackId);
						trackSet.getTrackSet().put(tif.getFileTag(), tif);
					}
					TrackIntervalFeature tif= (TrackIntervalFeature) trackSet.getTrackSet().get(trackId);
					tif.setGc(gch.current());
					try {
						if(tif.getyMaxLines() > 0){
							tif.update();
						}
					} catch(InvalidGenomicCoordsException e){
						e.printStackTrace();
					}
				} 
				/* Wiggles */
				if(Utils.getFileTypeFromName(inputFileName).equals(TrackFormat.BIGWIG) 
						|| Utils.getFileTypeFromName(inputFileName).equals(TrackFormat.TDF) 
						|| Utils.getFileTypeFromName(inputFileName).equals(TrackFormat.BEDGRAPH)){

					String trackId= new File(inputFileName).getName() + "#" + (idForTrack+1);
					idForTrack++;
					if(!trackSet.getTrackSet().containsKey(trackId)){
						TrackWiggles tw= new TrackWiggles(inputFileName, gch.current(), 4);
						tw.setFileTag(trackId);
						trackSet.getTrackSet().put(tw.getFileTag(), tw);
					}
					TrackWiggles tw= (TrackWiggles) trackSet.getTrackSet().get(trackId);
					tw.setGc(gch.current());
					if(tw.getyMaxLines() > 0){
						tw.update();
					}
					tw.printToScreen();
				}
			} // End loop through files 

			/* Print tracks */
			/* ************ */
			if(!nonInteractive){
				console.clearScreen();
				console.flush();
			}
			
			if(gch.current().getChromIdeogram(20) != null && execDone){
				Utils.printer(gch.current().getChromIdeogram(20) + "\n", snapshotFile, snapshotStripAnsi);
			}			
			for(Track tr : trackSet.getTrackSet().values()){
				if(tr.getFileTag() == gch.current().getGcProfileFileTag()){
					continue;
				}
				tr.setNoFormat(noFormat);
				if(execDone && tr.getyMaxLines() > 0){
					Utils.printer(tr.getTitle(), snapshotFile, snapshotStripAnsi);
					Utils.printer(tr.printToScreen() + "\n", snapshotFile, snapshotStripAnsi);
					Utils.printer(tr.getPrintableConsensusSequence(), snapshotFile, snapshotStripAnsi);
				}

				// Print features
				String printable= tr.printFeatures(windowSize);
				if(execDone){
					Utils.printer(printable, snapshotFile, snapshotStripAnsi);
				}
			}

			/* Footers and interactive prompt */
			/* ****************************** */
			// GC content profile
			if(trackSet.getTrackSet().containsKey(gch.current().getGcProfileFileTag())){
				// There is a bad hack here: We get yMaxLines from the TrackSet, but we don't plot the wiggle profile in the track set. 
				// Instead we get wiggle from GenomicCoords object and replace yMaxLines in GenomicCoords with the one from the TrackSet.
				// This is because the TrackSet profile is not updated.
				int yMaxLines= trackSet.getTrackSet().get(gch.current().getGcProfileFileTag()).getyMaxLines();
				double yLimitMin= trackSet.getTrackSet().get(gch.current().getGcProfileFileTag()).getYLimitMin();
				double yLimitMax= trackSet.getTrackSet().get(gch.current().getGcProfileFileTag()).getYLimitMax();
				String col= trackSet.getTrackSet().get(gch.current().getGcProfileFileTag()).getTitleColour();
				if(yMaxLines > 0){
					TrackWiggles tw= gch.current().getGCProfile();				
					tw.setyMaxLines(yMaxLines);
					tw.setYLimitMin(yLimitMin);
					tw.setYLimitMax(yLimitMax);
					tw.setTitleColour(col);
					tw.setNoFormat(noFormat);
					// tw.setHidden(trackSet.getTrackSet().get(GenomicCoords.gcProfileFileTag).isHidden());
					if(execDone && tw.getyMaxLines() > 0){
						Utils.printer(tw.getTitle(), snapshotFile, snapshotStripAnsi);
						String gcPrintable= tw.printToScreen();
						Utils.printer(gcPrintable + "\n", snapshotFile, snapshotStripAnsi);
					}
				}
			}
			// Track for matching regex
			int prevHeight= seqRegexTrack.getyMaxLines();
			seqRegexTrack = gch.current().findRegex(seqRegex);
			seqRegexTrack.setNoFormat(noFormat);
			seqRegexTrack.setyMaxLines(prevHeight);
			// See if we need to change height
			for(Pattern p : trackSet.getRegexForTrackHeight()){
				if(p.matcher(seqRegexTrack.getFileTag()).find()){
					seqRegexTrack.setyMaxLines(trackSet.getTrackHeightForRegex());
					break;
				}
			}
			String seqPattern= seqRegexTrack.printToScreen();
			if(!seqPattern.isEmpty()){
				seqPattern+="\n";
			} 
			if(execDone){
				Utils.printer(seqPattern, snapshotFile, snapshotStripAnsi);
			}
			// Ruler and sequence
			if(execDone){
				Utils.printer(gch.current().printableRefSeq(noFormat), snapshotFile, snapshotStripAnsi);
				String ruler= gch.current().printableRuler(10);
				Utils.printer(ruler.substring(0, ruler.length() <= windowSize ? ruler.length() : windowSize) + "\n", snapshotFile, snapshotStripAnsi);
			}
			String footer= gch.current().toString() + "; " + Math.rint(gch.current().getBpPerScreenColumn() * 10d)/10d + " bp/char; " + getMemoryStat();
			if(execDone){
				if(!noFormat){
					Utils.printer("\033[0;34m" + footer + "\033[0m; \n", snapshotFile, snapshotStripAnsi);
				} else {
					Utils.printer(footer + "\n", snapshotFile, snapshotStripAnsi);
				}
				// Optionally convert to png
				if(snapshotFile != null && snapshotFile.endsWith("png")){
					Utils.convertTextFileToGraphic(new File(snapshotFile), new File(snapshotFile));
				} 
			}
			
			/* Interactive input */
			/* ================= */
			if(nonInteractive && execDone){
				// Exit now if non-interctive mode is set and no string is passed to be executed. 
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
							Utils.printer(CommandList.briefHelp(), snapshotFile, snapshotStripAnsi);
							currentCmdConcatInput= cmdConcatInput;
							cmdConcatInput= null;	
							
						} else if(cmdInput.size() >= 2 && cmdInput.get(1).equals("-h")){ // Help on this command
							Utils.printer("\n" + CommandList.getHelpForCommand(cmdInput.get(0)) + "\n", snapshotFile, snapshotStripAnsi);
							currentCmdConcatInput= cmdConcatInput;
							cmdConcatInput= null;
						
						//} else if(cmdInput.get(0).equals("windowSize")){
						//	if(cmdInput.size() == 1){
						//		windowSize= jline.TerminalFactory.get().getWidth() - 1;
						//	} else {
						//		windowSize= Integer.parseInt(cmdInput.get(1));
						//	}
						//	// Replace the current genomicCoords obj with a new one having the same coordinates 
						//	// but different windowSize.
						//	// NB: The current genomic obj might not be the last one in the history list.
						//	String newRegion= gch.current().getChrom() + ":" + gch.current().getFrom() + "-" + gch.current().getTo(); 
						//	gch.getHistory().add(gch.getHistory().indexOf(gch.current()), new GenomicCoords(newRegion, samSeqDict, windowSize, fasta));							
							
						} else if(cmdInput.get(0).equals("history")){
							for(GenomicCoords xgc : gch.getHistory()){
								Utils.printer(xgc.toString() + "\n", snapshotFile, snapshotStripAnsi);
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
							Utils.addTrack(inputFileList, cmdInput);
				
							if(gch.current().getSamSeqDict().size() == 0){
								samSeqDict = GenomicCoords.getSamSeqDictFromAnyFile(inputFileList, null, null);
								GenomicCoords gc= gch.current();
								gc.setSamSeqDict(samSeqDict);
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
							if(cmdInput.size() == 1){
								seqRegex= "";
							} else {
								seqRegex= cmdInput.get(1);
								try{
									Pattern.compile(seqRegex);
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
						
					    else if(cmdInput.get(0).equals("save") || cmdInput.get(0).equals("savef")) {
							snapshotFile= Utils.parseCmdinputToGetSnapshotFile(Joiner.on(" ").join(cmdInput), gch.current());
							if(cmdInput.get(0).equals("save")){
								snapshotStripAnsi= true;
							} else if(cmdInput.get(0).equals("savef")){
								snapshotStripAnsi= false;
							}
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
			
			idForTrack= 0;
		} // End while loop keep going until quit or if no interactive input set
	}
	
}
