package samTextViewer;

import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

import com.google.common.base.Joiner;
import com.google.common.base.Splitter;

import commandHelp.Command;
import commandHelp.CommandList;
import exceptions.InvalidCommandLineException;
import exceptions.InvalidGenomicCoordsException;
import exceptions.InvalidRecordException;
import htsjdk.samtools.SAMSequenceDictionary;
import jline.console.ConsoleReader;
import jline.console.history.History.Entry;
import tracks.Track;

/** Class to process input from console
 * */
public class InteractiveInput {

	private boolean nonInteractive;
	private int interactiveInputExitCode= 0;
	private ConsoleReader console;
	public InteractiveInput(ConsoleReader console){
		this.console= console;
	}
	
	/** Parse the input list of commands to print information or modify the input TrackProcessor.  
	 * @throws IOException 
	 * @throws InvalidGenomicCoordsException 
	 * @throws InvalidCommandLineException 
	 * @throws SQLException 
	 * @throws InvalidRecordException 
	 * @throws ClassNotFoundException 
	 * */
	protected TrackProcessor processInput(String cmdConcatInput, TrackProcessor proc) throws InvalidGenomicCoordsException, IOException, ClassNotFoundException, InvalidRecordException, SQLException, InvalidCommandLineException{

		// ConsoleReader console = CommandList.initConsole();
		
		// cmdInputList: List of individual commands in tokens to be issued. 
		// E.g.: [ ["zi"], 
		//         ["-F", "16"], 
		//         ["mapq", "10"] ]
		// Don't check the validity of each cmd now. Execute one by one and if anything goes wrong
		// reset interactiveInputExitCode = 1 (or else other than 0) so that console input is asked again. Of course, what is executed is not 
		// rolled back.
		List<List<String>> cmdInputChainList= new ArrayList<List<String>>();
		
		for(String cmd : Splitter.on("&&").trimResults().omitEmptyStrings().split(cmdConcatInput)){
			cmdInputChainList.add(Utils.tokenize(cmd, " "));
		}

		String fasta= proc.getGenomicCoordsHistory().current().getFastaFile();
		SAMSequenceDictionary samSeqDict = proc.getGenomicCoordsHistory().current().getSamSeqDict();

		String messages= ""; // Messages that may be sent from the various methods.
		for(List<String> cmdInput : cmdInputChainList){
			
			this.interactiveInputExitCode= 0; // If something goes wrong this will change to != 0
			try {
				
				// * These commands only print info or do stuff without editing the GenomicCoordinates or the Tracks:
				if(cmdInput.get(0).equals("h") || cmdInput.get(0).equals("-h")){  
					System.err.println(CommandList.briefHelp());
					this.interactiveInputExitCode= 1;
					
				} else if(cmdInput.size() >= 2 && cmdInput.get(1).equals("-h")){ // Help on this command
					System.err.println("\n" + CommandList.getHelpForCommand(cmdInput.get(0)));
					this.interactiveInputExitCode= 1;
				
				} else if(cmdInput.get(0).equals("posHistory")){
					for(GenomicCoords xgc : proc.getGenomicCoordsHistory().getHistory()){
						System.err.println(xgc.toString());
					}
					this.interactiveInputExitCode= 1;

				} else if(cmdInput.get(0).equals("history")){
					System.err.println(this.cmdHistoryToString());
					this.interactiveInputExitCode= 1;
					
				} else if(cmdInput.get(0).equals("showGenome")) {
					this.showGenome(proc);
					this.interactiveInputExitCode= 1;
				
				} else if(cmdInput.get(0).equals("infoTracks")) {
					System.out.println(proc.getTrackSet().showTrackInfo());
					this.interactiveInputExitCode= 1;
				
				} else if(cmdInput.get(0).equals("recentlyOpened")) {
					System.out.println(proc.getTrackSet().showRecentlyOpened());
					this.interactiveInputExitCode= 1;
					
				} else if(cmdInput.get(0).equals("save")) {
					List<String> args= new ArrayList<String>(cmdInput);
					
					proc.setAppendToSnapshotFile(false); // Default: do not append
					
					if(args.contains(">>")){
						proc.setAppendToSnapshotFile(true);
						args.remove(">>");
					} else if (args.contains(">")){
						proc.setAppendToSnapshotFile(false);
						args.remove(">");
					}  
					proc.setSnapshotFile( Utils.parseCmdinputToGetSnapshotFile(Joiner.on(" ").join(args), proc.getGenomicCoordsHistory().current()) );
					
				} else if(cmdInput.get(0).equals("sessionSave")) {
					if(cmdInput.size() < 2){
						System.err.println("Output file name is missing!");
						this.interactiveInputExitCode= 1;
					} else {
						proc.exportTrackSetSettings(cmdInput.get(1));
					}

				} else if(cmdInput.get(0).equals("q")){
					System.exit(0);
				
				// * These commands change the GenomicCoordinates (navigate) but do not touch the tracks.
				} else if(cmdInput.get(0).equals("f")
						|| cmdInput.get(0).equals("b")
						|| cmdInput.get(0).equals("ff") 
						|| cmdInput.get(0).equals("bb")
						|| cmdInput.get(0).matches("^\\d+.*")
						|| cmdInput.get(0).matches("^\\-\\d+.*") 
						|| cmdInput.get(0).matches("^\\+\\d+.*")){ // No cmd line args either f/b ops or ints
						String newRegion= Utils.parseConsoleInput(cmdInput, proc.getGenomicCoordsHistory().current()).trim();
						proc.getGenomicCoordsHistory().add(new GenomicCoords(newRegion, samSeqDict, fasta));
						
				} else if(cmdInput.get(0).equals("goto") || cmdInput.get(0).startsWith(":")){
					String reg= Joiner.on(" ").join(cmdInput).replaceFirst("goto|:", "").trim();
					proc.getGenomicCoordsHistory().add(new GenomicCoords(reg, samSeqDict, fasta));
					
				} else if (cmdInput.get(0).equals("p")) {
					proc.getGenomicCoordsHistory().previous();
					
				} else if (cmdInput.get(0).equals("n")) {
					proc.getGenomicCoordsHistory().next();
					
				} else if(cmdInput.get(0).equals("zo")){
					int nz= Utils.parseZoom(Joiner.on(" ").join(cmdInput), 1);
					GenomicCoords gc = (GenomicCoords)proc.getGenomicCoordsHistory().current().clone();
					for(int i= 0; i < nz; i++){
						gc.zoomOut();
					}
					proc.getGenomicCoordsHistory().add(gc);
					
				} else if(cmdInput.get(0).equals("zi")){
					int nz= Utils.parseZoom(Joiner.on(" ").join(cmdInput), 1);
					GenomicCoords gc = (GenomicCoords)proc.getGenomicCoordsHistory().current().clone();
					for(int i= 0; i < nz; i++){
						gc.zoomIn();
					}
					proc.getGenomicCoordsHistory().add(gc);

				} else if(cmdInput.get(0).equals("extend")) {
					if(cmdInput.size() == 1){
						System.err.println("Expected at least one argument.");
					}
					GenomicCoords gc = (GenomicCoords)proc.getGenomicCoordsHistory().current().clone();
					gc.cmdInputExtend(cmdInput);
					proc.getGenomicCoordsHistory().add(gc);
				
				} else if(cmdInput.get(0).equals("trim")){
					GenomicCoords gc = proc.getTrackSet().trimCoordsForTrack(cmdInput);
					proc.getGenomicCoordsHistory().add(gc);
					
				} else if(cmdInput.get(0).equals("l")) {
					GenomicCoords gc = (GenomicCoords)proc.getGenomicCoordsHistory().current().clone();
					gc.left();
					proc.getGenomicCoordsHistory().add(gc);
				
				} else if(cmdInput.get(0).equals("r")) {
					GenomicCoords gc = (GenomicCoords)proc.getGenomicCoordsHistory().current().clone();
					gc.right();
					proc.getGenomicCoordsHistory().add(gc);
					
				} else if(cmdInput.get(0).equals("setGenome")){
					cmdInput.remove(0);
					proc.getGenomicCoordsHistory().setGenome(cmdInput);

				// * These commands change the Tracks but do not touch the GenomicCoordinates.
				} else if(cmdInput.get(0).equals("dataCol")){
					try{
						proc.getTrackSet().setDataColForRegex(cmdInput);
					} catch(Exception e){
						System.err.println("Error processing " + cmdInput + ". Perhaps a non-numeric column was selected?");
						this.interactiveInputExitCode= 1;
						continue;
					}
					
				} else if(cmdInput.get(0).equals("ylim") && cmdInput.size() > 1){
					proc.getTrackSet().setTrackYlimitsForRegex(cmdInput);
					
				} else if(cmdInput.get(0).equals("trackHeight") && cmdInput.size() > 1){
					proc.getTrackSet().setTrackHeightForRegex(cmdInput);
					// proc.getTrackSet().setGenomicCoordsAndUpdateTracks(proc.getGenomicCoordsHistory().current());
					
				} else if((cmdInput.get(0).equals("colorTrack") || cmdInput.get(0).equals("colourTrack")) && cmdInput.size() > 1){
					proc.getTrackSet().setTrackColourForRegex(cmdInput); 

				} else if(cmdInput.get(0).equals("hideTitle")){
					proc.getTrackSet().setHideTitleForRegex(cmdInput); 
					
				} else if(cmdInput.get(0).equals(Command.BSseq.getCmdDescr())) {
					if( proc.getGenomicCoordsHistory().current().getFastaFile() == null ){
						System.err.println("Cannot set BSseq mode without reference sequence");
						this.interactiveInputExitCode= 1;
						continue;
					}
					proc.getTrackSet().setBisulfiteModeForRegex(cmdInput);
					
				} else if (cmdInput.get(0).equals("squash") || cmdInput.get(0).equals(Command.featureDisplayMode.toString())){
					proc.getTrackSet().setFeatureDisplayModeForRegex(cmdInput);
					
				} else if (cmdInput.get(0).equals("gap")){
					proc.getTrackSet().setFeatureGapForRegex(cmdInput);
					
				} else if(cmdInput.get(0).equals("gffNameAttr")) {
					proc.getTrackSet().setAttributeForGFFName(cmdInput);
					
				} else if(cmdInput.get(0).equals("addTracks") && cmdInput.size() > 1){
					cmdInput.remove(0);
					
					List<String> globbed = Utils.globFiles(cmdInput);
					if(globbed.size() == 0){
					
						System.err.println(cmdInput + ": No file found.");
						this.interactiveInputExitCode= 1;
					
					} else {
						
						for(String sourceName : globbed){
	
							System.err.println("Adding: " + sourceName);
							try{
								proc.getTrackSet().addTrackFromSource(sourceName, proc.getGenomicCoordsHistory().current(), null);
							} catch(Exception e){
								try{
									// It may be that you are in position that doesn't exist in the sequence dictionary that
									// came with this new file. To recover, find an existing position, move there and try to reload the 
									// file. This fixes issue#23
									String region= Main.initRegion(null, globbed, null, null);
									proc.getGenomicCoordsHistory().add(new GenomicCoords(region, samSeqDict, fasta));
									proc.getTrackSet().addTrackFromSource(sourceName, proc.getGenomicCoordsHistory().current(), null);							
								} catch (Exception x){
									System.err.println("Failed to add: " + sourceName);
								}
							}
							
							if(proc.getGenomicCoordsHistory().current().getSamSeqDict() == null || proc.getGenomicCoordsHistory().current().getSamSeqDict().size() == 0){
								GenomicCoords gc= proc.getGenomicCoordsHistory().current();
								gc.setGenome(sourceName);
							}
							
						}
					}
					
				} else if(cmdInput.get(0).equals("dropTracks")){
					if(cmdInput.size() <= 1){
						System.err.println("List one or more tracks to drop or `dropTracks -h` for help.");
						this.interactiveInputExitCode= 1;
						continue;
					}
					messages += proc.getTrackSet().dropTracksWithRegex(cmdInput);
					
				} else if(cmdInput.get(0).equals("orderTracks")){
					cmdInput.remove(0);
					proc.getTrackSet().orderTracks(cmdInput);
				
				} else if(cmdInput.get(0).equals("editNames")){
					messages += proc.getTrackSet().editNamesForRegex(cmdInput);
					
				} else if(cmdInput.get(0).equals(Command.print.toString())){
					proc.getTrackSet().setPrintModeAndPrintFeaturesForRegex(cmdInput);

				} else if(cmdInput.get(0).equals("grep")){
					proc.getTrackSet().setFilterForTrackIntervalFeature(cmdInput);
					
				} else if(cmdInput.get(0).equals(Command.rpm.getCmdDescr())) {
					proc.getTrackSet().setRpmForRegex(cmdInput);

				} else if(cmdInput.get(0).equals("samtools")){
					proc.getTrackSet().setSamFilterForRegex(cmdInput);

				// * These commands change both the Tracks and the GenomicCoordinates
				} else if(cmdInput.get(0).equals("next")){
					GenomicCoords gc= (GenomicCoords)proc.getGenomicCoordsHistory().current().clone();
					String trackId= "";
					boolean start= false;
					if(cmdInput.contains("-start")){
						start= true;
						cmdInput.remove("-start");
					}
					boolean getPrevious= false;
					if(cmdInput.contains("-back")){
						getPrevious= true;
						cmdInput.remove("-back");
					}
					if(cmdInput.size() > 1){
						trackId= cmdInput.get(1);
					}
					if(start){
						proc.getGenomicCoordsHistory().add(proc.getTrackSet().goToNextFeatureOnFile(trackId, gc, -1.0, getPrevious));
					} else {
						proc.getGenomicCoordsHistory().add(proc.getTrackSet().goToNextFeatureOnFile(trackId, gc, 5.0, getPrevious));
					}
					
				} else if(cmdInput.get(0).equals("find")) {  
					// Determine whether we match first or all
					boolean all= false;
					if(cmdInput.contains("-all")){
						all= true;
						cmdInput.remove("-all");
					}

					if(cmdInput.size() < 2){
						System.err.println("Error in find command. Expected at least 1 argument got: " + cmdInput);
						this.interactiveInputExitCode= 1;
						continue;
					}
					if(cmdInput.size() == 2){
						cmdInput.add(""); // If track arg is missing use this placeholder.
					}
					GenomicCoords gc= (GenomicCoords)proc.getGenomicCoordsHistory().current().clone();
					proc.getGenomicCoordsHistory().add(proc.getTrackSet().findNextMatchOnTrack(cmdInput.get(1), cmdInput.get(2), gc, all));
					
				} else if (cmdInput.get(0).equals("seqRegex")){
					try{
						proc.getTrackSet().setRegexForTrackSeqRegex(cmdInput, proc.getGenomicCoordsHistory().current());
					} catch( InvalidCommandLineException e){
						System.err.println("Cannot find regex in sequence without fasta reference!");
						this.interactiveInputExitCode= 1;
						continue;						
					}
					
					
				} else if(cmdInput.get(0).equals("bookmark")){
					messages += proc.getTrackSet().bookmark(proc.getGenomicCoordsHistory().current(), cmdInput);
					
				} else {
					System.err.println("Unrecognized command: " + cmdInput);
					throw new InvalidCommandLineException();
				}
			
			} catch(Exception e){ // You shouldn't catch anything! Be more specific.
				System.err.println("\nError processing input: " + cmdInput + "\n");
				this.interactiveInputExitCode= 1; 
			} // END PARSING ONE COMMAND

			if(this.interactiveInputExitCode == 0){
				// Command has been parsed ok. Let's see if we can execute it without exceptions.
				try{
					console.clearScreen();
					console.flush();
					proc.iterateTracks();

				} catch (InvalidGenomicCoordsException e){
					
					String region= Main.initRegion(null, proc.getTrackSet().getFilenameList(), null, null);
					proc.getGenomicCoordsHistory().add(new GenomicCoords(region, samSeqDict, fasta));
					System.err.println("Invalid genomic coordinates found. Resetting to "  + region);
					
				} catch (Exception e){
					System.err.println("Error processing tracks with input " + cmdInput);
					this.interactiveInputExitCode= 1;
				}
			}
			if(this.interactiveInputExitCode != 0) {
				// If something goes wrong or help is invoked, stop executing commands and restart asking for input
				// Unless we are in non-interactive mode
				if(nonInteractive){
					System.exit(1);
				} 
				break;
			}
		} // END OF LOOP THROUGH CHAIN OF INPUT COMMANDS
		System.err.print(messages); 
		messages= "";
		return proc;
	}

	private void showGenome(TrackProcessor proc) {
		String genome= Utils.printSamSeqDict(proc.getGenomicCoordsHistory().current().getSamSeqDict(), 30);
		if(genome != null && ! genome.isEmpty()){
			System.err.println(genome);
			return;
		}
		Set<String> chroms= new TreeSet<String>();
		for(Track tr : proc.getTrackSet().getTrackList()){
			chroms.addAll(tr.getChromosomeNames());
		}
		System.err.println("Known contigs:");
		for(String x : chroms){
			System.err.println(x);	
		}
		return;
	}

	private String cmdHistoryToString() {
		List<String> cmd= new ArrayList<String>();
		int i = 1;
		for(Entry x : console.getHistory()){
			cmd.add(i + ": \t" + x.value().toString());
			i++;
		}
		List<String> cmdTab= Utils.tabulateList(cmd);
		String tab= "";
		for(String x : cmdTab){
			tab += (x + "\n");
		}
		return tab.trim();
	}

	public int getInteractiveInputExitCode() {
		return interactiveInputExitCode;
	}

	public void setInteractiveInputExitCode(int exitCode) {
		this.interactiveInputExitCode= exitCode;
	}

//	protected List<String> getCmdHistory() {
//		return cmdHistory;
//	}
//
//	protected void setCmdHistory(List<String> cmdHistory) {
//		this.cmdHistory = cmdHistory;
//	}
}
