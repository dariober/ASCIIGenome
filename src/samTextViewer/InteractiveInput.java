package samTextViewer;

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
import htsjdk.samtools.SAMSequenceDictionary;

/** Class to process input from console
 * */
public class InteractiveInput {

	private boolean nonInteractive;
	private int interactiveInputExitCode= 0;
	private List<String> cmdHistory= new ArrayList<String>();
	public InteractiveInput(){
		
	}
	
	/** Parse the input list of commands to print information or modify the input TrackProcessor.  
	 * @throws IOException 
	 * @throws InvalidGenomicCoordsException 
	 * */
	protected TrackProcessor processInput(String cmdConcatInput, TrackProcessor proc) throws InvalidGenomicCoordsException, IOException{

		// cmdInputList: List of individual commands in tokens to be issued. 
		// E.g.: [ ["zi"], 
		//         ["-F", "16"], 
		//         ["mapq", "10"] ]
		// Don't check the validity of each cmd now. Execute one by one and if anything goes wrong
		// reset interactiveInputExitCode = 1 (or else other than 0) so that console input is asked again. Of course, what is executed is not 
		// rolled back.
		List<List<String>> cmdInputList= new ArrayList<List<String>>();
		
		for(String cmd : Splitter.on("&&").trimResults().omitEmptyStrings().split(cmdConcatInput)){
			cmdInputList.add(Utils.tokenize(cmd, " "));
		}

		int windowSize= proc.getWindowSize();
		String fasta= proc.getGenomicCoordsHistory().current().getFastaFile();
		SAMSequenceDictionary samSeqDict = proc.getGenomicCoordsHistory().current().getSamSeqDict();
		
		for(List<String> cmdInput : cmdInputList){
			// REMEMBER TO CALL TrackSet.update() AFTER METHODS THAT CHANGE THE DATA!!
			
			this.interactiveInputExitCode= 0; // If something goes wrong this will change to != 0
			try {
				
				// * These commands only print info or do stuff without editing the GenomicCoordinates or the Tracks:
				if(cmdInput.get(0).equals("h")){  
					Utils.printer(CommandList.briefHelp(), proc.getSnapshotFile());
					this.interactiveInputExitCode= 1;
					
				} else if(cmdInput.size() >= 2 && cmdInput.get(1).equals("-h")){ // Help on this command
					Utils.printer("\n" + CommandList.getHelpForCommand(cmdInput.get(0)) + "\n", proc.getSnapshotFile());
					this.interactiveInputExitCode= 1;
				
				} else if(cmdInput.get(0).equals("history")){
					for(GenomicCoords xgc : proc.getGenomicCoordsHistory().getHistory()){
						Utils.printer(xgc.toString() + "\n", proc.getSnapshotFile());
					}
					this.interactiveInputExitCode= 1;

				} else if(cmdInput.get(0).equals("cmdHistory")){
					Utils.printer(this.cmdHistoryToString() + "\n", proc.getSnapshotFile());
					this.interactiveInputExitCode= 1;
					
				} else if(cmdInput.get(0).equals("showGenome")) {
					System.out.println(Utils.printSamSeqDict(proc.getGenomicCoordsHistory().current().getSamSeqDict(), 30));
					this.interactiveInputExitCode= 1;
				
				} else if(cmdInput.get(0).equals("infoTracks")) {
					System.out.println(proc.getTrackSet().showTrackInfo());
					this.interactiveInputExitCode= 1;
				
				} else if(cmdInput.get(0).equals("save")) {
					proc.setSnapshotFile( Utils.parseCmdinputToGetSnapshotFile(Joiner.on(" ").join(cmdInput), proc.getGenomicCoordsHistory().current()) );

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
						// FIXME: You shouldn't join the list of args back to string. You should refactor parseConsoleInput! 
						String newRegion= Utils.parseConsoleInput(Joiner.on(" ").join(cmdInput), proc.getGenomicCoordsHistory().current()).trim();
						proc.getGenomicCoordsHistory().add(new GenomicCoords(newRegion, samSeqDict, windowSize, fasta));
						proc.getTrackSet().setGenomicCoordsAndUpdateTracks(proc.getGenomicCoordsHistory().current());
						
				} else if(cmdInput.get(0).equals("goto") || cmdInput.get(0).startsWith(":")){
					String reg= Joiner.on(" ").join(cmdInput).replaceFirst("goto|:", "").trim();
					proc.getGenomicCoordsHistory().add(new GenomicCoords(reg, samSeqDict, windowSize, fasta));
					proc.getTrackSet().setGenomicCoordsAndUpdateTracks(proc.getGenomicCoordsHistory().current());
					
				} else if (cmdInput.get(0).equals("p")) {
					proc.getGenomicCoordsHistory().previous();
					proc.getTrackSet().setGenomicCoordsAndUpdateTracks(proc.getGenomicCoordsHistory().current());
					
				} else if (cmdInput.get(0).equals("n")) {
					proc.getGenomicCoordsHistory().next();
					proc.getTrackSet().setGenomicCoordsAndUpdateTracks(proc.getGenomicCoordsHistory().current());
					
				} else if(cmdInput.get(0).equals("zo")){
					int nz= Utils.parseZoom(Joiner.on(" ").join(cmdInput), 1);
					GenomicCoords gc = (GenomicCoords)proc.getGenomicCoordsHistory().current().clone();
					for(int i= 0; i < nz; i++){
						gc.zoomOut();
					}
					proc.getGenomicCoordsHistory().add(gc);
					proc.getTrackSet().setGenomicCoordsAndUpdateTracks(proc.getGenomicCoordsHistory().current());
					
				} else if(cmdInput.get(0).equals("zi")){
					int nz= Utils.parseZoom(Joiner.on(" ").join(cmdInput), 1);
					GenomicCoords gc = (GenomicCoords)proc.getGenomicCoordsHistory().current().clone();
					for(int i= 0; i < nz; i++){
						gc.zoomIn();
					}
					proc.getGenomicCoordsHistory().add(gc);
					proc.getTrackSet().setGenomicCoordsAndUpdateTracks(proc.getGenomicCoordsHistory().current());
					
				// * These commands change the Tracks but do not touch the GenomicCoordinates.
				} else if(cmdInput.get(0).equals("dataCol")){
					proc.getTrackSet().setDataColForRegex(cmdInput);
					proc.getTrackSet().setGenomicCoordsAndUpdateTracks(proc.getGenomicCoordsHistory().current());
					
				} else if(cmdInput.get(0).equals("ylim") && cmdInput.size() > 1){
					proc.getTrackSet().setTrackYlimitsForRegex(cmdInput);
					proc.getTrackSet().setGenomicCoordsAndUpdateTracks(proc.getGenomicCoordsHistory().current());
					
				} else if(cmdInput.get(0).equals("trackHeight") && cmdInput.size() > 1){
					proc.getTrackSet().setTrackHeightForRegex(cmdInput);
					proc.getTrackSet().setGenomicCoordsAndUpdateTracks(proc.getGenomicCoordsHistory().current());
					
				} else if((cmdInput.get(0).equals("colorTrack") || cmdInput.get(0).equals("colourTrack")) && cmdInput.size() > 1){
					proc.getTrackSet().setTrackColourForRegex(cmdInput); 

				} else if(cmdInput.get(0).equals("hideTitle")){
					proc.getTrackSet().setHideTitleForRegex(cmdInput); 
					
				} else if(cmdInput.get(0).equals("BSseq")) {
					if( proc.getGenomicCoordsHistory().current().getFastaFile() == null ){
						System.err.println("Cannot set BSseq mode without fasta");
						this.interactiveInputExitCode= 1;
						continue;
					}
					proc.getTrackSet().setBisulfiteModeForRegex(cmdInput);
					proc.getTrackSet().setGenomicCoordsAndUpdateTracks(proc.getGenomicCoordsHistory().current());
					
				} else if (cmdInput.get(0).equals("squash") || cmdInput.get(0).equals("merge")){
					proc.getTrackSet().setFeatureDisplayModeForRegex(cmdInput);
					proc.getTrackSet().setGenomicCoordsAndUpdateTracks(proc.getGenomicCoordsHistory().current());
					
				} else if (cmdInput.get(0).equals("gap")){
					proc.getTrackSet().setFeatureGapForRegex(cmdInput);
					proc.getTrackSet().setGenomicCoordsAndUpdateTracks(proc.getGenomicCoordsHistory().current());
					
				} else if(cmdInput.get(0).equals("gffNameAttr")) {
					proc.getTrackSet().setAttributeForGFFName(cmdInput);
					proc.getTrackSet().setGenomicCoordsAndUpdateTracks(proc.getGenomicCoordsHistory().current());
					
				} else if(cmdInput.get(0).equals("addTracks") && cmdInput.size() > 1){
					cmdInput.remove(0);
					
					List<String> inputFileList = proc.getTrackSet().getFilenameList();
					
					Utils.addSourceName(inputFileList, cmdInput);

					if(proc.getGenomicCoordsHistory().current().getSamSeqDict().size() == 0){
						samSeqDict = GenomicCoords.getSamSeqDictFromAnyFile(inputFileList, null, null);
						GenomicCoords gc= proc.getGenomicCoordsHistory().current();
						gc.setSamSeqDict(samSeqDict);
					}
					
					for(String sourceName : cmdInput){
						proc.getTrackSet().add(sourceName, proc.getGenomicCoordsHistory().current());
					}
					
				} else if(cmdInput.get(0).equals("orderTracks")){
					cmdInput.remove(0);
					proc.getTrackSet().orderTracks(cmdInput);
									
				} else if(cmdInput.get(0).equals("print") || cmdInput.get(0).equals("printFull")){
					proc.getTrackSet().setPrintModeForRegex(cmdInput);

				} else if(cmdInput.get(0).equals("filter")){
					proc.getTrackSet().setFilterForTrackIntervalFeature(cmdInput);
					proc.getTrackSet().setGenomicCoordsAndUpdateTracks(proc.getGenomicCoordsHistory().current());
					
				} else if(cmdInput.get(0).equals("rpm")) {
					proc.getTrackSet().setRpmForRegex(cmdInput);
					proc.getTrackSet().setGenomicCoordsAndUpdateTracks(proc.getGenomicCoordsHistory().current());
					
				} else if(cmdInput.get(0).equals("-F")) {               //
					proc.getTrackSet().setFilterFlagForRegex(cmdInput); //	
					proc.getTrackSet().setGenomicCoordsAndUpdateTracks(proc.getGenomicCoordsHistory().current());          	    //
																		//
				} else if(cmdInput.get(0).equals("-f")) {       		//
					proc.getTrackSet().setFilterFlagForRegex(cmdInput); //  -F -f mapq are processed by the same method!
					proc.getTrackSet().setGenomicCoordsAndUpdateTracks(proc.getGenomicCoordsHistory().current());					//
																		//
				} else if(cmdInput.get(0).equals("mapq")) {     		//
					proc.getTrackSet().setFilterFlagForRegex(cmdInput); //
					proc.getTrackSet().setGenomicCoordsAndUpdateTracks(proc.getGenomicCoordsHistory().current());					//	

				// * These commands change both the Tracks and the GenomicCoordinates
				} else if(cmdInput.get(0).equals("next_start") || cmdInput.get(0).equals("next")){
					GenomicCoords gc= (GenomicCoords)proc.getGenomicCoordsHistory().current().clone();
					String trackId= "";
					if(cmdInput.size() > 1){
						trackId= cmdInput.get(1);
					}
					if(cmdInput.get(0).equals("next_start")){
						proc.getGenomicCoordsHistory().add(proc.getTrackSet().goToNextFeatureOnFile(trackId, gc, -1.0));
					} else {
						proc.getGenomicCoordsHistory().add(proc.getTrackSet().goToNextFeatureOnFile(trackId, gc, 5.0));
					}
					proc.getTrackSet().setGenomicCoordsAndUpdateTracks(proc.getGenomicCoordsHistory().current());
					
				} else if(cmdInput.get(0).equals("find_first") || 
						  cmdInput.get(0).equals("find_all")) {  
					if(cmdInput.size() < 2){
						System.err.println("Error in find* subcommand. Expected at least 1 args got: " + cmdInput);
						this.interactiveInputExitCode= 1;
						continue;
					}
					if(cmdInput.size() == 2){
						cmdInput.add("");
					}
					GenomicCoords gc= (GenomicCoords)proc.getGenomicCoordsHistory().current().clone();
					// Determine whether we match first or all
					boolean all= (cmdInput.get(0).equals("find_all")) ? true : false;
					proc.getGenomicCoordsHistory().add(proc.getTrackSet().findNextMatchOnTrack(cmdInput.get(1), cmdInput.get(2), gc, all));
					proc.getTrackSet().setGenomicCoordsAndUpdateTracks(proc.getGenomicCoordsHistory().current());
					
				} else if (cmdInput.get(0).equals("seqRegex")){
					if( proc.getGenomicCoordsHistory().current().getFastaFile() == null ){
						System.err.println("Cannot find regex in sequence without fasta reference!");
						this.interactiveInputExitCode= 1;
						continue;
					}
					String seqRegex= null;
					if(cmdInput.size() == 1){
						seqRegex= "";
					} else {
						seqRegex= cmdInput.get(1);
						try{
							Pattern.compile(seqRegex);
						} catch(PatternSyntaxException e){
					    	System.err.println("Invalid seqRegex in: " + cmdInput);
					    	System.err.println(e.getDescription());
					    	this.interactiveInputExitCode= 1;
							continue;
						}
					}
					proc.getTrackSet().setSeqRegexForTracks(seqRegex);
					proc.getTrackSet().setGenomicCoordsAndUpdateTracks(proc.getGenomicCoordsHistory().current());
					
				} else if(cmdInput.get(0).equals("bookmark")){
					// TODO
					
				} else {
					System.err.println("Unrecognized argument: " + cmdInput);
					throw new InvalidCommandLineException();
				}
				
			} catch(Exception e){ // You shouldn't catch anything! Be more specific.
				System.err.println("\nError processing input: " + cmdInput + "\n");
				this.interactiveInputExitCode= 1; 
			} // END PARSING ONE COMMAND
			
			if(this.interactiveInputExitCode != 0){
				// If something goes wrong or help is invoked, stop executing commands and restart asking for input
				// Unless we are in non-interactive mode
				if(nonInteractive){
					System.exit(1);
				} 
				break;
			}
		} // END OF LOOP THOURGH LIST OF INPUT
		return proc;
	}

	private String cmdHistoryToString() {
		List<String> cmd= new ArrayList<String>();
		int i = 1;
		for(String x : cmdHistory){
			cmd.add(i + ": \t" + x);
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

	protected List<String> getCmdHistory() {
		return cmdHistory;
	}

	protected void setCmdHistory(List<String> cmdHistory) {
		this.cmdHistory = cmdHistory;
	}
}
