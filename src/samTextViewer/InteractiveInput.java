package samTextViewer;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;
import java.util.regex.Pattern;

import com.google.common.base.Joiner;
import com.google.common.base.Splitter;

import coloring.Config;
import commandHelp.Command;
import commandHelp.CommandList;
import exceptions.InvalidCommandLineException;
import exceptions.InvalidGenomicCoordsException;
import exceptions.InvalidRecordException;
import htsjdk.samtools.SAMSequenceDictionary;
import jline.console.ConsoleReader;
import jline.console.history.History.Entry;
import net.sourceforge.argparse4j.inf.ArgumentParserException;
import tracks.Track;

/** Class to process input from console
 * */
public class InteractiveInput {

	private boolean nonInteractive;
	private ExitCode interactiveInputExitCode= ExitCode.CLEAN;
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
	protected TrackProcessor processInput(String cmdConcatInput, TrackProcessor proc, boolean debug) throws InvalidGenomicCoordsException, IOException, ClassNotFoundException, InvalidRecordException, SQLException, InvalidCommandLineException{

		// cmdInputList: List of individual commands in tokens to be issued. 
		// E.g.: [ ["zi"], 
		//         ["-F", "16"], 
		//         ["mapq", "10"] ]
		// Don't check the validity of each cmd now. Execute one by one and if anything goes wrong
		// reset interactiveInputExitCode = 1 (or else other than 0) so that console input is asked again. Of course, what is executed is not 
		// rolled back.
		List<String> cmdInputChainList= new ArrayList<String>();
		
		// See http://stackoverflow.com/questions/1757065/java-splitting-a-comma-separated-string-but-ignoring-commas-in-quotes
		// For splitting at delimiter (&&) and ignore delimiters inside single quotes.
		for(String cmd : Splitter.on(Pattern.compile("&&(?=([^']*'[^']*')*[^']*$)")).trimResults().omitEmptyStrings().split(cmdConcatInput)){
			cmdInputChainList.add(cmd);
		}
		if(cmdInputChainList.size() >= 2 && cmdInputChainList.get(0).startsWith("setConfig ")){
			cmdInputChainList.add("+0"); // This is to refresh the screen and actually set the new color
		}

		String fasta= proc.getGenomicCoordsHistory().current().getFastaFile();
		SAMSequenceDictionary samSeqDict = proc.getGenomicCoordsHistory().current().getSamSeqDict();

		String messages= ""; // Messages that may be sent from the various methods.
		for(String cmdString : cmdInputChainList){
			
			List<String> cmdTokens= Utils.tokenize(cmdString, " ");
			
			this.interactiveInputExitCode= ExitCode.CLEAN; // If something goes wrong this will change
			try {
				
				// * These commands only print info or do stuff without editing the GenomicCoordinates or the Tracks:
				if(cmdTokens.get(0).equals("h") || cmdTokens.get(0).equals("-h")){  		
					System.err.println(Utils.padEndMultiLine(CommandList.briefHelp(), proc.getWindowSize()));
					this.interactiveInputExitCode= ExitCode.CLEAN_NO_FLUSH;
					
				} else if(cmdTokens.size() >= 2 && cmdTokens.get(1).equals("-h")){ // Help on this command
					String help= Utils.padEndMultiLine("\n" + CommandList.getHelpForCommand(cmdTokens.get(0)), proc.getWindowSize());
					System.err.println(help);
					this.interactiveInputExitCode= ExitCode.CLEAN_NO_FLUSH;
				
				} else if(cmdTokens.get(0).equals("posHistory")){
					this.posHistory(cmdTokens, proc.getGenomicCoordsHistory().getHistory(), proc.getWindowSize());
					this.interactiveInputExitCode= ExitCode.CLEAN_NO_FLUSH;

				} else if(cmdTokens.get(0).equals("history")){
					String hist= Utils.padEndMultiLine(this.cmdHistoryToString(cmdTokens), proc.getWindowSize());
					System.err.println(hist);
					this.interactiveInputExitCode= ExitCode.CLEAN_NO_FLUSH;
					
				} else if(cmdTokens.get(0).equals("showGenome")) {
					this.showGenome(proc);
					this.interactiveInputExitCode= ExitCode.CLEAN_NO_FLUSH;
				
				} else if(cmdTokens.get(0).equals("sys")) {
					this.execSysCmd(cmdString, proc.getWindowSize());
					this.interactiveInputExitCode= ExitCode.CLEAN_NO_FLUSH;
					
				} else if(cmdTokens.get(0).equals("infoTracks")) {
					String info= Utils.padEndMultiLine(proc.getTrackSet().showTrackInfo(), proc.getWindowSize());
					System.out.println(info);
					this.interactiveInputExitCode= ExitCode.CLEAN_NO_FLUSH;
				
				} else if(cmdTokens.get(0).equals("recentlyOpened")) {
					String opened= Utils.padEndMultiLine(proc.getTrackSet().showRecentlyOpened(cmdTokens), proc.getWindowSize());
					System.out.println(opened);
					this.interactiveInputExitCode= ExitCode.CLEAN_NO_FLUSH;
				
				} else if(cmdTokens.get(0).equals("setConfig")) {
					try{
						new Config(cmdTokens.get(1));
						this.interactiveInputExitCode= ExitCode.CLEAN;
					} catch(Exception e){
						System.err.println(Utils.padEndMultiLine("Unable to set configuration", proc.getWindowSize()));
						this.interactiveInputExitCode= ExitCode.ERROR;
						if(debug){
							e.printStackTrace();
						}
					}
					
				} else if(cmdTokens.get(0).equals("save")) {
					List<String> args= new ArrayList<String>(cmdTokens);
					
					proc.setStripAnsi( ! Utils.argListContainsFlag(args, "-color"));
					
					proc.setAppendToSnapshotFile(false); // Default: do not append
					
					if(args.contains(">>")){
						proc.setAppendToSnapshotFile(true);
						args.remove(">>");
					} else if (args.contains(">")){
						proc.setAppendToSnapshotFile(false);
						args.remove(">");
					}  
					proc.setSnapshotFile( Utils.parseCmdinputToGetSnapshotFile(Joiner.on(" ").join(args), proc.getGenomicCoordsHistory().current()) );
					
				} else if(cmdTokens.get(0).equals("sessionSave")) {
					if(cmdTokens.size() < 2){
						System.err.println(Utils.padEndMultiLine("Output file name is missing!", proc.getWindowSize()));
						this.interactiveInputExitCode= ExitCode.ERROR;
					} else {
						proc.exportTrackSetSettings(cmdTokens.get(1));
					}

				} else if(cmdTokens.get(0).equals("q")){
					System.exit(0);
				
				// * These commands change the GenomicCoordinates (navigate) but do not touch the tracks.
				} else if(cmdTokens.get(0).equals("f")
						|| cmdTokens.get(0).equals("b")
						|| cmdTokens.get(0).equals("ff") 
						|| cmdTokens.get(0).equals("bb")
						|| cmdTokens.get(0).matches("^\\d+.*")
						|| cmdTokens.get(0).matches("^\\-\\d+.*") 
						|| cmdTokens.get(0).matches("^\\+\\d+.*")){ // No cmd line args either f/b ops or ints
						String newRegion= Utils.parseConsoleInput(cmdTokens, proc.getGenomicCoordsHistory().current()).trim();
						proc.getGenomicCoordsHistory().add(new GenomicCoords(newRegion, samSeqDict, fasta));
						
				} else if(cmdTokens.get(0).equals("goto") || cmdTokens.get(0).startsWith(":")){
					String reg= Joiner.on(" ").join(cmdTokens).replaceFirst("goto|:", "").trim();
					proc.getGenomicCoordsHistory().add(new GenomicCoords(reg, samSeqDict, fasta));
					
				} else if (cmdTokens.get(0).equals("p")) {
					proc.getGenomicCoordsHistory().previous(); 
					
				} else if (cmdTokens.get(0).equals("n")) {
					proc.getGenomicCoordsHistory().next();
					
				} else if(cmdTokens.get(0).equals("zo")){
					int nz= Utils.parseZoom(Joiner.on(" ").join(cmdTokens), 1);
					GenomicCoords gc = (GenomicCoords)proc.getGenomicCoordsHistory().current().clone();
					for(int i= 0; i < nz; i++){
						gc.zoomOut();
					}
					proc.getGenomicCoordsHistory().add(gc);
					
				} else if(cmdTokens.get(0).equals("zi")){
					int nz= Utils.parseZoom(Joiner.on(" ").join(cmdTokens), 1);
					GenomicCoords gc = (GenomicCoords)proc.getGenomicCoordsHistory().current().clone();
					for(int i= 0; i < nz; i++){
						gc.zoomIn();
					}
					proc.getGenomicCoordsHistory().add(gc);

				} else if(cmdTokens.get(0).equals("extend")) {
					if(cmdTokens.size() == 1){
						System.err.println(Utils.padEndMultiLine("Expected at least one argument.", proc.getWindowSize()));
					}
					GenomicCoords gc = (GenomicCoords)proc.getGenomicCoordsHistory().current().clone();
					gc.cmdInputExtend(cmdTokens);
					proc.getGenomicCoordsHistory().add(gc);
				
				} else if(cmdTokens.get(0).equals("trim")){
					GenomicCoords gc = proc.getTrackSet().trimCoordsForTrack(cmdTokens);
					proc.getGenomicCoordsHistory().add(gc);
					
				} else if(cmdTokens.get(0).equals("l")) {
					GenomicCoords gc = (GenomicCoords)proc.getGenomicCoordsHistory().current().clone();
					gc.left();
					proc.getGenomicCoordsHistory().add(gc);
				
				} else if(cmdTokens.get(0).equals("r")) {
					GenomicCoords gc = (GenomicCoords)proc.getGenomicCoordsHistory().current().clone();
					gc.right();
					proc.getGenomicCoordsHistory().add(gc);
					
				} else if(cmdTokens.get(0).equals("setGenome")){
					this.setGenome(cmdTokens, proc);
					
				// * These commands change the Tracks but do not touch the GenomicCoordinates.
				} else if(cmdTokens.get(0).equals("dataCol")){
					try{
						proc.getTrackSet().setDataColForRegex(cmdTokens);
					} catch(Exception e){
						String msg= Utils.padEndMultiLine("Error processing " + cmdTokens + ". Perhaps a non-numeric column was selected?", proc.getWindowSize());
						System.err.println(msg);
						this.interactiveInputExitCode= ExitCode.ERROR;
						continue;
					}
					
				} else if(cmdTokens.get(0).equals("ylim")){
					proc.getTrackSet().setTrackYlimitsForRegex(cmdTokens);
					
				} else if(cmdTokens.get(0).equals("trackHeight")){
					proc.getTrackSet().setTrackHeightForRegex(cmdTokens);
					// proc.getTrackSet().setGenomicCoordsAndUpdateTracks(proc.getGenomicCoordsHistory().current());
					
				} else if((cmdTokens.get(0).equals("colorTrack") || cmdTokens.get(0).equals("colourTrack"))){
					proc.getTrackSet().setTrackColourForRegex(cmdTokens); 

				} else if(cmdTokens.get(0).equals("hideTitle")){
					proc.getTrackSet().setHideTitleForRegex(cmdTokens); 
					
				} else if(cmdTokens.get(0).equals(Command.BSseq.getCmdDescr())) {
					if( proc.getGenomicCoordsHistory().current().getFastaFile() == null ){
						String msg= Utils.padEndMultiLine("Cannot set BSseq mode without reference sequence", proc.getWindowSize());
						System.err.println(msg);
						this.interactiveInputExitCode= ExitCode.ERROR;
						continue;
					}
					proc.getTrackSet().setBisulfiteModeForRegex(cmdTokens);
					
				} else if (cmdTokens.get(0).equals("squash") || cmdTokens.get(0).equals(Command.featureDisplayMode.toString())){
					proc.getTrackSet().setFeatureDisplayModeForRegex(cmdTokens);
					
				} else if (cmdTokens.get(0).equals("gap")){
					proc.getTrackSet().setFeatureGapForRegex(cmdTokens);
					
				} else if(cmdTokens.get(0).equals("gffNameAttr")) {
					proc.getTrackSet().setAttributeForGFFName(cmdTokens);
					
				} else if(cmdTokens.get(0).equals("addTracks")){
					cmdTokens.remove(0);
					
					List<String> globbed = Utils.globFiles(cmdTokens);
					if(globbed.size() == 0){
						String msg= Utils.padEndMultiLine(cmdTokens + ": No file found.", proc.getWindowSize());
						System.err.println(msg);
						this.interactiveInputExitCode= ExitCode.ERROR;
					
					} else {
						
						for(String sourceName : globbed){
							String msg= Utils.padEndMultiLine("Adding: " + sourceName, proc.getWindowSize());
							System.err.println(msg);
							try{
								proc.getTrackSet().addTrackFromSource(sourceName, proc.getGenomicCoordsHistory().current(), null);
							} catch(Exception e){
								try{
									// It may be that you are in position that doesn't exist in the sequence dictionary that
									// came with this new file. To recover, find an existing position, move there and try to reload the 
									// file. This fixes issue#23
									String region= Main.initRegion(null, globbed, null, null, debug);
									proc.getGenomicCoordsHistory().add(new GenomicCoords(region, samSeqDict, fasta));
									proc.getTrackSet().addTrackFromSource(sourceName, proc.getGenomicCoordsHistory().current(), null);							
								} catch (Exception x){
									msg= Utils.padEndMultiLine("Failed to add: " + sourceName, proc.getWindowSize());
									System.err.println(msg);
								}
							}
							
							if(proc.getGenomicCoordsHistory().current().getSamSeqDict() == null || proc.getGenomicCoordsHistory().current().getSamSeqDict().size() == 0){
								GenomicCoords gc= proc.getGenomicCoordsHistory().current();
								// We are adding tracks. Check if we can set genome but do not treat these files as genome files.
								gc.setGenome(Arrays.asList(new String[] {sourceName}), false);
							}
							
						}
					}
					
				} else if(cmdTokens.get(0).equals("dropTracks")){
					if(cmdTokens.size() <= 1){
						System.err.println(Utils.padEndMultiLine("List one or more tracks to drop or `dropTracks -h` for help.", proc.getWindowSize()));
						this.interactiveInputExitCode= ExitCode.ERROR;
						continue;
					}
					messages += proc.getTrackSet().dropTracksWithRegex(cmdTokens);
					
				} else if(cmdTokens.get(0).equals("orderTracks")){
					cmdTokens.remove(0);
					proc.getTrackSet().orderTracks(cmdTokens);
				
				} else if(cmdTokens.get(0).equals("editNames")){
					messages += proc.getTrackSet().editNamesForRegex(cmdTokens);
					
				} else if(cmdTokens.get(0).equals(Command.print.toString())){
					proc.getTrackSet().setPrintModeAndPrintFeaturesForRegex(cmdTokens);

				} else if(cmdTokens.get(0).equals("grep")){
					proc.getTrackSet().setFilterForTrackIntervalFeature(cmdTokens);
				
				} else if(cmdTokens.get(0).equals("awk")){
					proc.getTrackSet().setAwkForTrackIntervalFeature(cmdTokens);
					
				} else if(cmdTokens.get(0).equals(Command.rpm.getCmdDescr())) {
					proc.getTrackSet().setRpmForRegex(cmdTokens);

				} else if(cmdTokens.get(0).equals("samtools")){
					proc.getTrackSet().setSamFilterForRegex(cmdTokens);

				// * These commands change both the Tracks and the GenomicCoordinates
				} else if(cmdTokens.get(0).equals("next")){
					
					this.next(cmdTokens, proc);
										
				} else if(cmdTokens.get(0).equals("find")) {  
					// Determine whether we match first or all
					boolean all= false;
					if(cmdTokens.contains("-all")){
						all= true;
						cmdTokens.remove("-all");
					}

					if(cmdTokens.size() < 2){
						System.err.println(Utils.padEndMultiLine("Error in find command. Expected at least 1 argument got: " + cmdTokens, proc.getWindowSize()));
						this.interactiveInputExitCode= ExitCode.ERROR;
						continue;
					}
					if(cmdTokens.size() == 2){
						cmdTokens.add(""); // If track arg is missing use this placeholder.
					}
					GenomicCoords gc= (GenomicCoords)proc.getGenomicCoordsHistory().current().clone();
					proc.getGenomicCoordsHistory().add(proc.getTrackSet().findNextMatchOnTrack(cmdTokens.get(1), cmdTokens.get(2), gc, all));
					
				} else if (cmdTokens.get(0).equals("seqRegex")){
					try{
						proc.getTrackSet().setRegexForTrackSeqRegex(cmdTokens, proc.getGenomicCoordsHistory().current());
					} catch( InvalidCommandLineException e){
						System.err.println(Utils.padEndMultiLine("Cannot find regex in sequence without fasta reference!", proc.getWindowSize()));
						this.interactiveInputExitCode= ExitCode.ERROR;
						continue;						
					}
					
				} else if(cmdTokens.get(0).equals("bookmark")){
					messages += proc.getTrackSet().bookmark(proc.getGenomicCoordsHistory().current(), cmdTokens);
					
				} else {
					System.err.println(Utils.padEndMultiLine("Unrecognized command: " + cmdTokens, proc.getWindowSize()));
					this.interactiveInputExitCode= ExitCode.ERROR;
					throw new InvalidCommandLineException();
				}
			} catch(ArgumentParserException e){
				this.interactiveInputExitCode= ExitCode.ERROR;
			
			} catch(Exception e){ // You shouldn't catch anything! Be more specific.
				System.err.println(Utils.padEndMultiLine("\nError processing input: " + cmdTokens, proc.getWindowSize()));
				System.err.println(Utils.padEndMultiLine("For help on command \"cmd\" execute 'cmd -h' or '-h' for list of commands.\n", proc.getWindowSize()));
				this.interactiveInputExitCode= ExitCode.ERROR; 
				if(debug){
					e.printStackTrace();
				}
			} // END PARSING ONE COMMAND

			if( this.interactiveInputExitCode.equals(ExitCode.CLEAN) || this.interactiveInputExitCode.equals(ExitCode.CLEAN_NO_FLUSH)){
				// Command has been parsed ok. Let's see if we can execute it without exceptions.
				try{
					if(this.interactiveInputExitCode.equals(ExitCode.CLEAN)){
						console.clearScreen();
						console.flush();
						proc.iterateTracks();
					} else {
						//
					}

				} catch (InvalidGenomicCoordsException e){
					
					String region= Main.initRegion(null, proc.getTrackSet().getFilenameList(), null, null, debug);
					proc.getGenomicCoordsHistory().add(new GenomicCoords(region, samSeqDict, fasta));
					System.err.println(Utils.padEndMultiLine("Invalid genomic coordinates found. Resetting to "  + region, proc.getWindowSize()));
					if(debug){
						e.printStackTrace();
					}
					
				} catch (Exception e){
					System.err.println(Utils.padEndMultiLine("Error processing tracks with input " + cmdTokens, proc.getWindowSize()));
					this.interactiveInputExitCode= ExitCode.ERROR;
					if(debug){
						e.printStackTrace();
					}
				}
			}
			if(this.interactiveInputExitCode.equals(ExitCode.ERROR)) {
				// If something goes wrong or help is invoked, stop executing commands and restart asking for input
				// Unless we are in non-interactive mode
				if(nonInteractive){
					System.exit(1);
				} 
				break;
			}
		} // END OF LOOP THROUGH CHAIN OF INPUT COMMANDS
		if( ! messages.isEmpty()){
			System.err.println(Utils.padEndMultiLine(messages.trim(), proc.getWindowSize())); 
		}
		messages= "";
		return proc;
	}

    
	private void setGenome(List<String> cmdTokens, TrackProcessor proc) throws InvalidGenomicCoordsException, IOException, InvalidCommandLineException {
		
		List<String> tokens = new ArrayList<String>(cmdTokens); 
		tokens.remove(0);
		
		if(tokens.size() == 0){
			throw new InvalidCommandLineException();
		}
		
//		if(tokens.get(0).equals("-a")){
//			System.err.println("HERE");
//			// Collect genomes in resource
//			List<String> files = IOUtils.readLines(Main.class.getResource("/genomes/"));		
//			for(String x : files){
//				if(x.endsWith(".genome")){
//					System.err.println(Utils.padEndMultiLine(x.replaceAll(".genome$", ""), proc.getWindowSize()));	
//				}
//			}
//			this.interactiveInputExitCode= 1;
//			return;
//		}

		GenomicCoords testSeqDict= new GenomicCoords("default", null, null); 
		testSeqDict.setGenome(tokens, true);
		if(testSeqDict.getSamSeqDict() != null){
			proc.getGenomicCoordsHistory().setGenome(tokens);
		} else {
			System.err.println(Utils.padEndMultiLine("Cannot set genome from " + tokens, proc.getWindowSize()));
			this.interactiveInputExitCode= ExitCode.ERROR;
		}
		
	}

	/** Execute arbitrary system command and print its output
	 *  
	 * @param cmdInput: String, in contrast to other coomands, process the raw string, not
	 * the tokenized version so you don't mess up with single quotes inside the 
	 * system command.
	 */
	private void execSysCmd(String cmdInput, int userWindowSize) {

		String rawSysCmd= cmdInput.trim().replaceAll("^sys +", ""); // Remove command name
		boolean isLiteral= false;
		if(rawSysCmd.trim().startsWith("-L ")){
			// Check if -L option is present and remove it if yes.
			isLiteral= true;
			rawSysCmd= rawSysCmd.trim().replaceAll("-L +", "");
		}
		if(rawSysCmd.isEmpty()){
			System.err.println(Utils.padEndMultiLine("Please provide a system comand to execute. Use `sys -h` for help.", userWindowSize));
			return;
		}
		
		List<String> tokens= new ArrayList<String>();
		
		if( isLiteral ){
			tokens.add(rawSysCmd);			
		} else { // w/o -L option execute command as bash string 
			tokens.add("bash");
			tokens.add("-c");
			tokens.add(rawSysCmd);
		}
		
		try {
			Process p = new ProcessBuilder().inheritIO().command(tokens).start();
			p.waitFor();
			BufferedReader reader= new BufferedReader(new InputStreamReader(p.getInputStream()));

            String line = "";
			while ((line = reader.readLine())!= null) {
				System.err.println(Utils.padEndMultiLine(line, userWindowSize));
			}

		} catch (Exception e) {
			System.err.println(Utils.padEndMultiLine(e.getMessage(), userWindowSize));
		}

	}

	private void posHistory(List<String> cmdInput, List<GenomicCoords> history, int userWindowSize) throws InvalidCommandLineException {
		
		List<String> args= new ArrayList<String>(cmdInput);
		args.remove(0); // Remove cmd name.
		
		int n= 10;
		if(args.contains("-n")){
			try{
				n= Integer.parseInt(args.get(args.indexOf("-n") + 1));
				args.remove(args.get(args.indexOf("n") + 1));
				args.remove("-n");
			} catch (Exception e){
				System.err.println(Utils.padEndMultiLine("Argument to -n parameter must be an integer", userWindowSize));
				throw new InvalidCommandLineException();
			}
		}

		int start= 0; // Start listing positions from this index
		if(history.size() > n && n > 0){
			start= history.size() - n;
		}
		
		for(int i= start; i < history.size(); i++){
			GenomicCoords xgc= history.get(i);
			System.err.println(Utils.padEndMultiLine(xgc.toString(), userWindowSize));
		}
	}

	/** Move to next feature using parameters in cmdInput. 
	 * First arg in cmdInput is command name itself. 
	 * The side effect is to modify the TrackProcessor obj to update the position.
	 * */
	private void next(List<String> cmdInput, TrackProcessor proc) throws InvalidGenomicCoordsException, IOException, InvalidCommandLineException {
		
		List<String> args= new ArrayList<String>(cmdInput); 
		args.remove(0); // Remove command name
		
		int zo= 5;
		if(args.contains("-zo")){
			try{
				zo= Integer.parseInt(args.get(args.indexOf("-zo") + 1));
				args.remove(args.get(args.indexOf("-zo") + 1));
				args.remove("-zo");
				if(zo < 0){
					zo= 0;
				}
			} catch (Exception e){
				System.err.println(Utils.padEndMultiLine("Argument to -zo parameter must be an integer", proc.getWindowSize()));
				throw new InvalidCommandLineException();
			}
		}
		
		GenomicCoords gc= (GenomicCoords)proc.getGenomicCoordsHistory().current().clone();
		String trackId= "";
		boolean start= false;
		if(args.contains("-start")){
			start= true;
			args.remove("-start");
		}
		boolean getPrevious= false;
		if(args.contains("-back")){
			getPrevious= true;
			args.remove("-back");
		}
		if(args.size() > 0){
			trackId= args.get(0);
		}
		if(start){
			proc.getGenomicCoordsHistory().add(proc.getTrackSet().goToNextFeatureOnFile(trackId, gc, -1.0, getPrevious));
		} else {
			proc.getGenomicCoordsHistory().add(proc.getTrackSet().goToNextFeatureOnFile(trackId, gc, zo, getPrevious));
		}
		
	}

	private void showGenome(TrackProcessor proc) throws InvalidGenomicCoordsException, IOException {
		String genome= Utils.printSamSeqDict(proc.getGenomicCoordsHistory().current().getSamSeqDict(), 30);
		if(genome != null && ! genome.isEmpty()){
			System.err.println(Utils.padEndMultiLine(genome, proc.getWindowSize()));
			return;
		}
		Set<String> chroms= new TreeSet<String>();
		for(Track tr : proc.getTrackSet().getTrackList()){
			chroms.addAll(tr.getChromosomeNames());
		}
		System.err.println(Utils.padEndMultiLine("Known contigs:", proc.getWindowSize()));
		for(String x : chroms){
			System.err.println(Utils.padEndMultiLine(x, proc.getWindowSize()));
		}
		return;
	}

	private String cmdHistoryToString(List<String> cmdInput) throws InvalidCommandLineException {
		
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
		List<String> cmd= new ArrayList<String>();
		int i = 1;
		for(Entry x : console.getHistory()){
			if(pattern.matcher(x.value().toString()).find()){
				cmd.add(i + ": \t" + x.value().toString());	
			}
			i++;
		}
		
		if(cmd.size() > nmax){ // Trim list to 
			cmd= cmd.subList(cmd.size() - nmax, cmd.size());
		}
		
		List<String> cmdTab= Utils.tabulateList(cmd);
		String tab= "";
		for(String x : cmdTab){
			tab += (x + "\n");
		}
		return tab.trim();
	}

	public ExitCode getInteractiveInputExitCode() {
		return interactiveInputExitCode;
	}

	public void setInteractiveInputExitCode(ExitCode exitCode) {
		this.interactiveInputExitCode= exitCode;
	}

}
