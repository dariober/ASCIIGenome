package samTextViewer;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.MalformedURLException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.ListIterator;
import java.util.Set;

import org.apache.commons.lang3.StringUtils;

import com.google.common.base.Joiner;
import com.itextpdf.text.DocumentException;

import coloring.Config;
import coloring.ConfigKey;
import coloring.Xterm256;
import commandHelp.CommandHelp;
import commandHelp.CommandList;
import exceptions.BamIndexNotFoundException;
import exceptions.InvalidColourException;
import exceptions.InvalidCommandLineException;
import exceptions.InvalidConfigException;
import exceptions.InvalidGenomicCoordsException;
import exceptions.InvalidRecordException;
import faidx.UnindexableFastaFileException;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import jline.console.ConsoleReader;
import jline.console.completer.StringsCompleter;
import jline.console.history.History;
import jline.console.history.History.Entry;
import jline.console.history.MemoryHistory;
import net.sourceforge.argparse4j.inf.Namespace;
import tracks.IntervalFeature;
import tracks.Track;
import tracks.TrackFormat;
import tracks.TrackPileup;
import tracks.TrackSet;

/**
 * @author berald01
 *
 */
public class Main {
	
	public static void main(String[] args) throws IOException, InvalidGenomicCoordsException, InvalidCommandLineException, InvalidRecordException, BamIndexNotFoundException, ClassNotFoundException, SQLException, DocumentException, UnindexableFastaFileException, InvalidColourException, InvalidConfigException {

		final String CMD_HISTORY_FILE= System.getProperty("user.home") + File.separator + ".asciigenome_history";  
		final String MARKER_FOR_HISTORY_FILE= "## file ##";
		final String MARKER_FOR_HISTORY_CMD= "## cmd ##";
		/* Start parsing arguments * 
		 * *** If you change something here change also in console input ***/
		Namespace opts= ArgParse.argParse(args);
		
		List<String> initFileList= opts.getList("input");
		String region= opts.getString("region");
		final String fasta= opts.getString("fasta");
		String exec= opts.getString("exec");
		String config= opts.getString("config");
		exec= parseExec(exec);
		int debug= opts.getInt("debug");

		// Get configuration. Note that we don't need to assign this to a variable. 
		new Config(config);
		new Xterm256();
		
		// Init console right at start so if something goes wrong the user's terminal is reset to 
		// initial defaults with the shutdown hook. This could be achieved in cleaner way probably.
		ConsoleReader console = initConsole();
		
		messageVersion(opts.getBoolean("noFormat"));
		
		/* Set up console */
		
		Utils.checkFasta(fasta, debug);
		
		/* Test input files exist */
		List<String> inputFileList= new ArrayList<String>();
		Utils.addSourceName(inputFileList, initFileList, debug);
		
		/* Initialize trackSet */
		/* ------------------- */
		// This part only prepares a dummy GenomicCoords object to initialize the start position:
		// ----------------------------
		
		if(region == null || region.isEmpty()){
			region= initRegion(inputFileList, fasta, null, debug);			
		}		
		int terminalWidth= Utils.getTerminalWidth();
		GenomicCoords initGc= new GenomicCoords(region, terminalWidth, null, null);
		
		List<String>initGenomeList= new ArrayList<String>();
		for(String x : inputFileList){
			initGenomeList.add(x);
		}
		initGenomeList.add(fasta);
		initGc.setGenome(initGenomeList, false);
		// ----------------------------
		// Genomic positions start here:
		final GenomicCoordsHistory gch= new GenomicCoordsHistory();
		GenomicCoords start= new GenomicCoords(initGc.toStringRegion(), terminalWidth, initGc.getSamSeqDict(), initGc.getFastaFile());
		gch.readHistory(new File(CMD_HISTORY_FILE), start);
		gch.add(start);

		final TrackSet trackSet= new TrackSet(inputFileList, gch.current());
		addHistoryFiles(trackSet, CMD_HISTORY_FILE, MARKER_FOR_HISTORY_FILE, debug);
		
		setDefaultTrackHeights(console.getTerminal().getHeight(), trackSet.getTrackList());
		
		final TrackProcessor proc= new TrackProcessor(trackSet, gch);
		
		proc.setNoFormat(opts.getBoolean("noFormat"));
		
		// Put here the previous command so that it is re-issued if no input is given
		// You have to initialize this var outside the while loop that processes input files.
		String currentCmdConcatInput= ""; 

		if(!proc.isNoFormat()){
			String str= String.format("\033[48;5;%sm", Config.get256Color(ConfigKey.background));
			System.out.print(str);
		}

		// Batch processing file of regions
		final String batchFile= opts.getString("batchFile");
		if(batchFile != null && ! batchFile.isEmpty()){

			console.clearScreen();
			console.flush();

			BufferedReader br= batchFileReader(batchFile); // new BufferedReader(new FileReader(new File(batchFile)));
			String line = null;  
			while ((line = br.readLine()) != null){
				// Start processing intervals one by one
				IntervalFeature target= new IntervalFeature(line, TrackFormat.BED, null);
				String reg= target.getChrom() + ":" + target.getFrom() + "-" + target.getTo();
				String gotoAndExec= ("goto " + reg + " && " + exec).trim().replaceAll("&&$", "");
				InteractiveInput itr = new InteractiveInput(console);
				itr.processInput(gotoAndExec, proc, debug);
				if (itr.getInteractiveInputExitCode().equals(ExitCode.ERROR)){
					System.err.println("Error processing '" + gotoAndExec + "' at line '" + line + "'");
					System.exit(1);
				}
			}
			br.close();
			System.exit(0);
		}

		// See if we need to process the exec arg before going to interactive mode. 
		// Also if we are in non-interactive mode, we process the track set now and later exit 
		console.clearScreen();
		console.flush();		
		proc.iterateTracks();
		if(!exec.isEmpty() || opts.getBoolean("nonInteractive")){

			InteractiveInput itr = new InteractiveInput(console);
			itr.processInput(exec, proc, debug);
			if(opts.getBoolean("nonInteractive")){
				System.out.print("\033[0m");
				System.exit(0);
			}
		}

		/* Set up done, start processing */
		/* ============================= */
		History cmdHistory= initCmdHistory(CMD_HISTORY_FILE, MARKER_FOR_HISTORY_CMD, debug);
		console.setHistory(cmdHistory);
		writeHistory(console.getHistory(), trackSet.getOpenedFiles(), CMD_HISTORY_FILE, MARKER_FOR_HISTORY_FILE, MARKER_FOR_HISTORY_CMD, gch);
		
		while(true){  
			// keep going until quit or if no interactive input set
			// *** START processing interactive input
			String cmdConcatInput= ""; // String like "zi && -F 16 && mapq 10"
			InteractiveInput interactiveInput= new InteractiveInput(console);
			ExitCode currentExitCode= ExitCode.NULL;
			interactiveInput.setInteractiveInputExitCode(currentExitCode);
			
			while( ! interactiveInput.getInteractiveInputExitCode().equals(ExitCode.ERROR) 
					||  interactiveInput.getInteractiveInputExitCode().equals(ExitCode.NULL)){
				
				console.setPrompt(
						StringUtils.repeat(' ', proc.getWindowSize()) + '\r' + "[h] for help: "
						);

				cmdConcatInput= console.readLine().trim();
				if (cmdConcatInput.isEmpty()) {
					// Empty inout: User only issued <ENTER> 
					if( interactiveInput.getInteractiveInputExitCode().equals(ExitCode.CLEAN)){
						// User only issued <ENTER>: Repeat previous command if the exit code was not an error.
						cmdConcatInput= currentCmdConcatInput;					
					} else {
						// Refresh screen if the exit code was not CLEAN.
						cmdConcatInput= "+0";
					}
				}
				interactiveInput.processInput(cmdConcatInput, proc, debug);
				currentCmdConcatInput= cmdConcatInput;
			}
			// *** END processing interactive input 
		}
	}

	/** Return a suitable region to start. If a region is already given, do nothing.
	 * This method is a mess and should be cleaned up together with GenomicCoords class.
	 * @throws InvalidGenomicCoordsException 
	 * */
	public static String initRegion(List<String> inputFileList, String fasta, String genome, int debug ) throws IOException, InvalidGenomicCoordsException{
		// Preferably we start from a position that has a feature rather than from the start of a 
		// random chrom.
		
		System.err.print("Initializing coordinates... ");
		// First search for files that can init chrom and position
		List<String> skipped= new ArrayList<String>();
		for(String x : inputFileList){
			TrackFormat fmt = Utils.getFileTypeFromName(x);
			if(fmt.equals(TrackFormat.TDF)){
				skipped.add(x);
				continue;
			}
			try {
				String region= Utils.initRegionFromFile(x);
				System.err.println("Done from: " + x);
				return region;
			} catch(Exception e){
				System.err.println("\nCould not initilize from file " + x);
				if(debug > 0){
					e.printStackTrace();
				}
			}
		}
		// Try to initialize from fasta
		if(fasta != null && ! fasta.trim().isEmpty()){ 
			IndexedFastaSequenceFile faSeqFile = new IndexedFastaSequenceFile(new File(fasta));
			String region= faSeqFile.nextSequence().getName();
			faSeqFile.close();
			return region;
		}
		// Try genome file
		if(genome != null && ! genome.trim().isEmpty()){
			GenomicCoords gc= new GenomicCoords(Utils.getTerminalWidth());
			gc.setGenome(Arrays.asList(new String[] {genome}), false);
			SAMSequenceDictionary samSeqDict = gc.getSamSeqDict();
			String region= samSeqDict.getSequence(0).getSequenceName();
			return region;
		}
		// Failing that, look for any file that gives at least chrom
		for(String x : skipped){
			try {
				String region= Utils.initRegionFromFile(x);
				System.err.println("Done from: " + x);
				return region;
			} catch(Exception e){
				System.err.println("\nCould not initilize from file " + x);
				if(debug > 0){
					e.printStackTrace();
				}
			}
		}
		// It appears everything failed to initialise...
		return "";
//		List<String>initGenomeList= new ArrayList<String>(inputFileList);
		
		/* Prepare genomic coordinates to fetch. This should probably be a function in itself */
		// Create a dummy gc object just to get the sequence dict.
//		GenomicCoords gc= new GenomicCoords(Utils.getTerminalWidth());
//		
//		if(genome != null && ! genome.trim().isEmpty()){
//			initGenomeList.add(genome);
//		}
//		gc.setGenome(initGenomeList, false);
//		
//		SAMSequenceDictionary samSeqDict = gc.getSamSeqDict();

//		if(samSeqDict != null && ! samSeqDict.isEmpty()){
//			region= samSeqDict.getSequence(0).getSequenceName();
//			System.err.println("");
//		} else {
//			for(String x : initGenomeList){
//				try {
//					region= Utils.initRegionFromFile(x);
//					System.err.println("Done from: " + x);
//					break;
//				} catch(Exception e){
//					System.err.println("\nCould not initilize from file " + x);
//					if(debug > 0){
//						e.printStackTrace();
//					}
//				}
//			}
//		}
//		return region;
	}
	
	/** If exec is a file, parse it to return a string suitable for execution.  
	 * @throws IOException 
	 * */
	private static String parseExec(String exec) throws IOException{
		
		if(exec == null){
			return "";
		}
		
		if(new File(exec).isFile()){
			BufferedReader br= new BufferedReader(new FileReader(new File(exec)));
			List<String> x = new ArrayList<String>();
			String line= "";
			while((line = br.readLine()) != null){
				x.add(line.trim());
			}
			br.close();
			return Joiner.on(" && ").skipNulls().join(x);
		} else {
			return exec;
		}
	}

	private static BufferedReader batchFileReader(String batchFile) throws FileNotFoundException {

		if(batchFile.equals("-")){
			return new BufferedReader(new InputStreamReader(System.in));
		}
		
		if(! new File(batchFile).exists()){
			System.err.print("\033[0m");
			System.err.println("File " + batchFile + " does not exist.");
			System.exit(1);
		}
		
		return new BufferedReader(new FileReader(new File(batchFile)));
		
	}

	/**Read the asciigenome history file and put it a list as current history. Or
	 * return empty history file does not exist or can't be read. 
	 * */
	private static History initCmdHistory(String cmdHistoryFile, String MARKER_FOR_HISTORY_CMD, int debug){
		History cmdHistory= new MemoryHistory();
		try{
			BufferedReader br= new BufferedReader(new FileReader(new File(cmdHistoryFile)));
			String line;
			while((line = br.readLine()) != null){
				if(line.startsWith(MARKER_FOR_HISTORY_CMD)){
					cmdHistory.add(line.replaceFirst(MARKER_FOR_HISTORY_CMD, ""));
				}
			}
			br.close();
		} catch(IOException e){
			if(debug > 0){
				e.printStackTrace();
			}
		}
		return cmdHistory;
	}	

	/** Merge set of opened files in trackSet with the files found in the history file.
	 * */
	private static void addHistoryFiles(TrackSet trackSet, String cmdHistoryFile, String MARKER_FOR_HISTORY_FILE, int debug){
		LinkedHashSet<String> union= initRecentlyOpenedFiles(cmdHistoryFile, MARKER_FOR_HISTORY_FILE, debug);
		LinkedHashSet<String> now= trackSet.getOpenedFiles();
		for(String file : now){ // If a file is in the current track set and in the history file, put it last. I.e. last opened. 
			if(union.contains(file)){
				union.remove(file);
			}
			union.add(file);
		}
		trackSet.setOpenedFiles(union);
	}
	
	/**Read the asciigenome history file to extract the list of opened files 
	 * */
	private static LinkedHashSet<String> initRecentlyOpenedFiles(String cmdHistoryFile, 
			String MARKER_FOR_HISTORY_FILE, int debug){
		LinkedHashSet<String> opened= new LinkedHashSet<String>();
		try{
			BufferedReader br= new BufferedReader(new FileReader(new File(cmdHistoryFile)));
			String line;
			while((line = br.readLine()) != null){
				if(line.startsWith(MARKER_FOR_HISTORY_FILE)){
					line= line.substring(MARKER_FOR_HISTORY_FILE.length(), line.length()).trim();
					opened.add(line);
				}
			}
			br.close();
		} catch(IOException e){
			if(debug > 0){
				e.printStackTrace();
			}
		}
		return opened;
	}	
	
	/** Set some sensible defaults for track heights 
	 * @throws SQLException 
	 * @throws InvalidRecordException 
	 * @throws InvalidGenomicCoordsException 
	 * @throws IOException 
	 * @throws ClassNotFoundException 
	 * @throws MalformedURLException 
	 * */
	private static void setDefaultTrackHeights(int consoleHeight, List<Track> trackList) throws MalformedURLException, ClassNotFoundException, IOException, InvalidGenomicCoordsException, InvalidRecordException, SQLException{
		
		if(trackList.size() == 0){
			return;
		}

		if (trackList.get(0).getGc().getFastaFile() != null && trackList.get(0).getGc().isSingleBaseResolution){
			consoleHeight -= 1; // Reference sequence
		}
		if (trackList.get(0).getGc().getSamSeqDict() != null){
			consoleHeight -= 1; // Chrom ideogram
		}
		consoleHeight -= 1; // Region info
		consoleHeight -= 1; // prompt
		consoleHeight -= trackList.size(); // Track headers

		int consensus= 0; // Additional line for possible consensus sequence in TrackCoverage
		if(trackList.get(0).getGc().isSingleBaseResolution){
			consensus= 1;
		}
		for(Track tr : trackList){
			if(tr instanceof TrackPileup){
				consoleHeight -= consensus;
			}
		}
		
		if(consoleHeight <= 0){
			for(int i= 0; i < trackList.size(); i++){
				trackList.get(i).setyMaxLines(1);
			}			
			return;
		}
				
		// This is the list of heights that will be set at the end:
		List<Integer> trackHeights= new ArrayList<Integer>();
		for(int i= 0; i < trackList.size(); i++){
			trackHeights.add(0);
		}
		
		while(true){
			for(int i= 0; i < trackList.size(); i++){
				
				int h= trackHeights.get(i);
				h += 1;
				consoleHeight -= 1;
				
				if(consoleHeight <= 0){
					break;
				}
				trackHeights.set(i, h);
			}
			if(consoleHeight <= 0){
				break;
			}
		}
		
		for(int i= 0; i < trackList.size(); i++){
			if(trackHeights.get(i) < 20){
				trackList.get(i).setyMaxLines(trackHeights.get(i));
			} else {
				trackList.get(i).setyMaxLines(20);
			}
		}
	}
	
	/**Write the history of commands to given file.
	 * Note that the existing history file is completely overwritten.
	 * */
	private static void writeHistory(final History cmdHistory, final Set<String> fileHistory, final String cmdHistoryFile, final String MARKER_FOR_HISTORY_FILE, 
				final String MARKER_FOR_HISTORY_CMD, final GenomicCoordsHistory gch){
		
		Runtime.getRuntime().addShutdownHook(new Thread(new Runnable() {
			public void run() {
	            
				// List of commands
				ListIterator<Entry> iter = cmdHistory.entries();
				List<String>lastHist= new ArrayList<String>();
				int max_cmds= 2000; // Maximum number of commands to write out to asciigenomo_history.
				while(iter.hasNext()){
					if(max_cmds == 0){
						break;
					}
					max_cmds--;
					lastHist.add(iter.next().value().toString());
				}
				
				// List of files
				List<String>lastFiles= new ArrayList<String>(fileHistory);
				int max_files= 200; // Maximum number of files to write out to asciigenomo_history.
				lastFiles= lastFiles.subList(Math.max(0, lastFiles.size() - max_files), lastFiles.size());

				try{
					
					StringBuilder sb= new StringBuilder();
					
					for(String cmd : lastHist){
						sb.append(MARKER_FOR_HISTORY_CMD + cmd + "\n");
					}
					for(String file : lastFiles){
						sb.append(MARKER_FOR_HISTORY_FILE + file + "\n");
					}
					for(String pos : gch.prepareHistoryForHistoryFile(new File(cmdHistoryFile), 100)){
						sb.append(pos + "\n"); // Positions
					}

					BufferedWriter wr= new BufferedWriter(new FileWriter(new File(cmdHistoryFile)));
					wr.write(sb.toString());
					wr.close();
				} catch (IOException e) {
					System.err.println("Unable to write history to " + cmdHistoryFile);
				}
	        }
	    }, "Shutdown-thread"));
	}
		
	/** On exit print a message informing a new version of ASCIIGenome is available
	 * */
	private static void messageVersion(final boolean noFormat) throws IOException, InvalidColourException{
		new Xterm256();
		Runtime.getRuntime().addShutdownHook(new Thread(new Runnable() {
			public void run() {
				try{
					List<String> up = Utils.checkUpdates(5000);
					int cmp= Utils.versionCompare(up.get(0), up.get(1));
					// cmp= -1; // For testing
					String msg= "";
					if(cmp == -1){
						msg= "NEW: ASCIIGenome version v" + up.get(1) + " is available at " + ArgParse.WEB_ADDRESS + "/releases\n";
					}
					if( ! noFormat){
						msg= "\033[48;5;231;38;5;" + Xterm256.colorNameToXterm256("red") + "m" +  msg + "\033[0m\n";
					}
					System.err.print(msg);
				} catch(Exception e){
					// e.printStackTrace();
				}
	        }
	    }, "Shutdown-thread"));		
	}

	public static ConsoleReader initConsole() throws IOException, InvalidColourException{
		
		Runtime.getRuntime().addShutdownHook(new Thread(new Runnable() {
		    public void run() {
		    	System.out.print("\033[0m"); // On exit turn off all formatting
		    }
		}));
		
		ConsoleReader console= new ConsoleReader(); 

		try {
			// Autcomplete commands with length > x 
			for(CommandHelp x : CommandList.commandHelpList()){
				if(x.getName().length() > 2){
					console.addCompleter(new StringsCompleter(x.getName()));
				}
			}
		} catch (InvalidCommandLineException e) {
			e.printStackTrace();
		}
		console.setExpandEvents(false);
		return console;
	}
	
}