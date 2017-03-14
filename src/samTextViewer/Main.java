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
import jline.console.history.History;
import jline.console.history.History.Entry;
import jline.console.history.MemoryHistory;
import net.sourceforge.argparse4j.inf.Namespace;
import tracks.IntervalFeature;
import tracks.Track;
import tracks.TrackCoverage;
import tracks.TrackFormat;
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
		final String genome= opts.getString("genome");
		final String fasta= opts.getString("fasta");
		String exec= opts.getString("exec");
		String config= opts.getString("config");
		exec= parseExec(exec);
		boolean debug= opts.getBoolean("debug");

		// Get configuration. Note that we don't need to assign this to a variable. 
		new Config(config);
		
		// Init console right at start so if something goes wrong the user's terminal is reset to 
		// initial defaults with the shutdown hook. This could be achieved in cleaner way probably.
		ConsoleReader console = CommandList.initConsole();
		
		messageVersion(opts.getBoolean("noFormat"));
		
		/* Set up console */
		
		Utils.checkFasta(fasta);
		
		/* Test input files exist */
		List<String> inputFileList= new ArrayList<String>();
		Utils.addSourceName(inputFileList, initFileList);
		
		/* Initialize trackSet */
		/* ------------------- */
		// This part only prepares a dummy GenomicCoords object to initialize the start position:
		// ----------------------------
		region= initRegion(region, inputFileList, fasta, genome, debug);
		
		GenomicCoords initGc= new GenomicCoords(region, null, null);
		
		List<String>initGenomeList= new ArrayList<String>();
		for(String x : inputFileList){
			initGenomeList.add(x);
		}
		initGenomeList.add(fasta);
		initGenomeList.add(genome);
		initGc.setGenome(initGenomeList);
		// ----------------------------
		// Genomic positions start here:
		final GenomicCoordsHistory gch= new GenomicCoordsHistory();
		GenomicCoords start= new GenomicCoords(initGc.toStringRegion(), initGc.getSamSeqDict(), initGc.getFastaFile());
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
			String str= String.format("\033[48;5;%sm", Config.getColor(ConfigKey.background));
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
				IntervalFeature target= new IntervalFeature(line, TrackFormat.BED);
				String reg= target.getChrom() + ":" + target.getFrom() + "-" + target.getTo();
				String gotoAndExec= ("goto " + reg + " && " + exec).trim().replaceAll("&&$", "");
				InteractiveInput itr = new InteractiveInput(console);
				itr.processInput(gotoAndExec, proc, debug);
				if (itr.getInteractiveInputExitCode() != 0){
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
			int currentExitCode= 9;
			interactiveInput.setInteractiveInputExitCode(currentExitCode);
			
			while(interactiveInput.getInteractiveInputExitCode() != 0){
				
				console.setPrompt(
						StringUtils.repeat(' ', proc.getWindowSize()) + '\r' + "[h] for help: "
						);

				cmdConcatInput= console.readLine().trim();
				if (cmdConcatInput.isEmpty()) {
					if(interactiveInput.getInteractiveInputExitCode() == 0 || 
					   interactiveInput.getInteractiveInputExitCode() == currentExitCode){
						// User only issued <ENTER>: Repeat previous command if the exit code was not an error.
						cmdConcatInput= currentCmdConcatInput;					
					} else {
						cmdConcatInput= "+0";
					}
				}
				interactiveInput.processInput(cmdConcatInput, proc, debug);
				currentCmdConcatInput= cmdConcatInput;
			}
			// *** END processing interactive input 
		}
	}

	// private static Config config= new Config(null);
	
//	public static Config getConfig(){
//		return new Config(main.);
//	}
		
	/** Return a suitable region to start. If a region is already given, do nothing.
	 * This method is a mess and should be cleaned up together with GenomicCoords class.
	 * */
	public static String initRegion(String region, List<String> inputFileList, String fasta, String genome, boolean debug) throws IOException{

		if( region != null && ! region.isEmpty() ){
			return region;
		}

		// Try to initialize from fasta
		if(fasta != null && ! fasta.trim().isEmpty()){ 
			IndexedFastaSequenceFile faSeqFile = new IndexedFastaSequenceFile(new File(fasta));
			region= faSeqFile.nextSequence().getName();
			faSeqFile.close();
			return region;
		}
		
		/* Prepare genomic coordinates to fetch. This should probably be a function in itself */
		// Create a dummy gc object just to get the sequence dict.
		GenomicCoords gc= new GenomicCoords();
		
		List<String>initGenomeList= new ArrayList<String>(inputFileList);
		
		if(genome != null && ! genome.trim().isEmpty()){
			initGenomeList.add(genome);
		}
		
		gc.setGenome(initGenomeList);
		SAMSequenceDictionary samSeqDict = gc.getSamSeqDict();
		
		System.err.print("Initializing coordinates... ");
		if(samSeqDict != null && ! samSeqDict.isEmpty()){
			region= samSeqDict.getSequence(0).getSequenceName();
			System.err.println("");
		} else {
			for(String x : initGenomeList){
				try {
					region= Utils.initRegionFromFile(x);
					System.err.println("Done from: " + x);
					break;
				} catch(Exception e){
					System.err.println("\nCould not initilize from file " + x);
					if(debug){
						e.printStackTrace();
					}
				}
			}
		}
		return region;
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
	private static History initCmdHistory(String cmdHistoryFile, String MARKER_FOR_HISTORY_CMD, boolean debug){
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
			if(debug){
				e.printStackTrace();
			}
		}
		return cmdHistory;
	}	

	/** Merge set of opened files in trackSet with the files found in the history file.
	 * */
	private static void addHistoryFiles(TrackSet trackSet, String cmdHistoryFile, String MARKER_FOR_HISTORY_FILE, boolean debug){
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
			String MARKER_FOR_HISTORY_FILE, boolean debug){
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
			if(debug){
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

		if (trackList.get(0).getGc().getFastaFile() != null && trackList.get(0).getGc().getBpPerScreenColumn() < 1.0001){
			consoleHeight -= 1; // Reference sequence
		}
		if (trackList.get(0).getGc().getSamSeqDict() != null){
			consoleHeight -= 1; // Chrom ideogram
		}
		consoleHeight -= 1; // Region info
		consoleHeight -= 1; // prompt
		consoleHeight -= trackList.size(); // Track headers

		int consensus= 0; // Additional line for possible consensus sequence in TrackCoverage
		if(trackList.get(0).getGc().getBpPerScreenColumn() < 1.0001){
			consensus= 1;
		}
		for(Track tr : trackList){
			if(tr instanceof TrackCoverage){
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

		Runtime.getRuntime().addShutdownHook(new Thread(new Runnable() {
			public void run() {
				try{
					List<String> up = Utils.checkUpdates(5000);
					int cmp= Utils.versionCompare(up.get(0), up.get(1));
					String msg= "";
					if(cmp == -1){
						msg= "NOTE: Newer version of ASCIIGenome is available: v" + up.get(1);
					}
					if( ! noFormat){
						msg= "\033[48;5;231;38;5;" + Xterm256.colorNameToXterm256("red") + "m" +  msg + "\033[0m";
					}
					System.err.println(msg);
				} catch(Exception e){
					// e.printStackTrace();
				}
	        }
	    }, "Shutdown-thread"));		
	}
	
}
