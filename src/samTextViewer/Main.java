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
import java.util.List;
import java.util.ListIterator;

import com.google.common.base.Joiner;
import com.itextpdf.text.DocumentException;

import commandHelp.CommandList;
import exceptions.BamIndexNotFoundException;
import exceptions.InvalidCommandLineException;
import exceptions.InvalidGenomicCoordsException;
import exceptions.InvalidRecordException;
import net.sourceforge.argparse4j.inf.Namespace;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import jline.console.ConsoleReader;
import jline.console.history.History;
import jline.console.history.History.Entry;
import jline.console.history.MemoryHistory;
import tracks.IntervalFeature;
import tracks.Track;
import tracks.TrackCoverage;
import tracks.TrackFormat;
import tracks.TrackReads;
import tracks.TrackSet;

/**
 * @author berald01
 *
 */
public class Main {
	
	public static void main(String[] args) throws IOException, InvalidGenomicCoordsException, InvalidCommandLineException, InvalidRecordException, BamIndexNotFoundException, ClassNotFoundException, SQLException, DocumentException {

		final String cmdHistoryFile= System.getProperty("user.home") + File.separator + ".asciigenome_history";  
		
		/* Start parsing arguments * 
		 * *** If you change something here change also in console input ***/
		Namespace opts= ArgParse.argParse(args);
		
		List<String> initFileList= opts.getList("input");
		String region= opts.getString("region");
		final String genome= opts.getString("genome");
		final String fasta= opts.getString("fasta");
		String exec= opts.getString("exec");
		exec= parseExec(exec);
		// Init console right at start so if something goes wrong the user's terminal is reset to 
		// initial defaults with the shutdown hook. This could be achieved in cleaner way probably.
		ConsoleReader console = CommandList.initConsole();
		
		/* Set up console */
		
		Utils.checkFasta(fasta);
		
		/* Test input files exist */
		List<String> inputFileList= new ArrayList<String>();
		try{
			Utils.addSourceName(inputFileList, initFileList);
		} catch(InvalidCommandLineException e){
			//
		}	

		/* Initialize trackSet */
		/* ------------------- */
		// This part only prepares a dummy GenomicCoords object to initialize the start position:
		// ----------------------------
		region= initRegion(region, inputFileList, fasta, genome);
		
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
		gch.add(new GenomicCoords(initGc.toStringRegion(), initGc.getSamSeqDict(), initGc.getFastaFile()));

		final TrackSet trackSet= new TrackSet(inputFileList, gch.current());
		
		setDefaultTrackHeights(console.getTerminal().getHeight(), trackSet.getTrackList());
		
		final TrackProcessor proc= new TrackProcessor(trackSet, gch);
		
		proc.setNoFormat(opts.getBoolean("noFormat"));
		
		// Put here the previous command so that it is re-issued if no input is given
		// You have to initialize this var outside the while loop that processes input files.
		String currentCmdConcatInput= ""; 

		if(!proc.isNoFormat()){
			System.out.print("\033[48;5;231m");
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
				itr.processInput(gotoAndExec, proc);
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
			itr.processInput(exec, proc);
			if(opts.getBoolean("nonInteractive")){
				System.out.print("\033[0m");
				System.exit(0);
			}
		}

		/* Set up done, start processing */
		/* ============================= */
		History cmdHistory= initCmdHistory(cmdHistoryFile);
		console.setHistory(cmdHistory);
		writeHistory(console.getHistory(), cmdHistoryFile, 500);
		
		while(true){  
			// keep going until quit or if no interactive input set
			// *** START processing interactive input
			String cmdConcatInput= ""; // String like "zi && -F 16 && mapq 10"
			InteractiveInput interactiveInput= new InteractiveInput(console);
			interactiveInput.setInteractiveInputExitCode(9);
			
			while(interactiveInput.getInteractiveInputExitCode() != 0){
				
				console.setPrompt("[h] for help: ");
				cmdConcatInput= console.readLine().trim();
				if (cmdConcatInput.isEmpty()){
					// Repeat previous command
					cmdConcatInput= currentCmdConcatInput;
				}
				interactiveInput.processInput(cmdConcatInput, proc);
				currentCmdConcatInput= cmdConcatInput;
				
			}
			// *** END processing interactive input 
		}
	}

	/** Return a suitable region to start. If a region is already given, do nothing.
	 * This method is a mess and should be cleaned up together with GenomicCoords class.
	 * */
	public static String initRegion(String region, List<String> inputFileList, String fasta, String genome) throws IOException{

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
	private static History initCmdHistory(String cmdHistoryFile){
		History cmdHistory= new MemoryHistory();
		try{
			BufferedReader br= new BufferedReader(new FileReader(new File(cmdHistoryFile)));
			String line;
			while((line = br.readLine()) != null){
				cmdHistory.add(line);
			}
			br.close();
		} catch(IOException e){
			//
		}
		return cmdHistory;
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
			trackList.get(i).setyMaxLines(trackHeights.get(i));
		}
	}
	
	/**Write the history of commands to given file.
	 * */
	private static void writeHistory(final History cmdHistory, final String cmdHistoryFile, final int nmax){
		
		Runtime.getRuntime().addShutdownHook(new Thread(new Runnable() {
			public void run() {
	            
				ListIterator<Entry> iter = cmdHistory.entries();
				List<String>lastHist= new ArrayList<String>();
				int i= nmax;
				while(iter.hasNext()){
					if(i == 0){
						break;
					}
					i--;
					lastHist.add(iter.next().value().toString());
				}
				
				try {
					BufferedWriter wr= new BufferedWriter(new FileWriter(new File(cmdHistoryFile)));
					for(String cmd : lastHist){
						wr.write(cmd + "\n");
					}
					wr.close();
				} catch (IOException e) {
					System.err.println("Unable to write history to " + cmdHistoryFile);
				}
	        }
	    }, "Shutdown-thread"));

	}

}
