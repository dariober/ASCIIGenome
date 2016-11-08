package samTextViewer;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.List;

import com.google.common.base.Joiner;

import commandHelp.CommandList;
import exceptions.BamIndexNotFoundException;
import exceptions.InvalidCommandLineException;
import exceptions.InvalidGenomicCoordsException;
import exceptions.InvalidRecordException;
import net.sourceforge.argparse4j.inf.Namespace;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import jline.console.ConsoleReader;
import tracks.IntervalFeature;
import tracks.TrackFormat;
import tracks.TrackSeqRegex;
import tracks.TrackSet;

/**
 * @author berald01
 *
 */
public class Main {
	
	public static void main(String[] args) throws IOException, InvalidGenomicCoordsException, InvalidCommandLineException, InvalidRecordException, BamIndexNotFoundException, ClassNotFoundException, SQLException {

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
		final GenomicCoordsHistory gch= new GenomicCoordsHistory();
		region= initRegion(region, inputFileList, fasta, genome);
		gch.add(new GenomicCoords(region, null, null));
		List<String>initGenomeList= new ArrayList<String>();
		for(String x : inputFileList){
			initGenomeList.add(x);
		}
		initGenomeList.add(fasta);
		initGenomeList.add(genome);
		gch.current().setGenome(initGenomeList);
		final TrackSet trackSet= new TrackSet(inputFileList, gch.current());
		final TrackProcessor proc= new TrackProcessor(trackSet, gch);
		
//		if(proc.getGenomicCoordsHistory().current().getFastaFile() != null){
//			TrackSeqRegex re= new TrackSeqRegex(proc.getGenomicCoordsHistory().current());
//			proc.getTrackSet().addTrack(re, "regex_seq_matches");
//		}

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
				InteractiveInput itr = new InteractiveInput();
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

			InteractiveInput itr = new InteractiveInput();
			itr.processInput(exec, proc);
			if(opts.getBoolean("nonInteractive")){
				System.out.print("\033[0m");
				System.exit(0);
			}
		}

		/* Set up done, start processing */
		/* ============================= */
		List<String> cmdHistory= new ArrayList<String>();
		while(true){ 
			
			// *** START processing interactive input
			String cmdConcatInput= ""; // String like "zi && -F 16 && mapq 10"
			InteractiveInput interactiveInput= new InteractiveInput();
			interactiveInput.setInteractiveInputExitCode(9);
			
			while(interactiveInput.getInteractiveInputExitCode() != 0){
				
				console.setPrompt("[h] for help: ");
				cmdConcatInput= console.readLine().trim();
				if (cmdConcatInput.isEmpty()){
					// Repeat previous command
					cmdConcatInput= currentCmdConcatInput;
				}
				cmdHistory.add(cmdConcatInput);
				interactiveInput.setCmdHistory(cmdHistory);
				interactiveInput.processInput(cmdConcatInput, proc);
				currentCmdConcatInput= cmdConcatInput;
			}
			// *** END processing interactive input 
			
		} // End while loop keep going until quit or if no interactive input set
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

}
