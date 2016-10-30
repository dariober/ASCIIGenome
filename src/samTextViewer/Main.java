package samTextViewer;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.List;

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
		
		if((region == null || region.isEmpty()) && fasta != null){ // Try to initilize from fasta
			IndexedFastaSequenceFile faSeqFile = new IndexedFastaSequenceFile(new File(fasta));
			region= faSeqFile.nextSequence().getName();
			faSeqFile.close();
		}

		SAMSequenceDictionary samSeqDict = GenomicCoords.getSamSeqDictFromAnyFile(inputFileList, fasta, genome);
		
		/* Prepare genomic coordinates to fetch. This should probably be a function in itself */
		if(region.isEmpty()){
			System.err.print("Initializing coordinates... ");
			if(!samSeqDict.isEmpty()){
				region= samSeqDict.getSequence(0).getSequenceName();
				System.err.println("");
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
		
		/* Initialize trackSet */
		/* ------------------- */
		GenomicCoordsHistory gch= new GenomicCoordsHistory();
		gch.add(new GenomicCoords(region, samSeqDict, fasta));

		TrackProcessor proc= new TrackProcessor(new TrackSet(inputFileList, gch.current()), gch);

		if(proc.getGenomicCoordsHistory().current().getFastaFile() != null){
			TrackSeqRegex re= new TrackSeqRegex(proc.getGenomicCoordsHistory().current());
			proc.getTrackSet().add(re, "regex_seq_matches");;
		}
		
		proc.setNoFormat(opts.getBoolean("noFormat"));
		
		// Put here the previous command so that it is re-issued if no imput is given
		// You have to initialize this var outside the while loop that processes input files.
		String currentCmdConcatInput= ""; 

		if(!proc.isNoFormat()){
			System.out.print("\033[48;5;231m");
		}

		// Batch processing file of regions
		final String batchFile= opts.getString("batchFile");
		if(! batchFile.isEmpty()){

			if(! new File(batchFile).exists()){
				System.err.print("\033[0m");
				System.err.println("File " + batchFile + " does not exist.");
				System.exit(1);
			}

			console.clearScreen();
			console.flush();

			BufferedReader br= new BufferedReader(new FileReader(new File(batchFile)));
			String line = null;  
			while ((line = br.readLine()) != null){

				IntervalFeature target= new IntervalFeature(line, TrackFormat.BED);
				String reg= target.getChrom() + ":" + target.getFrom() + "-" + target.getTo();
				String gotoAndExec= ("goto " + reg + " && " + exec).trim().replaceAll("&&$", "");
				InteractiveInput itr = new InteractiveInput();
				proc= itr.processInput(gotoAndExec, proc);
				if(itr.getInteractiveInputExitCode() == 0){
					proc.iterateTracks();
					Utils.printer("\n", proc.getSnapshotFile());
				}
			}
			br.close();
			System.exit(0);
		}

		// See if we need to process the exec arg before going to interactive mode. 
		// Also if we are in non-interactive mode, we process the track set now and later exit 
		if(!exec.isEmpty() || opts.getBoolean("nonInteractive")){

			console.clearScreen();
			console.flush();
			
			InteractiveInput itr = new InteractiveInput();
			proc= itr.processInput(exec, proc);
			if(itr.getInteractiveInputExitCode() == 0){
				proc.iterateTracks();
				Utils.printer("\n", proc.getSnapshotFile());
				proc.setSnapshotFile(null);
			}
			
			if(opts.getBoolean("nonInteractive")){
				System.out.print("\033[0m");
				System.exit(0);
			}
		}

		/* Set up done, start processing */
		/* ============================= */
		List<String> cmdHistory= new ArrayList<String>();
		while(true){ 
			
			console.clearScreen();
			console.flush();
			
			proc.iterateTracks();
			proc.setSnapshotFile(null); // This is to prevent taking screenshots one after another. It's a hack and should be changed

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
				proc= interactiveInput.processInput(cmdConcatInput, proc);
				currentCmdConcatInput= cmdConcatInput;
			}
			// *** END processing interactive input 
			
		} // End while loop keep going until quit or if no interactive input set
	}
	
}
