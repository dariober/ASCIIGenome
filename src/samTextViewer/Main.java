package samTextViewer;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.regex.Pattern;
import java.util.regex.PatternSyntaxException;

import org.apache.commons.lang3.text.StrTokenizer;

import exceptions.InvalidCommandLineException;
import exceptions.InvalidGenomicCoordsException;
import net.sourceforge.argparse4j.inf.Namespace;
import filter.FlagToFilter;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.filter.MappingQualityFilter;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import jline.console.ConsoleReader;
import jline.console.completer.StringsCompleter;
import tracks.Track;
import tracks.TrackCoverage;
import tracks.TrackFormat;
import tracks.TrackIntervalFeature;
import tracks.TrackMethylation;
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
		
		List<String> inputFileList= opts.getList("insam");
		String region= opts.getString("region");
		String genome= opts.getString("genome");
		//int windowSize= opts.getInt("windowSize");
		String fasta= opts.getString("fasta");
		boolean rpm= opts.getBoolean("rpm");
		int maxLines= opts.getInt("maxLines");
		// int trackHeight= opts.getInt("maxDepthLines");
		int maxMethylLines= opts.getInt("maxMethylLines");
		final int maxReadsStack= opts.getInt("maxReadsStack");
		int f_incl= opts.getInt("f");
		int F_excl= opts.getInt("F");
		int mapq= opts.getInt("mapq");
		boolean bs= opts.getBoolean("BSseq");
		boolean noFormat= opts.getBoolean("noFormat");
		boolean nonInteractive= opts.getBoolean("nonInteractive");
		boolean withReadName= false; // FIXME: Add to parser?
		
		if((F_excl & 4) != 4){ // Always filter out read unmapped
			F_excl += 4;
		}
			
		if(fasta == null && bs == true){
			System.err.println("Warning: Fasta reference not provided. Bisulfite mode will be disabled");
			bs= false;
		}
		
		int windowSize= 160;
		try{
			int terminalWidth = jline.TerminalFactory.get().getWidth();
			windowSize= (int) (terminalWidth * 0.999); 
		} catch(Exception e){
			e.printStackTrace();
		}
		
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

		/* Initialize console */
		ConsoleReader console = new ConsoleReader();
		// console.addCompleter(new FileNameCompleter());
		for(String x : inputFileList){
			console.addCompleter(new StringsCompleter(new File(x).getName()));
		}
		for(String x : "next next_start goto find_next find_all showGenome addTracks visible trackHeight ylim dataCol print printFull rNameOn rNameOff history".split(" ")){
			// Add options. Really you should use a dict for this.
			if(x.length() > 2){
				console.addCompleter(new StringsCompleter(x));
			}
		}
		
		String currentCmd = null; // Used to store the current interactive command and repeat it if no new cmd is given. 
		
		boolean printIntervalFeatures= false;
		boolean printIntervalFeaturesFull= false;
		TrackSet trackSet= new TrackSet();

		/* Initialize GC profile */
		if(fasta != null){
			TrackWiggles cgWiggle= gch.current().getGCProfile();
			trackSet.addOrReplace(cgWiggle);
		}
		
		int idForTrack= 0;
		while(true){ // Each loop processes the user's input files.

			/* Prepare filters */
			List<SamRecordFilter> filters= FlagToFilter.flagToFilterList(f_incl, F_excl);
			filters.add(new MappingQualityFilter(mapq));
			
			for(int i= 0; i < inputFileList.size(); i++){ 
				/* Iterate through each input file */
				String inputFileName= inputFileList.get(i);
				
				if(Utils.getFileTypeFromName(inputFileName).equals(TrackFormat.BAM)){
				
					/* Coverage and methylation track */
					String coverageTrackId= new File(inputFileName).getName() + "#" + (idForTrack+1);
					idForTrack++;
					if(!trackSet.getTrackSet().containsKey(coverageTrackId)){
						TrackCoverage trackCoverage= new TrackCoverage(inputFileName, gch.current(), filters, bs);
						trackCoverage.setFileTag(coverageTrackId);
						trackSet.addOrReplace(trackCoverage);
					}
					TrackCoverage trackCoverage= (TrackCoverage) trackSet.getTrackSet().get(coverageTrackId);
					trackCoverage.setGc(gch.current());
					trackCoverage.setFilters(filters);
					trackCoverage.setRpm(rpm);
					trackCoverage.update();
					trackCoverage.printToScreen();				
					
					if(bs && trackCoverage.getScreenLocusInfoList().size() > 0){
						if(maxMethylLines < 0){
							maxMethylLines= 0;
						}
						coverageTrackId= new File(inputFileName).getName() + "#" + (idForTrack+1);
						idForTrack++;
						if(!trackSet.getTrackSet().containsKey(coverageTrackId)){
							TrackMethylation trackMethylation= new TrackMethylation(inputFileName, trackCoverage.getScreenLocusInfoList());
							trackMethylation.setFileTag(coverageTrackId);
							trackSet.addOrReplace(trackMethylation);
						}
						TrackMethylation trackMethylation= (TrackMethylation) trackSet.getTrackSet().get(coverageTrackId);
						trackMethylation.setScreenLocusInfoList(trackCoverage.getScreenLocusInfoList());
						trackMethylation.setyMaxLines(maxMethylLines);
					}
										
					/* Reads */
					String trackId= new File(inputFileName).getName() + "#" + (idForTrack+1);
					idForTrack++;
					if(!trackSet.getTrackSet().containsKey(trackId)){
						TrackReads trackReads= new TrackReads(inputFileName, gch.current(), filters, maxReadsStack);
						trackReads.setFileTag(trackId);
						trackSet.addOrReplace(trackReads);
						trackReads.setFilename(inputFileName);
						trackReads.setFileTag(trackId);
					}
					TrackReads trackReads= (TrackReads) trackSet.getTrackSet().get(trackId);
					trackReads.setGc(gch.current());
					trackReads.setFilters(filters);
					trackReads.setyMaxLines(maxLines);
					trackReads.setBs(bs);
					trackReads.setWithReadName(withReadName);
					trackReads.update();
				} // End processing bam file
				
				/* Annotatation */
				if(Utils.getFileTypeFromName(inputFileName).equals(TrackFormat.BED) 
						|| Utils.getFileTypeFromName(inputFileName).equals(TrackFormat.GFF)){
					String trackId= new File(inputFileName).getName() + "#" + (idForTrack+1);
					idForTrack++;
					if(!trackSet.getTrackSet().containsKey(trackId)){
						TrackIntervalFeature tif= new TrackIntervalFeature(inputFileName, gch.current());
						tif.setFileTag(trackId);
						trackSet.addOrReplace(tif);
					}
					TrackIntervalFeature tif= (TrackIntervalFeature) trackSet.getTrackSet().get(trackId);
					tif.setGc(gch.current());
					//tif.setyMaxLines(trackHeight);
					tif.update();
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
						trackSet.addOrReplace(tw);
					}
					TrackWiggles tw= (TrackWiggles) trackSet.getTrackSet().get(trackId);
					tw.setGc(gch.current());
					tw.update();
					tw.printToScreen();
				}
			} // End loop through files 

			/* Print tracks */
			/* ************ */
			console.clearScreen();
			console.flush();
			
			if(gch.current().getChromIdeogram(20) != null){
				System.out.println(gch.current().getChromIdeogram(20));
			}			
			for(Track tr : trackSet.getTrackSet().values()){
				if(tr.getFileTag() == gch.current().getGcProfileFileTag()){
					continue;
				}
				tr.setNoFormat(noFormat);
				if(tr.isNoFormat()){
					System.out.print(tr.getTitle());
				} else {
					System.out.print("\033[0;34m" + tr.getTitle() + "\033[0m");
				}
				if(tr.getyMaxLines() > 0){
					System.out.println(tr.printToScreen());
				}
				if(printIntervalFeatures){
					String printable= tr.printFeatures(windowSize);
					System.out.print(printable);
				}
				if(printIntervalFeaturesFull){
					String printable= tr.printFeatures(Integer.MAX_VALUE);
					System.out.print(printable);
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
				if(yMaxLines > 0){
					TrackWiggles tw= gch.current().getGCProfile();				
					tw.setyMaxLines(yMaxLines);
					tw.setYLimitMin(yLimitMin);
					tw.setYLimitMax(yLimitMax);
					System.out.print(tw.getTitle());
					String gcPrintable= tw.printToScreen();
					if(!noFormat){
						gcPrintable= "\033[0;33m" + gcPrintable + "\033[0m";
					}
					System.out.print(gcPrintable + "\n");
				}
			}
			System.out.print(gch.current().printableRefSeq(noFormat));
			String ruler= gch.current().printableRuler(10);
			System.out.println(ruler.substring(0, ruler.length() <= windowSize ? ruler.length() : windowSize));

			String footer= gch.current().toString() + "; " + Math.rint(gch.current().getBpPerScreenColumn() * 10d)/10d + " bp/char; " 
					+ "Filters: -q " + mapq  + " -f " + f_incl + " -F " + F_excl
					+ "; " + getMemoryStat();
			if(!noFormat){
				System.out.println("\033[0;34m" + footer + "\033[0m; ");
			} else {
				System.out.println(footer);
			}
			/* Interactive input */
			/* ================= */
			if(!nonInteractive){
				break;
			}
			
			String cmdInput= null;
			while(cmdInput == null){ // Keep asking for input until you get something valid
				console.setPrompt("[h] for help: ");
				cmdInput = console.readLine().trim();
				
				if (cmdInput.trim().isEmpty()){
					// Repeat previous command
					cmdInput= currentCmd;
				}
				
				if(cmdInput == null || cmdInput.equals("h")){
String inline= "\n    N a v i g a t i o n\n\n"
+ "f / b \n"
+ "      Small step forward/backward 1/10 window\n"
+ "ff / bb\n"
+ "zi / zo [x]\n"
+ "      Zoom in / zoom out x times (default x= 1). Each zoom halves or doubles the window size\n"
+ "      Large step forward/backward 1/2 window\n"
+ "goto chrom:from-to\n"
+ "      Go to given region. E.g. \"goto chr1:1-1000\" or chr1:10 or chr1. goto keyword can be replaced with ':' (like goto in vim)\n"
+ "<from> [to]\n"
+ "      Go to position <from> or to region \"from to\" on current chromosome. E.g. 10 or \"10 1000\" or \"10-1000\"\n" 
+ "+/-<int>[k,m]\n"
+ "      Move forward/backward by <int> bases. Suffixes k (kilo) and M (mega) allowed. E.g. -2m or +10k\n"
+ "p / n\n"
+ "      Go to previous/next visited position\n"
+ "next / next_start [trackId]\n"
+ "      Move to the next feature in trackId on *current* chromosome\n"
+ "      'next' centers the window on the found feature while 'next_start' sets the window at the start of the feature.\n"

+ "\n    F i n d  \n\n"
+ "find_next <regex> [trackId]\n"
+ "      Find the next record in trackId matching regex. Use single quotes for strings containing spaces.\n"
+ "      For case insensitive matching prepend (?i) to regex. E.g. \"next '(?i).*actb.*' myTrack#1\"\n"
+ "find_all <regex> [trackId]\n"
+ "      Find all matches on chromosome. The search stops at the first chromosome returning hits\n"
+ "      starting with the current one. Useful to get all gtf records of a gene\n"

+ "\n    D i s p l a y  \n\n"
+ "visible [show regex] [hide regex] [track regex]\n"
+ "      In annotation tracks, only include rows captured by [show regex] and exclude [hide regex].\n"
+ "      Apply to tracks captured by [track regex]. With no optional arguments reset to default: \"'.*' '^$' '.*'\"\n"
+ "      Use '.*' to match everything and '^$' to hide nothing. E.g. \"visible .*exon.* .*CDS.* .*gtf#.*\"\n"       
+ "trackHeight <int> [track regex]\n"
+ "      Set track height to int lines for all tracks captured by regex. Default regex: '.*'\n"
+ "ylim <min> <max> [track regex]\n"
+ "      Set limits of y axis for all track IDs captured by regex. Use na to autoscale to min and/or max.\n"
+ "      E.g. ylim 0 na. If regex is omitted all tracks will be captured. Default: \"ylim na na .*\"\n"
+ "dataCol <idx> [regex]\n"
+ "      Select data column for all bedgraph tracks captured by regex. <idx>: 1-based column index.\n"
+ "print / printFull\n"
+ "      Turn on/off the printing of bed/gtf features.\n"
+ "      print clip lines to fit the screen, printFull will wrap the long lines\n"
+ "showGenome\n"
+ "      Print the genome file\n"
+ "addTracks [file or url]...\n"
+ "      Add tracks\n" 
+ "history\n"
+ "      Show visited positions\n";
					System.out.println(inline);
					System.out.println("    A l i g n m e n t s\n");
					System.out.println(ArgParse.getDocstrings());
					System.out.println("q      Quit");
					System.out.println("See also http://github.com/dariober/Java-cafe/tree/master/SamTextViewer");
					cmdInput= null;
					continue;
				} 
				/* Parse args */
				/* ---------- */
				try{
					if(cmdInput.equals("q")){
						System.exit(0);
					}
					
					if(cmdInput.equals("f") 
						|| cmdInput.equals("b")
						|| cmdInput.equals("ff") 
						|| cmdInput.equals("bb")
						|| cmdInput.matches("^\\d+.*")
						|| cmdInput.matches("^\\-\\d+.*") 
						|| cmdInput.matches("^\\+\\d+.*")){ // No cmd line args either f/b ops or ints
						String newRegion= Utils.parseConsoleInput(cmdInput, gch.current()).trim();
						GenomicCoords newGc= new GenomicCoords(newRegion, samSeqDict, windowSize, fasta);
						gch.add(newGc);	
					} else if(cmdInput.startsWith("goto") || cmdInput.startsWith(":")){
						String reg= cmdInput.replaceFirst("goto|:", "").trim();
						gch.add(new GenomicCoords(reg, samSeqDict, windowSize, fasta));
					} else if(cmdInput.startsWith("dataCol ")){
						
						StrTokenizer str= new StrTokenizer(cmdInput);
						str.setQuoteChar('\'');
						List<String> tokens= str.getTokenList();
						String trackIdRegex= ".*";
						if(tokens.size() >= 3){
							trackIdRegex= tokens.get(2);
							try{
								Pattern.compile(trackIdRegex);
							} catch(PatternSyntaxException e){
						    	System.err.println("Invalid regex in: " + cmdInput);
						    	System.err.println(e.getDescription());
								cmdInput= null;
								continue;
							}
						}
						trackSet.selectDataColumnForBedgraph(Integer.parseInt(tokens.get(1)), trackIdRegex);
					
					} else if(cmdInput.startsWith("ylim ")){
					
						try{
							trackSet.setTrackYlimitsForRegex(cmdInput);
						} catch(InvalidCommandLineException e){
							cmdInput= null;
							continue;
						} catch(PatternSyntaxException e) {
							cmdInput= null;
				        	continue;
						}
					} else if(cmdInput.startsWith("trackHeight ")){
						try{
							trackSet.setTrackHeightForRegex(cmdInput);
						} catch(InvalidCommandLineException e){
							cmdInput= null;
							continue;
						} catch(PatternSyntaxException e) {
							cmdInput= null;
				        	continue;
						}
					} else if(cmdInput.startsWith("addTracks ")){
						StrTokenizer str= new StrTokenizer(cmdInput);
						str.setQuoteChar('\'');
						List<String> newFileNames= str.getTokenList();
						newFileNames.remove(0);
						Utils.addTrack(inputFileList, newFileNames);
			
						if(gch.current().getSamSeqDict().size() == 0){
							samSeqDict = GenomicCoords.getSamSeqDictFromAnyFile(inputFileList, null, null);
							GenomicCoords gc= gch.current();
							gc.setSamSeqDict(samSeqDict);
						}
					} else if (cmdInput.equals("p")) {
						gch.previous();
					} else if (cmdInput.equals("n")) {
						gch.next();
					} else if(cmdInput.startsWith("zo")){
						int nz= Utils.parseZoom(cmdInput, 1);
						GenomicCoords gc = (GenomicCoords)gch.current().clone();
						for(int i= 0; i < nz; i++){
							gc.zoomOut();
						}
						gch.add(gc);
					} else if(cmdInput.startsWith("zi")){
						int nz= Utils.parseZoom(cmdInput, 1);
						GenomicCoords gc = (GenomicCoords)gch.current().clone();
						for(int i= 0; i < nz; i++){
							gc.zoomIn();
						}
						gch.add(gc);
					} else if(cmdInput.equals("history")){
						for(GenomicCoords x : gch.getHistory()){
							System.out.println(x);
						}
						cmdInput= null;
					} else if(cmdInput.toLowerCase().equals("print")){
						printIntervalFeatures= (printIntervalFeatures) ? false : true;
						printIntervalFeaturesFull= false;
						System.err.println("Print interval features: " + printIntervalFeatures);
					} else if(cmdInput.toLowerCase().equals("printfull")){
						printIntervalFeaturesFull= (printIntervalFeaturesFull) ? false : true;
						printIntervalFeatures= false;
						System.err.println("Print full interval features: " + printIntervalFeaturesFull);
					} else if(cmdInput.toLowerCase().equals("rnameon")){
						withReadName= true;
					} else if(cmdInput.toLowerCase().equals("rnameoff")) {
						withReadName= false;
					} else if(cmdInput.startsWith("next_start ") || cmdInput.equals("next_start")){
						GenomicCoords gc= (GenomicCoords)gch.current().clone();
						gch.add(trackSet.goToNextFeatureOnFile(cmdInput.replace("next_start", "").trim(), gc, -1.0));
					} else if(cmdInput.startsWith("next ") || cmdInput.equals("next")){
							GenomicCoords gc= (GenomicCoords)gch.current().clone();
							gch.add(trackSet.goToNextFeatureOnFile(cmdInput.replace("next", "").trim(), gc, 5.0));
					} else if(cmdInput.startsWith("find_next ") || cmdInput.startsWith("find_all ")) {  
						StrTokenizer str= new StrTokenizer(cmdInput);
						str.setQuoteChar('\'');
						List<String> tokens= str.getTokenList();
						if(tokens.size() < 2){
							System.err.println("Error in find* subcommand. Expected at least 2 args got: " + cmdInput);
							cmdInput= null;
							continue;
						}
						if(tokens.size() == 2){
							tokens.add("");
						}
						GenomicCoords gc= (GenomicCoords)gch.current().clone();
						boolean all= cmdInput.startsWith("find_all ") ? true : false; 
						gch.add(trackSet.findNextRegexOnTrack(tokens.get(1), tokens.get(2), gc, all));
						
					} else if(cmdInput.startsWith("visible ") || cmdInput.equals("visible")){
						try{
							trackSet.setVisibilityForTrackIntervalFeature(cmdInput);
						} catch (PatternSyntaxException e){
							System.err.println("Invalid pattern in " + cmdInput);
							cmdInput= null;
							continue;							
						}
					} else if(cmdInput.equals("showGenome")) {
						System.out.println(Utils.printSamSeqDict(gch.current().getSamSeqDict(), 30));
						cmdInput= null;
						continue;
					// Command line options from Argparse
					} else { 
						List<String> clArgs= Arrays.asList(cmdInput.split("\\s+"));
						if(clArgs.indexOf("-f") != -1){
							int i= clArgs.indexOf("-f") + 1;
							f_incl= Integer.parseInt(clArgs.get(i));
						} else if(clArgs.indexOf("-F") != -1){
							int i= clArgs.indexOf("-F") + 1;
							F_excl= Integer.parseInt(clArgs.get(i));
							if((F_excl & 4) != 4){ // Always filter out read unmapped
								F_excl += 4;
							}
						} else if(clArgs.indexOf("-q") != -1){
							int i= clArgs.indexOf("-q") + 1;
							mapq= Integer.parseInt(clArgs.get(i));
						} else if(clArgs.indexOf("-m") != -1){
							int i= clArgs.indexOf("-m") + 1;
							maxLines= Integer.parseInt(clArgs.get(i));
						} else if(clArgs.indexOf("-ml") != -1){
							int i= clArgs.indexOf("-ml") + 1;
							maxMethylLines= Integer.parseInt(clArgs.get(i));
						} else if(clArgs.indexOf("-rpm") != -1){
							rpm= (rpm) ? false : true; // Invert rpm set.
						// END OF CMD LINE ARGS
						} else {
							System.err.println("Unrecognized argument: " + cmdInput);
							cmdInput= null;
						}
					} // END Command line options from Argparse
				} catch(Exception e){
					System.err.println("\nError processing input: " + cmdInput + "\n");
					e.printStackTrace(); // Print trace for debugging
					System.err.println("");
					cmdInput= null;
				}
				currentCmd= cmdInput;
			} // END while loop to get cmLine args
			idForTrack= 0;
		} // End while loop keep going until quit or if no interactive input set
	}
}
