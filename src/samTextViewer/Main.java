package samTextViewer;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Pattern;
import java.util.regex.PatternSyntaxException;

import org.apache.commons.lang3.text.StrMatcher;
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
import jline.console.completer.AggregateCompleter;
import jline.console.completer.StringsCompleter;
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
		String genome= opts.getString("genome");
		String fasta= opts.getString("fasta");
		boolean rpm= false; // opts.getBoolean("rpm");
		final int maxReadsStack= opts.getInt("maxReadsStack");
		int f_incl= 0; // opts.getInt("f");
		int F_excl= 0; // opts.getInt("F");
		int mapq= 0; // opts.getInt("mapq");
		boolean noFormat= opts.getBoolean("noFormat");
		boolean nonInteractive= opts.getBoolean("nonInteractive");
		boolean withReadName= false; // FIXME: Add to parser?
		
		if((F_excl & 4) != 4){ // Always filter out read unmapped
			F_excl += 4;
		}
			
		//if(fasta == null && bs == true){
		//	System.err.println("Warning: Fasta reference not provided. Bisulfite mode will be disabled");
		//	bs= false;
		//}
		
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

		String currentCmd = null; // Used to store the current interactive command and repeat it if no new cmd is given. 
		
		String printIntervalFeaturesRegex= "x^"; // false;
		String printIntervalFeaturesFullRegex= "x^";
		String rpmRegex= "x^";
		
		TrackSet trackSet= new TrackSet();

		/* Initialize GC profile */
		if(fasta != null){
			TrackWiggles cgWiggle= gch.current().getGCProfile();
			trackSet.getTrackSet().put(cgWiggle.getFileTag(), cgWiggle);
		}
		
		String seqRegex= null;
		int idForTrack= 0;
		String snapshotFile= null;
		boolean snapshotStripAnsi= true;
		ConsoleReader console = InlineHelp.initConsole();
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
						TrackCoverage trackCoverage= new TrackCoverage(inputFileName, gch.current(), filters, false);
						trackCoverage.setFileTag(coverageTrackId);
						trackSet.getTrackSet().put(trackCoverage.getFileTag(), trackCoverage);
					}
					TrackCoverage trackCoverage= (TrackCoverage) trackSet.getTrackSet().get(coverageTrackId);
					trackCoverage.setGc(gch.current());
					trackCoverage.setFilters(filters);
					if(Pattern.compile(rpmRegex).matcher(trackCoverage.getFileTag()).find()){
						rpm= (trackCoverage.getRpm()) ? false : true; // Invert rpm set.
						System.err.println("Setting RPM to " + rpm + " for " + trackCoverage.getFileTag());
						trackCoverage.setRpm(rpm);
					}
					trackCoverage.update();
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
					String trackId= new File(inputFileName).getName() + "#" + (idForTrack+1);
					idForTrack++;
					if(!trackSet.getTrackSet().containsKey(trackId)){
						TrackReads trackReads= new TrackReads(inputFileName, gch.current(), filters, maxReadsStack);
						trackReads.setFileTag(trackId);
						trackSet.getTrackSet().put(trackReads.getFileTag(), trackReads);
						trackReads.setFilename(inputFileName);
						trackReads.setFileTag(trackId);
					}
					TrackReads trackReads= (TrackReads) trackSet.getTrackSet().get(trackId);
					trackReads.setGc(gch.current());
					trackReads.setFilters(filters);
					trackReads.setWithReadName(withReadName);
					trackReads.update();
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
					//tif.setyMaxLines(trackHeight);
					try {
						tif.update();
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
					tw.update();
					tw.printToScreen();
				}
			} // End loop through files 
			rpmRegex= "x^"; // Reset to match nothing so that nothing changes until rpm opt is called again.  
			
			/* Print tracks */
			/* ************ */
			console.clearScreen();
			console.flush();
			
			if(gch.current().getChromIdeogram(20) != null){
				Utils.printer(gch.current().getChromIdeogram(20) + "\n", snapshotFile, snapshotStripAnsi);
			}			
			for(Track tr : trackSet.getTrackSet().values()){
				if(tr.getFileTag() == gch.current().getGcProfileFileTag()){
					continue;
				}
				tr.setNoFormat(noFormat);
				Utils.printer(tr.getTitle(), snapshotFile, snapshotStripAnsi);
				if(tr.getyMaxLines() > 0){
					Utils.printer(tr.printToScreen() + "\n", snapshotFile, snapshotStripAnsi);
				}

				// Print features
				if(Pattern.compile(printIntervalFeaturesRegex).matcher(tr.getFileTag()).find()){
					String printable= tr.printFeatures(windowSize);
					Utils.printer(printable, snapshotFile, snapshotStripAnsi);
				}
				if(Pattern.compile(printIntervalFeaturesFullRegex).matcher(tr.getFileTag()).find()){
					String printable= tr.printFeatures(Integer.MAX_VALUE);
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
					Utils.printer(tw.getTitle(), snapshotFile, snapshotStripAnsi);
					String gcPrintable= tw.printToScreen();
					Utils.printer(gcPrintable + "\n", snapshotFile, snapshotStripAnsi);
				}
			}
			// Track for matching regex
			TrackIntervalFeature seqRegexTrack = gch.current().findRegex(seqRegex);
			seqRegexTrack.setNoFormat(noFormat);
			if(trackSet.getRegexForTrackHeight().matcher(seqRegexTrack.getFileTag()).find()){
				if(trackSet.getTrackHeightForRegex() < 0){
					seqRegexTrack.setyMaxLines(10); // Sensible default if trackHeightForRegex is unset 
				} else {
					seqRegexTrack.setyMaxLines(trackSet.getTrackHeightForRegex());
				}
			}
			// seqRegexTrack.setyMaxLines(Integer.MAX_VALUE);
			String seqPattern= seqRegexTrack.printToScreen();
			if(!seqPattern.isEmpty()){
				seqPattern+="\n";
			} 
			Utils.printer(seqPattern, snapshotFile, snapshotStripAnsi); 
			// Sequence 
			Utils.printer(gch.current().printableRefSeq(noFormat), snapshotFile, snapshotStripAnsi);
			String ruler= gch.current().printableRuler(10);
			Utils.printer(ruler.substring(0, ruler.length() <= windowSize ? ruler.length() : windowSize) + "\n", snapshotFile, snapshotStripAnsi);

			String footer= gch.current().toString() + "; " + Math.rint(gch.current().getBpPerScreenColumn() * 10d)/10d + " bp/char; " 
					+ "Filters: -q " + mapq  + " -f " + f_incl + " -F " + F_excl
					+ "; " + getMemoryStat();
			if(!noFormat){
				Utils.printer("\033[0;34m" + footer + "\033[0m; \n", snapshotFile, snapshotStripAnsi);
			} else {
				Utils.printer(footer + "\n", snapshotFile, snapshotStripAnsi);
			}
			
			// Optionally convert to png
			if(snapshotFile != null && snapshotFile.endsWith("png")){
				Utils.convertTextFileToGraphic(new File(snapshotFile), new File(snapshotFile));
			}
			
			/* Interactive input */
			/* ================= */
			if(!nonInteractive){
				break;
			}

/*
			boolean needValidInput= true;
			while(needValidInput){
				console.setPrompt("[h] for help: ");
				String cmdInputLong= console.readLine().trim();

				if (cmdInputLong.trim().isEmpty()){
					// Repeat previous command(s)
					cmdInputLong= currentCmd;
				}
				
				ArrayList<String> cmdList= Utils.tokenize(cmdInputLong, "&&");
				
				for(String cmdInput : cmdList){
					
				}
				needValidInput= false;
			}			
*/			
			String cmdInput= null;
			snapshotFile= null; 
			while(cmdInput == null){ // Keep asking for input until you get something valid
				console.setPrompt("[h] for help: ");
				cmdInput = console.readLine().trim();
				
				if (cmdInput.trim().isEmpty()){
					// Repeat previous command
					cmdInput= currentCmd;
				}
				
				if(cmdInput == null || cmdInput.equals("h")){
					Utils.printer(InlineHelp.getHelp() + "\n", snapshotFile, snapshotStripAnsi);
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
						str.setTrimmerMatcher(StrMatcher.spaceMatcher());
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
					} else if(cmdInput.startsWith("colorTrack") || cmdInput.startsWith("colourTrack")){
						try{
							trackSet.setTrackColourForRegex(cmdInput);
						} catch(InvalidCommandLineException e){
							cmdInput= null;
							continue;
						} catch(PatternSyntaxException e) {
							cmdInput= null;
				        	continue;
						} 
					} else if(cmdInput.startsWith("BSseq")) {
						if( fasta == null ){
							System.err.println("Cannot set BSseq mode without fasta");
							cmdInput= null;
							continue;
						}
						trackSet.setBisulfiteModeForRegex(cmdInput);
					} else if (cmdInput.startsWith("squash")){
						trackSet.setFeatureSquashForRegex(cmdInput);
					} else if(cmdInput.startsWith("gffNameAttr")) {
						trackSet.setAttributeForGFFName(cmdInput);
					} else if(cmdInput.startsWith("addTracks ")){
						StrTokenizer str= new StrTokenizer(cmdInput);
						str.setTrimmerMatcher(StrMatcher.spaceMatcher());
						str.setQuoteChar('\'');
						List<String> newFileNames= str.getTokenList();
						newFileNames.remove(0);
						Utils.addTrack(inputFileList, newFileNames);
			
						if(gch.current().getSamSeqDict().size() == 0){
							samSeqDict = GenomicCoords.getSamSeqDictFromAnyFile(inputFileList, null, null);
							GenomicCoords gc= gch.current();
							gc.setSamSeqDict(samSeqDict);
						}
					} else if(cmdInput.startsWith("orderTracks ")){
						StrTokenizer str= new StrTokenizer(cmdInput);
						str.setTrimmerMatcher(StrMatcher.spaceMatcher());
						str.setQuoteChar('\'');
						List<String> newOrder= str.getTokenList();
						newOrder.remove(0);
						trackSet.orderTracks(newOrder);
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
							Utils.printer(x.toString() + "\n", snapshotFile, snapshotStripAnsi);
						}
						cmdInput= null;
					} else if(cmdInput.equals("print") 
							|| cmdInput.startsWith("print ") 
							|| cmdInput.equals("printFull")
							|| cmdInput.startsWith("printFull ")){
						String regex= cmdInput.replaceAll("^printFull|^print", "").trim(); // Memo: ^printFull before ^print 
						try{
							Pattern.compile(regex); // Validate regex
						} catch(Exception e){
					    	System.err.println("Invalid regex in: " + cmdInput);
					    	System.err.println("regex: " + regex);
					    	cmdInput= null;
					    	continue;
						}
						if(regex.isEmpty()){
							regex= ".*";
						}
						if(cmdInput.toLowerCase().equals("print") || cmdInput.toLowerCase().startsWith("print ")){
							printIntervalFeaturesFullRegex= "x^";	// Turn off print/Full
							printIntervalFeaturesRegex= regex;
						} else {
							printIntervalFeaturesRegex= "x^";	// Turn off print/Full
							printIntervalFeaturesFullRegex= regex;
						}
					} else if(cmdInput.toLowerCase().equals("rnameon")){
						withReadName= true;
					} else if(cmdInput.toLowerCase().equals("rnameoff")) {
						withReadName= false;
					} else if(cmdInput.startsWith("next_start ") || cmdInput.equals("next_start")){
						GenomicCoords gc= (GenomicCoords)gch.current().clone();
						GenomicCoords nextGc= trackSet.goToNextFeatureOnFile(cmdInput.replace("next_start", "").trim(), gc, -1.0);
						gch.add(nextGc);
					} else if(cmdInput.startsWith("next ") || cmdInput.equals("next")){
							GenomicCoords gc= (GenomicCoords)gch.current().clone();
							gch.add(trackSet.goToNextFeatureOnFile(cmdInput.replace("next", "").trim(), gc, 5.0));
					} else if(cmdInput.startsWith("find_first ") || 
							  cmdInput.startsWith("find_all ")) {  
						StrTokenizer str= new StrTokenizer(cmdInput);
						str.setTrimmerMatcher(StrMatcher.spaceMatcher());
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
						// Determine whether we match first or all
						boolean all= (cmdInput.startsWith("find_all ")) ? true : false;
						gch.add(trackSet.findNextMatchOnTrack(tokens.get(1), tokens.get(2), gc, all));
					} else if (cmdInput.startsWith("seqRegex")){
						if( fasta == null ){
							System.err.println("Cannot find regex in sequence without fasta reference!");
							cmdInput= null;
							continue;
						}
						StrTokenizer str= new StrTokenizer(cmdInput);
				    	str.setTrimmerMatcher(StrMatcher.spaceMatcher());
						str.setQuoteChar('\'');
						List<String> tokens= str.getTokenList();
						if(tokens.size() == 1){
							seqRegex= "";
						} else {
							seqRegex= tokens.get(1).trim();
							try{
								Pattern.compile(seqRegex);
							} catch(PatternSyntaxException e){
						    	System.err.println("Invalid seqRegex in: " + cmdInput);
						    	System.err.println(e.getDescription());
								cmdInput= null;
								continue;
							}							
						}
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
					} else if(cmdInput.equals("rpm") || cmdInput.startsWith("rpm ")) {
						rpmRegex= cmdInput.replaceAll("^rpm", "").trim();
					} else if(cmdInput.startsWith("-f")) { 
						f_incl= Integer.parseInt(cmdInput.replaceAll("^-f", "").trim());
					} else if(cmdInput.startsWith("-F")) { 
						F_excl= Integer.parseInt(cmdInput.replaceAll("^-F", "").trim());
						if((F_excl & 4) != 4){ // Always filter out read unmapped
							F_excl += 4;
						}
					} else if(cmdInput.startsWith("mapq")) { 
						mapq= Integer.parseInt(cmdInput.replaceAll("^mapq", "").trim());
					}
					// else if(cmdInput.startsWith("maxLines")){
					// maxLines= Integer.parseInt(cmdInput.replaceAll("^maxLines", "").trim());
					//} 
				    else if(cmdInput.startsWith("save")) {
						snapshotFile= Utils.parseCmdinputToGetSnapshotFile(cmdInput, gch.current());
						if(cmdInput.startsWith("save ") || cmdInput.equals("save")){
							snapshotStripAnsi= true;
						} else if(cmdInput.startsWith("savef ") || cmdInput.equals("savef")){
							snapshotStripAnsi= false;
						}
					} else if(cmdInput.startsWith("bookmark")){
						String name= cmdInput.replaceAll("^bookmark", "").trim();
						//GenomicCoords gc = (GenomicCoords)gch.current().clone();
						//trackSet.addBookmark_IN_PREP(gc, name);
					} else {
						System.err.println("Unrecognized argument: " + cmdInput);
						cmdInput= null;
					} // END OF CMD LINE ARGS
				} catch(Exception e){
					System.err.println("\nError processing input: " + cmdInput + "\n");
					e.printStackTrace(); // Print trace for debugging
					System.err.println("");
					cmdInput= null;
				}
				currentCmd= cmdInput;
			} // END while loop to parse cmdInput
			idForTrack= 0;
		} // End while loop keep going until quit or if no interactive input set
	}
}
