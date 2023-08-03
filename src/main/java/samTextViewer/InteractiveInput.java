package samTextViewer;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.regex.Pattern;
import java.util.regex.PatternSyntaxException;

import com.google.common.base.Joiner;
import com.google.common.base.Splitter;
import com.google.common.collect.Lists;

import coloring.Config;
import coloring.ConfigKey;
import commandHelp.Command;
import commandHelp.CommandList;
import exceptions.InvalidColourException;
import exceptions.InvalidCommandLineException;
import exceptions.InvalidConfigException;
import exceptions.InvalidGenomicCoordsException;
import exceptions.InvalidRecordException;
import htsjdk.samtools.SAMSequenceDictionary;
import jline.console.ConsoleReader;
import jline.console.history.History.Entry;
import net.sourceforge.argparse4j.inf.ArgumentParserException;
import utils.Tokenizer;

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
    protected TrackProcessor processInput(String cmdConcatInput, TrackProcessor proc, int debug) throws InvalidGenomicCoordsException, IOException, ClassNotFoundException, InvalidRecordException, SQLException, InvalidCommandLineException{

        cmdConcatInput= cmdConcatInput.replaceAll("//.*", "").trim();
        
        int terminalWidth= Utils.getTerminalWidth();
        
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
            
            //List<String> cmdTokens= Utils.tokenize(cmdString, " ");
            List<String> cmdTokens= new Tokenizer(cmdString).tokenize();
            
            this.interactiveInputExitCode= ExitCode.CLEAN; // If something goes wrong this will change
            try {
                
                // * These commands only print info or do stuff without editing the GenomicCoordinates or the Tracks:
                if(cmdTokens.size() == 1 && 
                        (cmdTokens.get(0).equals("h") || 
                        cmdTokens.get(0).equals("-h") || 
                        cmdTokens.get(0).equals("help") || 
                        cmdTokens.get(0).equals("?"))){          
                    System.err.println(Utils.padEndMultiLine(CommandList.briefHelp(), proc.getWindowSize()));
                    this.interactiveInputExitCode= ExitCode.CLEAN_NO_FLUSH;
                    
                } else if( (cmdTokens.size() >= 2 && cmdTokens.get(1).equals("-h")) || 
                            (cmdTokens.size() >= 2 && cmdTokens.get(0).equals("help")) ||
                            cmdTokens.get(0).startsWith("?") ){ 
                    // Help on this command
                    String cmd;
                    if(cmdTokens.size() >= 2 && cmdTokens.get(0).equals("-h")){
                        cmd= cmdTokens.get(0);
                    } else if( cmdTokens.size() >= 2 && cmdTokens.get(0).equals("help") ){
                        cmd= cmdTokens.get(1);
                    } else {
                        cmd= cmdTokens.get(0).replaceAll("^\\?", "");    
                    }
                    
                    String help= Utils.padEndMultiLine("\n" + CommandList.getHelpForCommand(cmd), proc.getWindowSize());
                    System.err.println(help);
                    this.interactiveInputExitCode= ExitCode.CLEAN_NO_FLUSH;
                    
                } else if(cmdTokens.get(0).equals("posHistory")){
                    this.posHistory(cmdTokens, proc.getGenomicCoordsHistory().getCurrentSessionHistory(), proc.getWindowSize());
                    this.interactiveInputExitCode= ExitCode.CLEAN_NO_FLUSH;

                } else if(cmdTokens.get(0).equals("history")){
                    String hist= Utils.padEndMultiLine(this.cmdHistoryToString(cmdTokens), proc.getWindowSize());
                    System.err.println(hist);
                    this.interactiveInputExitCode= ExitCode.CLEAN_NO_FLUSH;
                
                } else if(cmdTokens.get(0).equals("show")){
                    this.interactiveInputExitCode= this.show(cmdTokens, proc);
                                    
                } else if(cmdTokens.get(0).equals("explainSamFlag")){
                    this.interactiveInputExitCode= this.explainSamFlag(cmdTokens, proc);
                    
                } else if(cmdTokens.get(0).equals("sys")) {
                    this.execSysCmd(cmdString, proc.getWindowSize());
                    this.interactiveInputExitCode= ExitCode.CLEAN_NO_FLUSH;
                    
                } else if(cmdTokens.get(0).equals("recentlyOpened")) {
                    String opened= Utils.padEndMultiLine(proc.getTrackSet().showRecentlyOpened(cmdTokens), proc.getWindowSize());
                    System.out.println(opened);
                    this.interactiveInputExitCode= ExitCode.CLEAN_NO_FLUSH;
                
                } else if(cmdTokens.get(0).equals("setConfig")) {
                    try{
                        this.setConfigOpt(cmdTokens);
                        this.interactiveInputExitCode= ExitCode.CLEAN;
                    } catch(Exception e){
                        System.err.println(Utils.padEndMultiLine("Unable to set configuration", proc.getWindowSize()));
                        this.interactiveInputExitCode= ExitCode.ERROR;
                        if(debug > 0){
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
                    
                } else if(cmdTokens.get(0).equals("q")){
                    System.out.print("\033[0m");
                    console.clearScreen();
                    console.flush();
                    System.exit(0);
                
                // * These commands change the GenomicCoordinates (navigate) but do not touch the tracks.
                } else if(cmdTokens.get(0).equals("f")
                        || cmdTokens.get(0).equals("b")
                        || cmdTokens.get(0).equals("ff") 
                        || cmdTokens.get(0).equals("bb")
                        || cmdTokens.get(0).matches("^\\-\\d+.*") 
                        || cmdTokens.get(0).matches("^\\+\\d+.*")){ // No cmd line args either f/b ops or ints
                        String newRegion= Utils.parseConsoleInput(cmdTokens, proc.getGenomicCoordsHistory().current()).trim();
                        proc.getGenomicCoordsHistory().add(new GenomicCoords(newRegion, terminalWidth, samSeqDict, fasta));
                
                } else if(cmdTokens.get(0).matches("(\\[+|\\]+)\\d*")) {
                        
                        int times= this.countBrackets(cmdTokens);
                        String newRegion= this.moveWindowByColumns(proc.getGenomicCoordsHistory().current(), times);
                        proc.getGenomicCoordsHistory().add(new GenomicCoords(newRegion, terminalWidth, samSeqDict, fasta));
                        
                } else if(cmdTokens.get(0).matches("^\\d+.*") || cmdTokens.get(0).matches("^\\.\\d+.*")){
                    String newRegion;
                    try{
                        newRegion= this.gotoOnCurrentChrom(cmdTokens, proc.getGenomicCoordsHistory().current());
                    } catch(IndexOutOfBoundsException e){
                        System.err.append("Column coordinates must be >= 1 and <= the screen width");
                        throw new InvalidCommandLineException();
                    }
                    proc.getGenomicCoordsHistory().add(new GenomicCoords(newRegion, terminalWidth, samSeqDict, fasta));
                    
                } else if(cmdTokens.get(0).equals("goto") || cmdTokens.get(0).startsWith(":")){
                    String reg= Joiner.on(" ").join(cmdTokens).replaceFirst("goto|:", "").trim();
                    proc.getGenomicCoordsHistory().add(new GenomicCoords(reg, terminalWidth, samSeqDict, fasta));

                } else if(cmdTokens.get(0).equals("nextChrom")){
                    this.interactiveInputExitCode= this.nextChrom(cmdTokens, proc);
                    
                } else if (cmdTokens.get(0).equals("p")) {
                    proc.getGenomicCoordsHistory().previous(); 
                    
                } else if (cmdTokens.get(0).equals("n")) {
                    proc.getGenomicCoordsHistory().next();
                    
                } else if(cmdTokens.get(0).equals("zo")){
                    int nz= Utils.parseZoom(Joiner.on(" ").join(cmdTokens), 1);
                    GenomicCoords gc = (GenomicCoords)proc.getGenomicCoordsHistory().current().clone();
                    gc.setTerminalWidth(terminalWidth);
                    for(int i= 0; i < nz; i++){
                        gc.zoomOut();
                    }
                    proc.getGenomicCoordsHistory().add(gc);
                    
                } else if(cmdTokens.get(0).equals("zi")){
                    int nz= Utils.parseZoom(Joiner.on(" ").join(cmdTokens), 1);
                    GenomicCoords gc = (GenomicCoords)proc.getGenomicCoordsHistory().current().clone();
                    gc.setTerminalWidth(terminalWidth);
                    for(int i= 0; i < nz; i++){
                        gc.zoomIn();
                    }
                    proc.getGenomicCoordsHistory().add(gc);

                } else if(cmdTokens.get(0).equals("extend")) {
                    if(cmdTokens.size() == 1){
                        System.err.println(Utils.padEndMultiLine("Expected at least one argument.", proc.getWindowSize()));
                    }
                    GenomicCoords gc = (GenomicCoords)proc.getGenomicCoordsHistory().current().clone();
                    gc.setTerminalWidth(terminalWidth);
                    gc.cmdInputExtend(cmdTokens);
                    proc.getGenomicCoordsHistory().add(gc);
                
                } else if(cmdTokens.get(0).equals("trim")){
                    GenomicCoords gc = proc.getTrackSet().trimCoordsForTrack(cmdTokens);
                    proc.getGenomicCoordsHistory().add(gc);
                    
                } else if(cmdTokens.get(0).equals("l")) {
                    GenomicCoords gc = (GenomicCoords)proc.getGenomicCoordsHistory().current().clone();
                    gc.setTerminalWidth(terminalWidth);
                    gc.left();
                    proc.getGenomicCoordsHistory().add(gc);
                
                } else if(cmdTokens.get(0).equals("r")) {
                    GenomicCoords gc = (GenomicCoords)proc.getGenomicCoordsHistory().current().clone();
                    gc.setTerminalWidth(terminalWidth);
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
                    
                } else if((cmdTokens.get(0).equals("colorTrack") || cmdTokens.get(0).equals("colourTrack"))){
                    proc.getTrackSet().setTrackColourForRegex(cmdTokens);
                    
                } else if(cmdTokens.get(0).equals("bedToBedgraph")) {
                    proc.getTrackSet().setTrackFormatForRegex(cmdTokens);

                } else if((cmdTokens.get(0).equals("featureColor"))){
                    proc.getTrackSet().setFeatureColorForRegex(cmdTokens); 
                    
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
                
                } else if (cmdTokens.get(0).equals("readsAsPairs")){
                    proc.getTrackSet().setReadsAsPairsForRegex(cmdTokens);
                    
                } else if(cmdTokens.get(0).equals("nameForFeatures")) {
                    proc.getTrackSet().setNameAttribute(cmdTokens);
                    
                } else if(cmdTokens.get(0).equals("open") || cmdTokens.get(0).equals("addTracks")){
                    cmdTokens.remove(0);
                    
                    List<String> globbed = Utils.globFiles(cmdTokens);
                    if(globbed.size() == 0){
                        globbed= this.openFilesFromIndexes(proc.getTrackSet().getOpenedFiles(), cmdTokens);
                    }
                    
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
                                    String region= Main.initRegion(globbed, null, null, debug);
                                    proc.getGenomicCoordsHistory().add(new GenomicCoords(region, terminalWidth, samSeqDict, fasta));
                                    proc.getTrackSet().addTrackFromSource(sourceName, proc.getGenomicCoordsHistory().current(), null);
                                } catch (Exception x){
                                    x.printStackTrace();
                                    msg= Utils.padEndMultiLine("Failed to add: " + sourceName, proc.getWindowSize());
                                    System.err.println(msg);
                                    this.interactiveInputExitCode= ExitCode.CLEAN_NO_FLUSH;
                                }
                            }

                            if(proc.getGenomicCoordsHistory().current().getSamSeqDict() == null || proc.getGenomicCoordsHistory().current().getSamSeqDict().size() == 0){
                                GenomicCoords testSeqDict= new GenomicCoords("default", Utils.getTerminalWidth(), null, null); 
                                List<String> candidateSourceGenome= new ArrayList<String>();
                                candidateSourceGenome.add(sourceName);
                                testSeqDict.setGenome(candidateSourceGenome, false);
                                if(testSeqDict.getSamSeqDict() != null){
                                    candidateSourceGenome.add(0, "cmd");
                                    proc.getGenomicCoordsHistory().setGenome(candidateSourceGenome);
                                }
                            }
                        }
                    }
                } 
                else if(cmdTokens.get(0).equals("reload")){
                    proc.getTrackSet().reload(cmdTokens);
                }
                else if(cmdTokens.get(0).equals("dropTracks")){
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
                    proc.getTrackSet().setAwkForTrack(cmdTokens);
                    
                } else if(cmdTokens.get(0).equals("filterVariantReads")){
                    proc.getTrackSet().setFilterVariantReads(cmdTokens);
                    
                } else if(cmdTokens.get(0).equals(Command.rpm.getCmdDescr())) {
                    proc.getTrackSet().setRpmForRegex(cmdTokens);

                } else if(cmdTokens.get(0).equals("samtools")){
                    proc.getTrackSet().setSamFilterForRegex(cmdTokens);
                    
                } else if(cmdTokens.get(0).equals("genotype")){
                    proc.getTrackSet().setGenotypeMatrix(cmdTokens);
            
                // * These commands change both the Tracks and the GenomicCoordinates
                } else if(cmdTokens.get(0).equals("next")){
                    
                    this.next(cmdTokens, proc);
                                        
                } else if(cmdTokens.get(0).equals("find")) {  
                    
                    boolean all= Utils.argListContainsFlag(cmdTokens, "-all");
                    boolean fixedPattern= Utils.argListContainsFlag(cmdTokens, "-F");
                    boolean caseIns= Utils.argListContainsFlag(cmdTokens, "-c");
                    
                    if(cmdTokens.size() < 2){
                        System.err.println(Utils.padEndMultiLine("Error in find command. Expected at least 1 argument got: " + cmdTokens, proc.getWindowSize()));
                        this.interactiveInputExitCode= ExitCode.ERROR;
                        continue;
                    }
                    if(cmdTokens.size() == 2){
                        cmdTokens.add(""); // If track arg is missing use this placeholder.
                    }
                    GenomicCoords gc= (GenomicCoords)proc.getGenomicCoordsHistory().current().clone();
                    gc.setTerminalWidth(terminalWidth);
                    
                    int flag= 0;
                    if(fixedPattern){
                        flag |= Pattern.LITERAL;
                    }
                    if( ! caseIns){
                        flag |= Pattern.CASE_INSENSITIVE;
                    }
                    Pattern pattern;
                    try{
                        pattern= Pattern.compile(cmdTokens.get(1), flag);
                    } catch(PatternSyntaxException e){
                        System.err.println("Invalid regex");
                        throw new InvalidCommandLineException();
                    }
                    GenomicCoords nextGc= proc.getTrackSet().findNextMatchOnTrack(pattern, cmdTokens.get(2), gc, all);
                    if(nextGc.equalCoords(gc)){
                        System.err.println("No match found outside of this window for query '" + cmdTokens.get(1) + "'");
                        this.interactiveInputExitCode= ExitCode.CLEAN_NO_FLUSH;
                    } else {
                        proc.getGenomicCoordsHistory().add(nextGc);
                    }
                    
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
                    System.err.println(Utils.padEndMultiLine("Unrecognized command: " + cmdTokens.get(0), proc.getWindowSize()));
                    String suggestions= Joiner.on(" or ").join(Utils.suggestCommand(cmdTokens.get(0).trim(), CommandList.cmds()));
                    if( ! suggestions.isEmpty()){
                        System.err.println(Utils.padEndMultiLine("Maybe you mean " + suggestions + "?", proc.getWindowSize()));
                    }
                    this.interactiveInputExitCode= ExitCode.ERROR;
                    // throw new InvalidCommandLineException();
                }
            } catch(ArgumentParserException e){
                this.interactiveInputExitCode= ExitCode.ERROR;
            
            } catch(Exception e){ // You shouldn't catch anything! Be more specific.
                System.err.println(Utils.padEndMultiLine("\nError processing input: " + cmdTokens, proc.getWindowSize()));
                System.err.println(Utils.padEndMultiLine("For help on command \"cmd\" execute 'cmd -h' or '-h' for list of commands.\n", proc.getWindowSize()));
                this.interactiveInputExitCode= ExitCode.ERROR; 
                if(debug == 1){
                    e.printStackTrace();
                } else if(debug == 2){
                    e.printStackTrace();
                    System.exit(1);
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
                    String region= Main.initRegion(proc.getTrackSet().getFilenameList(), null, null, debug);
                    proc.getGenomicCoordsHistory().add(new GenomicCoords(region, terminalWidth, samSeqDict, fasta));
                    System.err.println(Utils.padEndMultiLine("Invalid genomic coordinates found. Resetting to "  + region, proc.getWindowSize()));
                    if(debug > 0){
                        e.printStackTrace();
                    }
                    
                } catch (Exception e){
                    System.err.println(Utils.padEndMultiLine("Error processing tracks with input " + cmdTokens, proc.getWindowSize()));
                    this.interactiveInputExitCode= ExitCode.ERROR;
                    if(debug > 0){
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

    private ExitCode nextChrom(List<String> cmdTokens, TrackProcessor proc) throws IOException, NumberFormatException, InvalidCommandLineException {
        
        int minSize = Integer.parseInt(Utils.getArgForParam(cmdTokens, "-m", "-1"));
        int maxSize = Integer.parseInt(Utils.getArgForParam(cmdTokens, "-M", "-1"));
        String regex = Utils.getArgForParam(cmdTokens, "-r", ".*");
        String so = Utils.getArgForParam(cmdTokens, "-s", "u");
        ContigOrder sortOrder;
        if(so.equals("u")) {
            sortOrder = ContigOrder.UNSORTED;
        } else if(so.equals("s")) {
            sortOrder = ContigOrder.SIZE_ASC;
        } else if(so.equals("S")) {
            sortOrder = ContigOrder.SIZE_DESC;
        } else {
            System.err.println("Invalid option for sort order: " + so);
            return ExitCode.CLEAN_NO_FLUSH;
        }
        
        try {
            GenomicCoords gc= (GenomicCoords)proc.getGenomicCoordsHistory().current().clone();
            gc.nextChrom(proc.getTrackSet().getKnownContigs(), minSize, maxSize, regex, sortOrder);
            proc.getGenomicCoordsHistory().add(gc);
            return ExitCode.CLEAN;
        } catch(InvalidGenomicCoordsException e) {
            System.err.println(e.getMessage());
            return ExitCode.CLEAN_NO_FLUSH;
        }
    }

    private int countBrackets(List<String> cmdTokens) throws InvalidCommandLineException {
        
        String xtimes= cmdTokens.get(0).replaceAll("\\[|\\]", "");
        if(! xtimes.isEmpty() && cmdTokens.size() > 1) {
            throw new InvalidCommandLineException();
        }
        int times= 1;
        if( ! xtimes.isEmpty()) {
            times= Integer.valueOf(xtimes);
        }
        else if(cmdTokens.size() > 1) {
            try {
                times= Integer.valueOf(cmdTokens.get(1));
            } catch(NumberFormatException e) {
                System.err.append(cmdTokens.get(1) + " cannot be converted to integer");
                throw new InvalidCommandLineException();
            }
        }
        String brackets= cmdTokens.get(0).replaceAll("\\d", "");
        times = times * brackets.length(); 

        if(cmdTokens.get(0).startsWith("]")) {
            times *= -1;
        }
        return times;
    }

    /** Return a string where the current genomic coordinates are moved forward or 
     * backwards "times" screen columns. Move backwards if times is negative.
     * @throws IOException 
     * @throws InvalidGenomicCoordsException 
     * @throws InvalidCommandLineException 
     * */
    private String moveWindowByColumns(GenomicCoords gc, int times) throws InvalidGenomicCoordsException, IOException, InvalidCommandLineException {
        String x = String.valueOf(Math.round(gc.getBpPerScreenColumn() * times));
        if(Integer.valueOf(x) >= 0) {
            x= "+" + x;
        }
        List<String> tokens= new ArrayList<String>();
        tokens.add(String.valueOf(x));
        return Utils.parseConsoleInput(tokens, gc);
    }

    private ExitCode explainSamFlag(List<String> cmdTokens, TrackProcessor proc) throws InvalidCommandLineException, IOException {
        List<String> args= new ArrayList<String>(cmdTokens);
        args.remove(0);
        if(args.size() == 0){
            System.err.println("One argument is required.");
            throw new InvalidCommandLineException();
        }
        List<Integer> flagsToDecode= new ArrayList<Integer>();
        for(String x : args){
            try{
                flagsToDecode.add(Integer.valueOf(x));
            } catch (NumberFormatException e){
                System.err.println(Utils.padEndMultiLine("Expected an integer, got " + x, Utils.getTerminalWidth()));
                return ExitCode.ERROR;
            }            
        }
        final Map<Integer, String> FLAGS= new LinkedHashMap<Integer, String>();
        
        FLAGS.put(1, "read paired [1]");
        FLAGS.put(2, "read mapped in proper pair [2]");
        FLAGS.put(4, "read unmapped [4]");
        FLAGS.put(8, "mate unmapped [8]");
        FLAGS.put(16, "read reverse strand [16]");
        FLAGS.put(32, "mate reverse strand [32]");
        FLAGS.put(64, "first in pair [64]");
        FLAGS.put(128, "second in pair [128]");
        FLAGS.put(256, "not primary alignment [256]");
        FLAGS.put(512, "read fails platform/vendor quality checks [512]");
        FLAGS.put(1024, "read is PCR or optical duplicate [1024]");
        FLAGS.put(2048, "supplementary alignment [2048]");
        
        String LEFT_PAD= " ";
        List<String> table= new ArrayList<String>();
        table.add(LEFT_PAD + Joiner.on("\t").join(flagsToDecode));
        
        System.out.println(Utils.padEndMultiLine("", Utils.getTerminalWidth()));
        for(Integer i : FLAGS.keySet()){
            String line= LEFT_PAD;
            for(int f : flagsToDecode){
                if((i & f) == i) {
                    line += "X";
                } else {
                    line += ".";
                }
                line += "\t";
            }
            line += FLAGS.get(i);
            table.add(line);
        }
        List<String> xtable= Utils.tabulateList(table, Utils.getTerminalWidth(), " ");
        System.out.println(Utils.padEndMultiLine(Joiner.on("\n").join(xtable), Utils.getTerminalWidth()));
        System.out.println(Utils.padEndMultiLine("", Utils.getTerminalWidth()));
        return ExitCode.CLEAN_NO_FLUSH;
    }

    /**Get the items (files) corresponding to the indexes. Errors are silently ignored.
     * */
    private List<String> openFilesFromIndexes(LinkedHashSet<String> openedFiles, List<String> indexes) {
        List<String> files= new ArrayList<String>();
        List<Integer> idxs= new ArrayList<Integer>();
        for(String x : indexes){
            try{
                idxs.add(Integer.parseInt(x));
            } catch(NumberFormatException e){
                // Return empty list
                return files;
                //
            }
        }
        for(int i : idxs){
            try{
                String x= Lists.reverse(new ArrayList<String>(openedFiles)).get(i - 1);
                files.add(x);
            } catch (Exception e){
                //
            }
        }
        return files;
        }

    /** Parse the given list of options and move to new coordinates.
     * Visibility set to protected only for testing.
     * */
    protected String gotoOnCurrentChrom(List<String> cmdTokens, GenomicCoords gc) throws InvalidGenomicCoordsException, IOException {
        List<String> args= new ArrayList<String>(cmdTokens);
        int regFrom;
        int regTo;
        
        boolean center= false;
        if(args.get(0).endsWith("c")){
            center= true;
            args.set(0, args.get(0).replaceAll("c$", ""));
        }
        if(args.size() > 1 && args.get(1).equals("c")){
            center= true;
            args.remove(1);
        }
        
        if(Float.valueOf(args.get(0)) < 1){
            // Switch to screen percent coordinates
            float pctFrom= Float.parseFloat(args.get(0));
            // Which is the screen column matching this percent?
            int screenIdx= (int) Math.rint(pctFrom * gc.getUserWindowSize());
            // Which is the genomic coordinate corresponding to this screen index?
            regFrom= (int) Math.rint(gc.getMapping().get(screenIdx));

            // Same for end position. Accounting for possibility that only one pct value is given
            if(args.size() > 1 && ! center){
                float pctTo= Float.parseFloat(args.get(1));
                if(pctTo > 1){
                    System.err.println("Coordinate interval given as percent of screen width must be between 0 and 1. Got " + args);
                    throw new InvalidGenomicCoordsException();
                }
                screenIdx= (int) Math.rint(pctTo * gc.getUserWindowSize());
                if(screenIdx >= gc.getMapping().size()){
                    screenIdx= gc.getMapping().size() - 1;
                }
                regTo= (int) Math.rint(gc.getMapping().get(screenIdx));
            
            } 
            else if(center){
                regFrom= regFrom - (int)Math.floor(gc.getGenomicWindowSize()/2) - 1;
                regFrom= regFrom < 1 ? 1 : regFrom;
                regTo= regFrom + gc.getGenomicWindowSize() - 1;
            } 
            else {
                regTo= regFrom + gc.getGenomicWindowSize() - 1;
            }
            if(regTo < regFrom){
                System.err.println("Invalid coordinates: end < start for argument(s): " + args);
                throw new InvalidGenomicCoordsException();                
            }
            //return gc.getChrom() + ":" + regFrom + "-" + regTo; 
        }
        else if(args.size() == 1 && ! center){
            regFrom= Integer.valueOf(args.get(0));
            regTo= regFrom + gc.getUserWindowSize() - 1;
        }
        else if(! center){
            regFrom= Integer.valueOf(args.get(0));
            regTo= Integer.valueOf(args.get(args.size()-1));
        }
        else if(center){
            regFrom= Integer.valueOf(args.get(0)) - (int)Math.floor(gc.getGenomicWindowSize()/2);
            regFrom= regFrom < 1 ? 1 : regFrom;
            regTo= regFrom + gc.getGenomicWindowSize() - 1;
        }
        else {
            throw new IllegalArgumentException();
        }
        String newRegion= gc.getChrom()+ ":" + regFrom + "-" + regTo;
        return newRegion;
    }

    private void setConfigOpt(List<String> cmdTokens) throws IOException, InvalidConfigException, InvalidCommandLineException, InvalidColourException {
        List<String> args=  new ArrayList<String>(cmdTokens);
        args.remove(0);
        if(args.size() == 0){
            throw new InvalidCommandLineException();
        } 
        if(args.size() == 1){
            ConfigKey key= ConfigKey.getConfigKeyFromShort(args.get(0));
            if(ConfigKey.booleanKeys().contains(key)){
                // If configkey is a type boolean, just flip the boolean
                key = ConfigKey.valueOf(key.toString());
                boolean value= ! Utils.asBoolean(Config.get(key));
                Config.set(key, String.valueOf(value));
            } else {
                // configKey is expected to be the name of a configuration file
                new Config(args.get(0));
            }
        } else {
            ConfigKey key= ConfigKey.getConfigKeyFromShort(args.get(0));
            String value= args.get(1);
            try{
                Config.set(key, value);
            } catch(Exception e){
                throw new InvalidConfigException();
            }
        }
    }

    private void setGenome(List<String> cmdTokens, TrackProcessor proc) throws InvalidGenomicCoordsException, IOException, InvalidCommandLineException {
        
        List<String> tokens = new ArrayList<String>(cmdTokens); 
        tokens.remove(0);
        
        if(tokens.size() == 0){
            // Try to read fasta from history file
            ASCIIGenomeHistory ag= new ASCIIGenomeHistory();
            try{
                tokens.add(ag.getReference().get(0));
                System.err.println("Using " + ag.getReference().get(0));
            } catch(Exception e){
                System.err.println("A previous reference file was not found.");
                throw new InvalidCommandLineException();    
            }
        }
        
        GenomicCoords testSeqDict= new GenomicCoords("default", Utils.getTerminalWidth(), null, null); 
        testSeqDict.setGenome(tokens, true);
        if(testSeqDict.getSamSeqDict() != null){
            proc.getGenomicCoordsHistory().setGenome(tokens);
        } else {
            System.err.println(Utils.padEndMultiLine("Cannot set genome from " + tokens, Utils.getTerminalWidth()));
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
            BufferedReader reader= new BufferedReader(new InputStreamReader(p.getInputStream()));

            String line = "";
            while ((line = reader.readLine())!= null) {
                System.err.println(Utils.padEndMultiLine(line, userWindowSize));
            }
            p.waitFor();

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
        boolean start= false;
        if(args.contains("-start")){
            start= true;
            args.remove("-start");
        }
        boolean center= false;
        if(args.contains("-c")){
            center= true;
            args.remove("-c");
        }
        boolean getPrevious= false;
        if(args.contains("-back")){
            getPrevious= true;
            args.remove("-back");
        }
        String trackId= "";
        if(args.size() > 0){
            trackId= args.get(0);
        }
        if(start){
            proc.getGenomicCoordsHistory().add(proc.getTrackSet().goToNextFeatureOnFile(trackId, gc, -1.0, getPrevious));
        } else if(center){
            proc.getGenomicCoordsHistory().add(proc.getTrackSet().goToNextFeatureOnFile(trackId, gc, 0, getPrevious));
        } else {
            proc.getGenomicCoordsHistory().add(proc.getTrackSet().goToNextFeatureOnFile(trackId, gc, zo, getPrevious));
        }
        
    }
    
    /**Edit visualization setting in TrackProcessor as appropriate.
     * @return 
     * @throws InvalidCommandLineException 
     * @throws IOException 
     * @throws InvalidGenomicCoordsException 
     * */
    private ExitCode show(List<String> cmdTokens, TrackProcessor proc) throws InvalidCommandLineException, InvalidGenomicCoordsException, IOException {
        List<String> args= new ArrayList<String>(cmdTokens);
        args.remove(0);
        if(args.size() == 0){
            System.err.println("At least one argument is required.");
            throw new InvalidCommandLineException();
        }
        // With .startsWith() we allow partial matching of input to argument. I.e. "ge" will be enough to 
        // recognize "genome".
        if("genome".startsWith(args.get(0))){
            this.showGenome(proc);
            return ExitCode.CLEAN_NO_FLUSH;
        
        } else if("trackInfo".startsWith(args.get(0))){
            String info = Utils.padEndMultiLine(proc.getTrackSet().showTrackInfo(), proc.getWindowSize());
            System.err.println(info);
            return ExitCode.CLEAN_NO_FLUSH;
        
        } else if("gruler".startsWith(args.get(0))){
            proc.setShowGruler(! proc.isShowGruler());
            return ExitCode.CLEAN;

        } else if("pctRuler".startsWith(args.get(0))){
            proc.setShowCruler(! proc.isShowCruler());
            return ExitCode.CLEAN;
            
        } else {
            System.err.println("Unrecognized option: " + args.get(0));
            throw new InvalidCommandLineException();
        }
    }
    
    /**Print to screen the genome file*/
    private void showGenome(TrackProcessor proc) throws InvalidGenomicCoordsException, IOException {
        String genome= Utils.printSamSeqDict(proc.getGenomicCoordsHistory().current().getSamSeqDict(), 30);
        if(genome != null && ! genome.isEmpty()){
            System.err.println(Utils.padEndMultiLine(genome, proc.getWindowSize()));
            return;
        }
        
        List<String> chroms= proc.getTrackSet().getKnownContigs();
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
        
        List<String> cmdTab= Utils.tabulateList(cmd, -1, " ");
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
