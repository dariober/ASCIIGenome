package samTextViewer;

import static session.SessionHandler.writeSessions;

import colouring.Config;
import colouring.ConfigKey;
import com.google.common.base.Joiner;
import com.google.common.base.Splitter;
import com.google.common.collect.Lists;
import commandHelp.Command;
import commandHelp.CommandList;
import exceptions.*;
import htsjdk.samtools.SAMSequenceDictionary;
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.*;
import java.util.regex.Pattern;
import java.util.regex.PatternSyntaxException;
import jline.console.ConsoleReader;
import jline.console.history.History.Entry;
import net.sourceforge.argparse4j.inf.ArgumentParserException;
import org.biojava.nbio.core.sequence.io.IUPACParser;
import org.biojava.nbio.core.sequence.transcription.Frame;
import session.Session;
import session.SessionHandler;
import tracks.Track;
import utils.Tokenizer;

/** Class to process input from console */
public class InteractiveInput {

  private boolean nonInteractive;
  private ExitCode interactiveInputExitCode = ExitCode.CLEAN;
  private final ConsoleReader console;
  private static int debug;
  private SAMSequenceDictionary samSeqDict;
  private String fasta;
  private List<String> messages =
      new ArrayList<>(); // Messages that may be sent from the various methods.
  private SessionHandler sessionHandler;
  private String currentSessionName;

  public InteractiveInput(ConsoleReader console, int debug) {
    InteractiveInput.debug = debug;
    this.console = console;
  }

  /** Parse the input list of commands to print information or modify the input TrackProcessor. */
  protected TrackProcessor processInput(String cmdConcatInput, TrackProcessor proc)
      throws InvalidGenomicCoordsException, IOException {
    cmdConcatInput = cmdConcatInput.replaceAll("//.*", "").trim();
    int terminalWidth = Utils.getTerminalWidth();
    // cmdInputList: List of individual commands in tokens to be issued.
    // E.g.: [ ["zi"],
    //         ["-F", "16"],
    //         ["mapq", "10"] ]
    // Don't check the validity of each cmd now. Execute one by one and if anything goes wrong
    // reset interactiveInputExitCode = 1 (or else other than 0) so that console input is asked
    // again. Of course, what is executed is not
    // rolled back.
    List<String> cmdInputChainList = new ArrayList<String>();

    // See
    // http://stackoverflow.com/questions/1757065/java-splitting-a-comma-separated-string-but-ignoring-commas-in-quotes
    // For splitting at delimiter (&&) and ignore delimiters inside single quotes.
    for (String cmd :
        Splitter.on(Pattern.compile("&&(?=([^']*'[^']*')*[^']*$)"))
            .trimResults()
            .omitEmptyStrings()
            .split(cmdConcatInput)) {
      cmdInputChainList.add(cmd);
    }
    if (cmdInputChainList.size() >= 2 && cmdInputChainList.get(0).startsWith("setConfig ")) {
      cmdInputChainList.add("+0"); // This is to refresh the screen and actually set the new colour
    }

    this.fasta = proc.getGenomicCoordsHistory().current().getFastaFile();
    this.samSeqDict = proc.getGenomicCoordsHistory().current().getSamSeqDict();

    for (String cmdString : cmdInputChainList) {

      List<String> cmdTokens = new Tokenizer(cmdString).tokenize();
      if (!cmdTokens.isEmpty()) {
        cmdTokens.set(0, cmdTokens.get(0).replaceAll("color", "colour"));
        cmdTokens.set(0, cmdTokens.get(0).replaceAll("Color", "Colour"));
      }

      this.interactiveInputExitCode = ExitCode.CLEAN; // If something goes wrong this will change
      try {

        // * These commands only print info or do stuff without editing the GenomicCoordinates or
        // the Tracks:
        if (cmdTokens.size() == 1
            && (cmdTokens.get(0).equals("h")
                || cmdTokens.get(0).equals("-h")
                || cmdTokens.get(0).equals("help")
                || cmdTokens.get(0).equals("?"))) {
          System.err.println(Utils.padEndMultiLine(CommandList.briefHelp(), proc.getWindowSize()));
          this.interactiveInputExitCode = ExitCode.CLEAN_NO_FLUSH;

        } else if ((cmdTokens.size() >= 2 && cmdTokens.get(1).equals("-h"))
            || (cmdTokens.size() >= 2 && cmdTokens.get(0).equals("help"))
            || cmdTokens.get(0).startsWith("?")) {
          // Help on this command
          String cmd;
          if (cmdTokens.size() >= 2 && cmdTokens.get(0).equals("-h")) {
            cmd = cmdTokens.get(0);
          } else if (cmdTokens.size() >= 2 && cmdTokens.get(0).equals("help")) {
            cmd = cmdTokens.get(1);
          } else {
            cmd = cmdTokens.get(0).replaceAll("^\\?", "");
          }

          String help =
              Utils.padEndMultiLine(
                  "\n" + CommandList.getHelpForCommand(cmd), proc.getWindowSize());
          System.err.println(help);
          this.interactiveInputExitCode = ExitCode.CLEAN_NO_FLUSH;

        } else if (cmdTokens.get(0).equals("posHistory")) {
          this.posHistory(
              cmdTokens,
              proc.getGenomicCoordsHistory().getCurrentSessionHistory(),
              proc.getWindowSize());
          this.interactiveInputExitCode = ExitCode.CLEAN_NO_FLUSH;

        } else if (cmdTokens.get(0).equals("history")) {
          String hist =
              Utils.padEndMultiLine(this.cmdHistoryToString(cmdTokens), proc.getWindowSize());
          System.err.println(hist);
          this.interactiveInputExitCode = ExitCode.CLEAN_NO_FLUSH;

        } else if (cmdTokens.get(0).equals("show")) {
          this.interactiveInputExitCode = this.show(cmdTokens, proc);

        } else if (cmdTokens.get(0).equals("explainSamFlag")) {
          this.interactiveInputExitCode = this.explainSamFlag(cmdTokens, proc);

        } else if (cmdTokens.get(0).equals("sys")) {
          this.execSysCmd(cmdString, proc.getWindowSize());
          this.interactiveInputExitCode = ExitCode.CLEAN_NO_FLUSH;

        } else if (cmdTokens.get(0).equals("recentlyOpened")) {
          String opened =
              Utils.padEndMultiLine(
                  proc.getTrackSet().showRecentlyOpened(cmdTokens), proc.getWindowSize());
          System.out.println(opened);
          this.interactiveInputExitCode = ExitCode.CLEAN_NO_FLUSH;

        } else if (cmdTokens.get(0).equals("setConfig")) {
          try {
            this.interactiveInputExitCode = this.setConfigOpt(cmdTokens, proc);
          } catch (Exception e) {
            if (debug > 0) {
              e.printStackTrace();
            }
          }

        } else if (cmdTokens.get(0).equals("save")) {
          List<String> args = new ArrayList<String>(cmdTokens);

          proc.setStripAnsi(!Utils.argListContainsFlag(args, "-colour"));

          proc.setAppendToSnapshotFile(false); // Default: do not append

          if (args.contains(">>")) {
            proc.setAppendToSnapshotFile(true);
            args.remove(">>");
          } else if (args.contains(">")) {
            proc.setAppendToSnapshotFile(false);
            args.remove(">");
          }
          proc.setSnapshotFile(
              Utils.parseCmdinputToGetSnapshotFile(
                  Joiner.on(" ").join(args), proc.getGenomicCoordsHistory().current()));

        } else if (cmdTokens.get(0).equals("q")) {
          System.out.print("\033[0m");
          console.clearScreen();
          console.flush();
          try {
            SessionHandler.saveAs(SessionHandler.DEFAULT_SESSION_FILE, proc.toSession(), "");
          } catch (Exception e) {
            System.err.println(
                "Error saving session to file: " + SessionHandler.DEFAULT_SESSION_FILE);
            System.err.println(e.getMessage());
            System.exit(1);
          }
          System.exit(0);

          // * These commands change the GenomicCoordinates (navigate) but do not touch the tracks.
        } else if (cmdTokens.get(0).equals("f")
            || cmdTokens.get(0).equals("b")
            || cmdTokens.get(0).equals("ff")
            || cmdTokens.get(0).equals("bb")
            || cmdTokens.get(0).matches("^\\-\\d+.*")
            || cmdTokens.get(0).matches("^\\+\\d+.*")) { // No cmd line args either f/b ops or ints
          String newRegion =
              Utils.parseConsoleInput(cmdTokens, proc.getGenomicCoordsHistory().current()).trim();
          GenomicCoords gc = (GenomicCoords) proc.getGenomicCoordsHistory().current().clone();
          this.repositionGenomicCoords(gc, newRegion, terminalWidth);
          proc.getGenomicCoordsHistory().add(gc);

        } else if (cmdTokens.get(0).matches("(\\[+|\\]+)\\d*")) {
          this.interactiveInputExitCode = this.bracketsCommand(cmdTokens, proc, terminalWidth);
        } else if (cmdTokens.get(0).matches("^\\d+.*") || cmdTokens.get(0).matches("^\\.\\d+.*")) {
          String newRegion;
          try {
            newRegion =
                this.gotoOnCurrentChrom(cmdTokens, proc.getGenomicCoordsHistory().current());
          } catch (IndexOutOfBoundsException e) {
            System.err.append("Column coordinates must be >= 1 and <= the screen width");
            throw new InvalidCommandLineException();
          }
          GenomicCoords gc = (GenomicCoords) proc.getGenomicCoordsHistory().current().clone();
          this.repositionGenomicCoords(gc, newRegion, terminalWidth);
          proc.getGenomicCoordsHistory().add(gc);

        } else if (cmdTokens.get(0).equals("goto") || cmdTokens.get(0).startsWith(":")) {
          this.interactiveInputExitCode = this.goTo(cmdTokens, proc, terminalWidth);

        } else if (cmdTokens.get(0).equals("nextChrom")) {
          cmdTokens.remove(0);
          this.interactiveInputExitCode = this.nextChrom(cmdTokens, proc);

        } else if (cmdTokens.get(0).equals("p")) {
          proc.getGenomicCoordsHistory().previous();

        } else if (cmdTokens.get(0).equals("n")) {
          proc.getGenomicCoordsHistory().next();

        } else if (cmdTokens.get(0).equals("zo")) {
          int nz = Utils.parseZoom(Joiner.on(" ").join(cmdTokens), 1);
          GenomicCoords gc = (GenomicCoords) proc.getGenomicCoordsHistory().current().clone();
          gc.setTerminalWidth(terminalWidth);
          for (int i = 0; i < nz; i++) {
            gc.zoomOut();
          }
          proc.getGenomicCoordsHistory().add(gc);

        } else if (cmdTokens.get(0).equals("zi")) {
          int nz = Utils.parseZoom(Joiner.on(" ").join(cmdTokens), 1);
          GenomicCoords gc = (GenomicCoords) proc.getGenomicCoordsHistory().current().clone();
          gc.setTerminalWidth(terminalWidth);
          for (int i = 0; i < nz; i++) {
            gc.zoomIn();
          }
          proc.getGenomicCoordsHistory().add(gc);

        } else if (cmdTokens.get(0).equals("extend")) {
          if (cmdTokens.size() == 1) {
            System.err.println(
                Utils.padEndMultiLine("Expected at least one argument.", proc.getWindowSize()));
          }
          GenomicCoords gc = (GenomicCoords) proc.getGenomicCoordsHistory().current().clone();
          gc.setTerminalWidth(terminalWidth);
          gc.cmdInputExtend(cmdTokens);
          proc.getGenomicCoordsHistory().add(gc);

        } else if (cmdTokens.get(0).equals("trim")) {
          GenomicCoords gc = proc.getTrackSet().trimCoordsForTrack(cmdTokens);
          proc.getGenomicCoordsHistory().add(gc);

        } else if (cmdTokens.get(0).equals("l")) {
          GenomicCoords gc = (GenomicCoords) proc.getGenomicCoordsHistory().current().clone();
          gc.setTerminalWidth(terminalWidth);
          gc.left();
          proc.getGenomicCoordsHistory().add(gc);

        } else if (cmdTokens.get(0).equals("r")) {
          GenomicCoords gc = (GenomicCoords) proc.getGenomicCoordsHistory().current().clone();
          gc.setTerminalWidth(terminalWidth);
          gc.right();
          proc.getGenomicCoordsHistory().add(gc);

        } else if (cmdTokens.get(0).equals("sessionOpen")) {
          List<String> args = new ArrayList<>(cmdTokens);
          args.remove(0);
          this.openSession(args, proc);

        } else if (cmdTokens.get(0).equals("sessionSave")) {
          List<String> args = new ArrayList<>(cmdTokens);
          args.remove(0);
          this.saveSession(args, proc);
          this.setInteractiveInputExitCode(ExitCode.CLEAN_NO_FLUSH);

        } else if (cmdTokens.get(0).equals("sessionList")) {
          List<String> args = new ArrayList<>(cmdTokens);
          args.remove(0);
          this.setInteractiveInputExitCode(this.listSessions(args, proc));

        } else if (cmdTokens.get(0).equals("sessionDelete")) {
          List<String> args = new ArrayList<>(cmdTokens);
          args.remove(0);
          this.setInteractiveInputExitCode(this.deleteSession(args, proc));

        } else if (cmdTokens.get(0).equals("setGenome")) {
          this.setGenome(cmdTokens, proc);

          // * These commands change the Tracks but do not touch the GenomicCoordinates.
        } else if (cmdTokens.get(0).equals("dataCol")) {
          try {
            proc.getTrackSet().setDataColForRegex(cmdTokens);
          } catch (Exception e) {
            String msg =
                Utils.padEndMultiLine(
                    "Error processing "
                        + cmdTokens
                        + ". Perhaps a non-numeric column was selected?",
                    proc.getWindowSize());
            System.err.println(msg);
            this.interactiveInputExitCode = ExitCode.ERROR;
            continue;
          }

        } else if (cmdTokens.get(0).equals("ylim")) {
          proc.getTrackSet().setTrackYlimitsForRegex(cmdTokens);

        } else if (cmdTokens.get(0).equals("trackHeight")) {
          proc.getTrackSet().setTrackHeightForRegex(cmdTokens);

        } else if (cmdTokens.get(0).equals("colourTrack")) {
          proc.getTrackSet().setTrackColourForRegex(cmdTokens);

        } else if (cmdTokens.get(0).equals("bedToBedgraph")) {
          proc.getTrackSet().setTrackFormatForRegex(cmdTokens);

        } else if (cmdTokens.get(0).equals("featureColour")) {
          proc.getTrackSet().setFeatureColourForRegex(cmdTokens);

        } else if (cmdTokens.get(0).equals("hideTitle")) {
          proc.getTrackSet().setHideTitleForRegex(cmdTokens);

        } else if (cmdTokens.get(0).equals("addHeader")) {
          this.interactiveInputExitCode = proc.getTrackSet().addHeader(cmdTokens, terminalWidth);

        } else if (cmdTokens.get(0).equals(Command.BSseq.getCmdDescr())) {
          if (proc.getGenomicCoordsHistory().current().getFastaFile() == null) {
            String msg =
                Utils.padEndMultiLine(
                    "Cannot set BSseq mode without reference sequence", proc.getWindowSize());
            System.err.println(msg);
            this.interactiveInputExitCode = ExitCode.ERROR;
            continue;
          }
          proc.getTrackSet().setBisulfiteModeForRegex(cmdTokens);

        } else if (cmdTokens.get(0).equals("squash")
            || cmdTokens.get(0).equals(Command.featureDisplayMode.toString())) {
          proc.getTrackSet().setFeatureDisplayModeForRegex(cmdTokens);

        } else if (cmdTokens.get(0).equals("gap")) {
          proc.getTrackSet().setFeatureGapForRegex(cmdTokens);

        } else if (cmdTokens.get(0).equals("readsAsPairs")) {
          proc.getTrackSet().setReadsAsPairsForRegex(cmdTokens);
        } else if (cmdTokens.get(0).equals("nameForFeatures")) {
          proc.getTrackSet().setNameAttribute(cmdTokens);
        } else if (cmdTokens.get(0).equals("open") || cmdTokens.get(0).equals("addTracks")) {
          cmdTokens.remove(0);
          this.addTracks(cmdTokens, proc);
        } else if (cmdTokens.get(0).equals("reload")) {
          proc.getTrackSet().reload(cmdTokens);
        } else if (cmdTokens.get(0).equals("dropTracks")) {
          if (cmdTokens.size() <= 1) {
            System.err.println(
                Utils.padEndMultiLine(
                    "List one or more tracks to drop or `dropTracks -h` for help.",
                    proc.getWindowSize()));
            this.interactiveInputExitCode = ExitCode.ERROR;
            continue;
          }
          messages.add(proc.getTrackSet().dropTracksWithRegex(cmdTokens));

        } else if (cmdTokens.get(0).equals("orderTracks")) {
          cmdTokens.remove(0);
          proc.getTrackSet().orderTracks(cmdTokens);

        } else if (cmdTokens.get(0).equals("editNames")) {
          messages.add(proc.getTrackSet().editNamesForRegex(cmdTokens));

        } else if (cmdTokens.get(0).equals(Command.print.toString())) {
          proc.getTrackSet().setPrintModeAndPrintFeaturesForRegex(cmdTokens);

        } else if (cmdTokens.get(0).equals("grep")) {
          proc.getTrackSet().setFilterForTrackIntervalFeature(cmdTokens);

        } else if (cmdTokens.get(0).equals("awk")) {
          proc.getTrackSet().setAwkForTrack(cmdTokens);

        } else if (cmdTokens.get(0).equals("filterVariantReads")) {
          proc.getTrackSet().setFilterVariantReads(cmdTokens);

        } else if (cmdTokens.get(0).equals(Command.rpm.getCmdDescr())) {
          proc.getTrackSet().setRpmForRegex(cmdTokens);

        } else if (cmdTokens.get(0).equals("samtools")) {
          proc.getTrackSet().setSamFilterForRegex(cmdTokens);

        } else if (cmdTokens.get(0).equals("genotype")) {
          proc.getTrackSet().setGenotypeMatrix(cmdTokens);

          // * These commands change both the Tracks and the GenomicCoordinates
        } else if (cmdTokens.get(0).equals("next")) {
          this.next(cmdTokens, proc);

        } else if (cmdTokens.get(0).equals("find")) {

          boolean all = Utils.argListContainsFlag(cmdTokens, "-all");
          boolean fixedPattern = Utils.argListContainsFlag(cmdTokens, "-F");
          boolean caseIns = Utils.argListContainsFlag(cmdTokens, "-c");

          if (cmdTokens.size() < 2) {
            System.err.println(
                Utils.padEndMultiLine(
                    "Error in find command. Expected at least 1 argument got: " + cmdTokens,
                    proc.getWindowSize()));
            this.interactiveInputExitCode = ExitCode.ERROR;
            continue;
          }
          if (cmdTokens.size() == 2) {
            cmdTokens.add(""); // If track arg is missing use this placeholder.
          }
          GenomicCoords gc = (GenomicCoords) proc.getGenomicCoordsHistory().current().clone();
          gc.setTerminalWidth(terminalWidth);

          int flag = 0;
          if (fixedPattern) {
            flag |= Pattern.LITERAL;
          }
          if (!caseIns) {
            flag |= Pattern.CASE_INSENSITIVE;
          }
          Pattern pattern;
          try {
            pattern = Pattern.compile(cmdTokens.get(1), flag);
          } catch (PatternSyntaxException e) {
            System.err.println("Invalid regex");
            throw new InvalidCommandLineException();
          }
          GenomicCoords nextGc =
              proc.getTrackSet().findNextMatchOnTrack(pattern, cmdTokens.get(2), gc, all);
          if (nextGc.equalCoords(gc)) {
            System.err.println(
                "No match found outside of this window for query '" + cmdTokens.get(1) + "'");
            this.interactiveInputExitCode = ExitCode.CLEAN_NO_FLUSH;
          } else {
            proc.getGenomicCoordsHistory().add(nextGc);
          }

        } else if (cmdTokens.get(0).equals("seqRegex")) {
          try {
            proc.getTrackSet()
                .setRegexForTrackSeqRegex(cmdTokens, proc.getGenomicCoordsHistory().current());
          } catch (InvalidCommandLineException e) {
            System.err.println(
                Utils.padEndMultiLine(
                    "Cannot find regex in sequence without fasta reference!",
                    proc.getWindowSize()));
            this.interactiveInputExitCode = ExitCode.ERROR;
            continue;
          }
        } else if (cmdTokens.get(0).equals("translate")) {
          List<String> args = new ArrayList<>(cmdTokens);
          args.remove(0);
          this.interactiveInputExitCode = this.translate(args, proc);
        } else if (cmdTokens.get(0).equals("bookmark")) {
          messages.add(
              proc.getTrackSet().bookmark(proc.getGenomicCoordsHistory().current(), cmdTokens));

        } else {
          System.err.println(
              Utils.padEndMultiLine(
                  "Unrecognized command: " + cmdTokens.get(0), proc.getWindowSize()));
          String suggestions =
              Joiner.on(" or ")
                  .join(Utils.suggestCommand(cmdTokens.get(0).trim(), CommandList.cmds()));
          if (!suggestions.isEmpty()) {
            System.err.println(
                Utils.padEndMultiLine("Maybe you mean " + suggestions + "?", proc.getWindowSize()));
          }
          this.interactiveInputExitCode = ExitCode.ERROR;
        }
      } catch (ArgumentParserException e) {
        this.interactiveInputExitCode = ExitCode.ERROR;

      } catch (Exception e) { // You shouldn't catch anything! Be more specific.
        System.err.println(
            Utils.padEndMultiLine("\nError processing input: " + cmdTokens, proc.getWindowSize()));
        System.err.println(
            Utils.padEndMultiLine(
                "For help on command \"cmd\" execute 'cmd -h' or '-h' for list of commands.\n",
                proc.getWindowSize()));
        this.interactiveInputExitCode = ExitCode.ERROR;
        if (debug == 1) {
          e.printStackTrace();
        } else if (debug == 2) {
          e.printStackTrace();
          System.exit(1);
        }
      } // END PARSING ONE COMMAND

      if (this.interactiveInputExitCode.equals(ExitCode.CLEAN)
          || this.interactiveInputExitCode.equals(ExitCode.CLEAN_NO_FLUSH)) {
        // Command has been parsed ok. Let's see if we can execute it without exceptions.
        try {
          if (this.interactiveInputExitCode.equals(ExitCode.CLEAN)) {
            console.clearScreen();
            console.flush();
            proc.iterateTracks();
          } else {
            //
          }

        } catch (InvalidGenomicCoordsException e) {
          String newRegion =
              Main.initRegion(proc.getTrackSet().getFilenameList(), null, null, debug);
          GenomicCoords gc = (GenomicCoords) proc.getGenomicCoordsHistory().current().clone();
          this.repositionGenomicCoords(gc, newRegion, terminalWidth);
          proc.getGenomicCoordsHistory().add(gc);
          System.err.println(
              Utils.padEndMultiLine(
                  "Invalid genomic coordinates found. Resetting to " + newRegion,
                  proc.getWindowSize()));
          if (debug > 0) {
            e.printStackTrace();
          }

        } catch (Exception e) {
          System.err.println(
              Utils.padEndMultiLine(
                  "Error processing tracks with input " + cmdTokens, proc.getWindowSize()));
          this.interactiveInputExitCode = ExitCode.ERROR;
          if (debug > 0) {
            e.printStackTrace();
          }
        }
      }
      if (this.interactiveInputExitCode.equals(ExitCode.ERROR)) {
        // If something goes wrong or help is invoked, stop executing commands and restart asking
        // for input
        // Unless we are in non-interactive mode
        if (nonInteractive) {
          System.exit(1);
        }
        break;
      }
    } // END OF LOOP THROUGH CHAIN OF INPUT COMMANDS
    if (!this.messages.isEmpty()) {
      System.err.println(
          Utils.padEndMultiLine(Joiner.on("\n").join(this.messages), proc.getWindowSize()));
    }
    this.messages = new ArrayList<>();
    return proc;
  }

  private ExitCode bracketsCommand(List<String> cmdTokens, TrackProcessor proc, int terminalWidth)
      throws InvalidGenomicCoordsException, IOException, InvalidCommandLineException {
    int times = 0;
    try {
      times = this.countBrackets(cmdTokens);
    } catch (InvalidCommandLineException e) {
      return ExitCode.ERROR;
    }
    String newRegion =
        this.moveWindowByColumns(proc.getGenomicCoordsHistory().current(), times);
    GenomicCoords gc = (GenomicCoords) proc.getGenomicCoordsHistory().current().clone();
    this.repositionGenomicCoords(gc, newRegion, terminalWidth);
    proc.getGenomicCoordsHistory().add(gc);
    return ExitCode.CLEAN;
  }

  private ExitCode goTo(List<String> cmdTokens, TrackProcessor proc, int terminalWidth)
      throws InvalidGenomicCoordsException, IOException {
    String tmpRegion = Joiner.on(" ").join(cmdTokens).replaceFirst("goto|:", "").trim();
    List<String> reg;
    try {
      reg = Utils.parseStringCoordsToList(tmpRegion, 1, null);
    } catch (InvalidGenomicCoordsException e) {
      System.err.println(e.getMessage());
      return ExitCode.ERROR;
    }
    if (!proc.getGenomicCoordsHistory().current().isChromInSequenceDictionary(reg.get(0))) {
      System.err.println("Cannot find chromosome '" + reg.get(0) + "' in sequence dictionary");
      return ExitCode.ERROR;
    }
    int from = Integer.parseInt(reg.get(1));
    int to;
    if (reg.get(2) == null) {
      to = from + terminalWidth;
    } else {
      to = Integer.parseInt(reg.get(2));
    }
    if (from == to) {
      to = to + terminalWidth;
    }
    String newRegion = Utils.coordinatesToString(reg.get(0), from, to);
    GenomicCoords gc = (GenomicCoords) proc.getGenomicCoordsHistory().current().clone();
    this.repositionGenomicCoords(gc, newRegion, terminalWidth);
    proc.getGenomicCoordsHistory().add(gc);
    this.interactiveInputExitCode = ExitCode.CLEAN;
    return ExitCode.CLEAN;
  }

  private void repositionGenomicCoords(GenomicCoords gc, String newRegion, int terminalWidth)
      throws InvalidGenomicCoordsException, IOException {
    List<String> reg = Utils.parseStringCoordsToList(newRegion);
    String chrom = reg.get(0);
    int from = Integer.parseInt(reg.get(1));
    int to = Integer.parseInt(reg.get(2));
    int span = to - from + 1;
    if (this.samSeqDict != null ) {
      int seqLen = this.samSeqDict.getSequence(chrom).getSequenceLength();
      if (to > seqLen) {
        to = seqLen;
        from = to - span +  1;
        from = Math.max(from, 1);
      }
    }
    GenomicCoords tmp = new GenomicCoords(Utils.coordinatesToString(chrom, from, to), terminalWidth, this.samSeqDict, this.fasta);
    gc.setChrom(tmp.getChrom());
    gc.setTo(tmp.getTo());
    gc.setFrom(tmp.getFrom());
    gc.setTerminalWidth(terminalWidth);
    gc.setSamSeqDict(this.samSeqDict);
    gc.setFastaFile(this.fasta);
  }

  private void repositionGenomicCoords(GenomicCoords gc)
      throws InvalidGenomicCoordsException, IOException {
    this.repositionGenomicCoords(gc, gc.toStringRegion(), gc.getTerminalWidth());
  }

  private ExitCode translate(List<String> args, TrackProcessor proc)
      throws InvalidCommandLineException, InvalidGenomicCoordsException, IOException {
    GenomicSequence gs = proc.getGenomicCoordsHistory().current().getGenomicSequence();
    if (gs.getSequence() == null) {
      System.err.println(
          Utils.padEndMultiLine("Unable to translate without a reference sequence", proc.getWindowSize()));
      return ExitCode.ERROR;
    }
    String frame;
    if (args.isEmpty()){
      // With no arguments: toggle on/off
      if (gs.getFrames().isEmpty()) {
        frame = "all";
      } else {
        frame = "none";
      }
    } else {
      frame = Utils.getArgForParam(args, "-frame", "all");
    }
    String codon = Utils.getArgForParam(args, "-codon", "all");

    String currentCode = gs.getGeneticCode();
    if (currentCode == null) {
      currentCode = "universal";
    }
    String geneticCode = Utils.getArgForParam(args, "-geneticCode", currentCode);

    IUPACParser.IUPACTable table = IUPACParser.getInstance().getTable(geneticCode.toUpperCase());
    if (table == null) {
      List<String> tables = new ArrayList<>();
      for (IUPACParser.IUPACTable x : IUPACParser.getInstance().getTables()) {
        tables.add("  " + x.getName());
      }
      System.err.println(
          Utils.padEndMultiLine(
              "Invalid translation table: '"
                  + geneticCode.toUpperCase()
                  + "'\n"
                  + "Valid tables are:\n"
                  + Joiner.on('\n').join(tables),
              proc.getWindowSize()));
      return ExitCode.ERROR;
    }

    PrintCodon printCodon;
    try {
      printCodon = PrintCodon.valueOf(codon.toUpperCase());
    } catch (IllegalArgumentException ex) {
      System.err.println("Invalid option for codon: '" + codon + "'");
      return ExitCode.ERROR;
    }

    if (frame.equalsIgnoreCase("none")) {
      gs.setFrames(new Frame[] {});
    } else if (frame.equalsIgnoreCase("all")) {
      gs.setFrames(Frame.getAllFrames());
    } else if (frame.equalsIgnoreCase("forward")) {
      gs.setFrames(Frame.getForwardFrames());
    } else if (frame.equalsIgnoreCase("reverse")) {
      gs.setFrames(Frame.getReverseFrames());
    } else {
      System.err.println("Invalid option for frame: '" + frame + "'");
      return ExitCode.ERROR;
    }
    gs.setGeneticCode(geneticCode);
    gs.setPrintCodon(printCodon);
    this.repositionGenomicCoords(proc.getGenomicCoordsHistory().current());
    return ExitCode.CLEAN;
  }

  private ExitCode nextChrom(List<String> cmdTokens, TrackProcessor proc)
      throws IOException, NumberFormatException, InvalidCommandLineException {
    int minSize = Integer.parseInt(Utils.getArgForParam(cmdTokens, "-min", "-1"));
    int maxSize = Integer.parseInt(Utils.getArgForParam(cmdTokens, "-max", "-1"));
    String so = Utils.getArgForParam(cmdTokens, "-s", "s");
    ContigOrder sortOrder;
    if (so.equals("u")) {
      sortOrder = ContigOrder.UNSORTED;
    } else if (so.equals("s")) {
      sortOrder = ContigOrder.SIZE_ASC;
    } else if (so.equals("S")) {
      sortOrder = ContigOrder.SIZE_DESC;
    } else {
      System.err.println("Invalid option for sort order: " + so);
      return ExitCode.CLEAN_NO_FLUSH;
    }

    String regex = ".*";
    if (cmdTokens.size() == 1) {
      regex = cmdTokens.get(0);
    } else if (cmdTokens.size() > 1) {
      System.err.println("At most only one positional arguments is allowed. Got: " + cmdTokens);
      return ExitCode.CLEAN_NO_FLUSH;
    }

    try {
      GenomicCoords gc = (GenomicCoords) proc.getGenomicCoordsHistory().current().clone();
      gc.nextChrom(proc.getTrackSet().getKnownContigs(), minSize, maxSize, regex, sortOrder);
      proc.getGenomicCoordsHistory().add(gc);
      return ExitCode.CLEAN;
    } catch (InvalidGenomicCoordsException e) {
      System.err.println(e.getMessage());
      return ExitCode.CLEAN_NO_FLUSH;
    }
  }

  private int countBrackets(List<String> cmdTokens) throws InvalidCommandLineException {

    String xtimes = cmdTokens.get(0).replaceAll("\\[|\\]", "");
    if (!xtimes.isEmpty() && cmdTokens.size() > 1) {
      throw new InvalidCommandLineException();
    }
    int times = 1;
    if (!xtimes.isEmpty()) {
      times = Integer.valueOf(xtimes);
    } else if (cmdTokens.size() > 1) {
      try {
        times = Integer.valueOf(cmdTokens.get(1));
      } catch (NumberFormatException e) {
        System.err.append(cmdTokens.get(1) + " cannot be converted to integer");
        throw new InvalidCommandLineException();
      }
    }
    String brackets = cmdTokens.get(0).replaceAll("\\d", "");
    times = times * brackets.length();

    if (cmdTokens.get(0).startsWith("]")) {
      times *= -1;
    }
    return times;
  }

  /**
   * Return a string where the current genomic coordinates are moved forward or backwards "times"
   * screen columns. Move backwards if times is negative.
   *
   * @throws IOException
   * @throws InvalidGenomicCoordsException
   * @throws InvalidCommandLineException
   */
  private String moveWindowByColumns(GenomicCoords gc, int times)
      throws InvalidGenomicCoordsException, IOException, InvalidCommandLineException {
    String x = String.valueOf(Math.round(gc.getBpPerScreenColumn() * times));
    if (Integer.valueOf(x) >= 0) {
      x = "+" + x;
    }
    List<String> tokens = new ArrayList<String>();
    tokens.add(String.valueOf(x));
    return Utils.parseConsoleInput(tokens, gc);
  }

  private ExitCode explainSamFlag(List<String> cmdTokens, TrackProcessor proc)
      throws InvalidCommandLineException, IOException {
    List<String> args = new ArrayList<String>(cmdTokens);
    args.remove(0);
    if (args.size() == 0) {
      System.err.println("One argument is required.");
      throw new InvalidCommandLineException();
    }
    List<Integer> flagsToDecode = new ArrayList<Integer>();
    for (String x : args) {
      try {
        flagsToDecode.add(Integer.valueOf(x));
      } catch (NumberFormatException e) {
        System.err.println(
            Utils.padEndMultiLine("Expected an integer, got " + x, Utils.getTerminalWidth()));
        return ExitCode.ERROR;
      }
    }
    final Map<Integer, String> FLAGS = new LinkedHashMap<Integer, String>();

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

    String LEFT_PAD = " ";
    List<String> table = new ArrayList<String>();
    table.add(LEFT_PAD + Joiner.on("\t").join(flagsToDecode));

    System.out.println(Utils.padEndMultiLine("", Utils.getTerminalWidth()));
    for (Integer i : FLAGS.keySet()) {
      String line = LEFT_PAD;
      for (int f : flagsToDecode) {
        if ((i & f) == i) {
          line += "X";
        } else {
          line += ".";
        }
        line += "\t";
      }
      line += FLAGS.get(i);
      table.add(line);
    }
    List<String> xtable = Utils.tabulateList(table, Utils.getTerminalWidth(), " ");
    System.out.println(
        Utils.padEndMultiLine(Joiner.on("\n").join(xtable), Utils.getTerminalWidth()));
    System.out.println(Utils.padEndMultiLine("", Utils.getTerminalWidth()));
    return ExitCode.CLEAN_NO_FLUSH;
  }

  /** Get the items (files) corresponding to the indexes. Errors are silently ignored. */
  private List<String> openFilesFromIndexes(
      LinkedHashSet<String> openedFiles, List<String> indexes) {
    List<String> files = new ArrayList<String>();
    List<Integer> idxs = new ArrayList<Integer>();
    for (String x : indexes) {
      try {
        idxs.add(Integer.parseInt(x));
      } catch (NumberFormatException e) {
        // Return empty list
        return files;
        //
      }
    }
    for (int i : idxs) {
      try {
        String x = Lists.reverse(new ArrayList<String>(openedFiles)).get(i - 1);
        files.add(x);
      } catch (Exception e) {
        //
      }
    }
    return files;
  }

  /**
   * Parse the given list of options and move to new coordinates. Visibility set to protected only
   * for testing.
   */
  protected String gotoOnCurrentChrom(List<String> cmdTokens, GenomicCoords gc)
      throws InvalidGenomicCoordsException, IOException {
    List<String> args = new ArrayList<String>(cmdTokens);
    int regFrom;
    int regTo;

    boolean center = false;
    if (args.get(0).endsWith("c")) {
      center = true;
      args.set(0, args.get(0).replaceAll("c$", ""));
    }
    if (args.size() > 1 && args.get(1).equals("c")) {
      center = true;
      args.remove(1);
    }

    if (Float.valueOf(args.get(0)) < 1) {
      // Switch to screen percent coordinates
      float pctFrom = Float.parseFloat(args.get(0));
      // Which is the screen column matching this percent?
      int screenIdx = (int) Math.rint(pctFrom * gc.getUserWindowSize());
      // Which is the genomic coordinate corresponding to this screen index?
      regFrom = (int) Math.rint(gc.getMapping().get(screenIdx));

      // Same for end position. Accounting for possibility that only one pct value is given
      if (args.size() > 1 && !center) {
        float pctTo = Float.parseFloat(args.get(1));
        if (pctTo > 1) {
          System.err.println(
              "Coordinate interval given as percent of screen width must be between 0 and 1. Got "
                  + args);
          throw new InvalidGenomicCoordsException();
        }
        screenIdx = (int) Math.rint(pctTo * gc.getUserWindowSize());
        if (screenIdx >= gc.getMapping().size()) {
          screenIdx = gc.getMapping().size() - 1;
        }
        regTo = (int) Math.rint(gc.getMapping().get(screenIdx));

      } else if (center) {
        regFrom = regFrom - (int) Math.floor((double) gc.getGenomicWindowSize() / 2) - 1;
        regFrom = Math.max(regFrom, 1);
        regTo = regFrom + gc.getGenomicWindowSize() - 1;
      } else {
        regTo = regFrom + gc.getGenomicWindowSize() - 1;
      }
      if (regTo < regFrom) {
        System.err.println("Invalid coordinates: end < start for argument(s): " + args);
        throw new InvalidGenomicCoordsException();
      }
      // return gc.getChrom() + ":" + regFrom + "-" + regTo;
    } else if (args.size() == 1 && !center) {
      regFrom = Integer.parseInt(args.get(0));
      regTo = regFrom + gc.getUserWindowSize() - 1;
    } else if (!center) {
      regFrom = Integer.parseInt(args.get(0));
      regTo = Integer.parseInt(args.get(args.size() - 1));
    } else if (center) {
      regFrom = Integer.parseInt(args.get(0)) - (int) Math.floor((double) gc.getGenomicWindowSize() / 2);
      regFrom = Math.max(regFrom, 1);
      regTo = regFrom + gc.getGenomicWindowSize() - 1;
    } else {
      throw new IllegalArgumentException();
    }
    return gc.getChrom() + ":" + regFrom + "-" + regTo;
  }

  private ExitCode setConfigOpt(List<String> cmdTokens, TrackProcessor proc)
      throws IOException,
          InvalidConfigException,
          InvalidCommandLineException,
          InvalidColourException, InvalidGenomicCoordsException {
    List<String> args = new ArrayList<String>(cmdTokens);
    args.remove(0);
    if (args.isEmpty()) {
      throw new InvalidCommandLineException();
    }
    if (args.size() == 1) {
      ConfigKey key = ConfigKey.getConfigKeyFromShort(args.get(0));
      if (ConfigKey.booleanKeys().contains(key)) {
        // If configkey is a type boolean, just flip the boolean
        key = ConfigKey.valueOf(key.toString());
        boolean value = !Utils.asBoolean(Config.get(key));
        Config.set(key, String.valueOf(value));
      } else {
        // configKey is expected to be the name of a configuration file
        new Config(args.get(0));
      }
    } else {
      ConfigKey key = ConfigKey.getConfigKeyFromShort(args.get(0));
      String value = args.get(1);
      try {
        Config.set(key, value);
      } catch (Exception e) {
        System.err.println(
            Utils.padEndMultiLine("Unable to set configuration", proc.getWindowSize()));
        return ExitCode.ERROR;
      }
    }
    return ExitCode.CLEAN;
  }

  private void saveSession(List<String> args, TrackProcessor proc)
      throws InvalidCommandLineException, IOException, SessionException, InvalidGenomicCoordsException {
    File sessionYamlFile;
    String sf = Utils.getArgForParam(args, "-f", null);
    if (sf == null && this.sessionHandler == null) {
      sessionYamlFile = SessionHandler.DEFAULT_SESSION_FILE;
    } else if (sf == null) {
      sessionYamlFile = sessionHandler.getSessionFile();
    } else {
      sessionYamlFile = new File(sf);
    }

    if (!sessionYamlFile.getAbsoluteFile().getParentFile().exists()) {
      this.messages.add("Directory '" + sessionYamlFile.getParentFile() + "' does not exist");
      this.interactiveInputExitCode = ExitCode.ERROR;
      return;
    }
    if (!sessionYamlFile.exists()) {
      try {
        sessionYamlFile.createNewFile();
        sessionYamlFile.delete();
      } catch (IOException e) {
        this.messages.add("Cannot write to file '" + sessionYamlFile + "'");
        this.interactiveInputExitCode = ExitCode.ERROR;
        return;
      }
    } else if (!sessionYamlFile.canWrite()) {
      this.messages.add("Cannot write to file '" + sessionYamlFile + "'");
      this.interactiveInputExitCode = ExitCode.ERROR;
      return;
    }
    String sessionName = this.getCurrentSessionName();
    if (args.isEmpty() && sessionName == null) {
      this.messages.add("Please provide a session name. See `session -h` for help");
      this.interactiveInputExitCode = ExitCode.ERROR;
      return;
    } else if (args.size() > 1) {
      this.messages.add("Too many arguments. See `session -h` for help");
      this.interactiveInputExitCode = ExitCode.ERROR;
      return;
    } else if (args.size() == 1) {
      sessionName = args.get(0);
    }
    SessionHandler.saveAs(sessionYamlFile, proc.toSession(), sessionName);
    this.setCurrentSessionName(sessionName);
    this.sessionHandler = new SessionHandler(sessionYamlFile);
    this.messages.add(
        "Session '" + sessionName + "' saved to '" + sessionYamlFile.getAbsolutePath() + "'");
  }

  private ExitCode listSessions(List<String> args, TrackProcessor proc)
      throws InvalidGenomicCoordsException, IOException, InvalidCommandLineException {
    File sessionYamlFile =
        new File(
            Utils.getArgForParam(
                args, "-f", SessionHandler.DEFAULT_SESSION_FILE.getAbsolutePath()));
    SessionHandler sh;
    if (this.sessionHandler == null) {
      try {
        sh = new SessionHandler(sessionYamlFile);
      } catch (Exception e) {
        this.messages.add(
            "Failed to process session file '" + sessionYamlFile + "':\n" + e.getMessage());
        return ExitCode.ERROR;
      }
    } else {
      sh = this.sessionHandler;
    }
    int nsessions = Integer.parseInt(Utils.getArgForParam(args, "-n", "10"));
    String sessionString = sh.print(nsessions, true);
    sessionString += "\nCurrent session: ";
    sessionString += this.getCurrentSessionName() == null ? "n/a" : this.getCurrentSessionName();
    System.err.println(Utils.padEndMultiLine(sessionString, proc.getWindowSize()));
    return ExitCode.CLEAN_NO_FLUSH;
  }

  private ExitCode deleteSession(List<String> args, TrackProcessor proc)
      throws InvalidGenomicCoordsException, IOException, InvalidCommandLineException {
    File sessionYamlFile =
        new File(
            Utils.getArgForParam(
                args, "-f", SessionHandler.DEFAULT_SESSION_FILE.getAbsolutePath()));
    SessionHandler sh;
    if (this.sessionHandler == null) {
      try {
        sh = new SessionHandler(sessionYamlFile);
      } catch (Exception e) {
        this.messages.add(
            "Failed to process session file '" + sessionYamlFile + "':\n" + e.getMessage());
        return ExitCode.ERROR;
      }
    } else {
      sh = this.sessionHandler;
    }
    if (args.size() == 0) {
      this.messages.add("Please provide the name of the session to delete");
      return ExitCode.ERROR;
    }
    if (args.size() > 1) {
      this.messages.add("Please provide only one session to delete");
      return ExitCode.ERROR;
    }
    String sessionName = args.get(0);
    boolean found = sh.hasSessionName(sessionName);
    if (!found) {
      this.messages.add(
          "No session with name '" + sessionName + "' found in file '" + sessionYamlFile + "'");
      return ExitCode.ERROR;
    }
    boolean deleted = sh.deleteSession(sessionName);
    if (deleted) {
      System.err.println(
          Utils.padEndMultiLine("Session '" + sessionName + "' deleted", proc.getWindowSize()));
    } else {
      this.messages.add("Failed to delete session '" + sessionName + "'");
      return ExitCode.ERROR;
    }
    writeSessions(sh.getSessions(), sessionYamlFile);
    return ExitCode.CLEAN_NO_FLUSH;
  }

  private void openSession(List<String> args, TrackProcessor proc)
      throws IOException, InvalidCommandLineException {
    File sessionYamlFile =
        new File(
            Utils.getArgForParam(
                args, "-f", SessionHandler.DEFAULT_SESSION_FILE.getAbsolutePath()));
    String sessionNameOrIndex;
    if (args.isEmpty()) {
      sessionNameOrIndex = "1";
    } else if (args.size() == 1) {
      sessionNameOrIndex = args.get(0);
    } else {
      this.messages.add("Too many arguments. See `session -h` for details");
      this.interactiveInputExitCode = ExitCode.ERROR;
      return;
    }
    Session session;
    try {
      this.sessionHandler = new SessionHandler(sessionYamlFile);
      session = this.sessionHandler.get(sessionNameOrIndex);
      if (session.getGenome().samSeqDictSource != null
          && !new File(session.getGenome().samSeqDictSource).exists()) {
        this.messages.add(
            "Warning: Sequence dictionary file '"
                + session.getGenome().samSeqDictSource
                + "' does not exist");
      }
      if (session.getGenome().fastaFile != null
          && !new File(session.getGenome().fastaFile).exists()) {
        this.messages.add(
            "Warning: Fasta file '" + session.getGenome().fastaFile + "' does not exist");
      }
    } catch (SessionException e) {
      this.messages.add(e.getMessage());
      this.interactiveInputExitCode = ExitCode.ERROR;
      return;
    }
    try {
      GenomicCoords gc = session.getGenome().toGenomicCoords();
      proc.getGenomicCoordsHistory().setGenome(Collections.singletonList(gc.getFastaFile()));
      proc.getGenomicCoordsHistory().add(gc);
      proc.setTrackSet(session.toTrackSet());
      this.setCurrentSessionName(session.getSessionName());
    } catch (Exception e) {
      this.messages.add(
          "Unable to open session "
              + sessionNameOrIndex
              + ". Session is:\n"
              + session
              + "\n"
              + e.getMessage());
    }
    List<String> passed = new ArrayList<>();
    for (Track x : proc.getTrackSet().getTrackList()) {
      passed.add(x.getTrackTag());
    }
    for (String tag : session.getTracks().keySet()) {
      if (!passed.contains(tag)) {
        this.messages.add("Warning: Failed to add track '" + tag + "'");
      }
    }
  }

  private void setGenome(List<String> cmdTokens, TrackProcessor proc)
      throws InvalidGenomicCoordsException, IOException, InvalidCommandLineException {

    List<String> tokens = new ArrayList<String>(cmdTokens);
    tokens.remove(0);

    if (tokens.isEmpty()) {
      // Try to read fasta from history file
      ASCIIGenomeHistory ag = new ASCIIGenomeHistory();
      try {
        tokens.add(ag.getReference().get(0));
        System.err.println("Using " + ag.getReference().get(0));
      } catch (Exception e) {
        System.err.println("A previous reference file was not found.");
        throw new InvalidCommandLineException();
      }
    }

    GenomicCoords testSeqDict = new GenomicCoords("default", Utils.getTerminalWidth(), null, null);
    testSeqDict.setGenome(tokens, true);
    if (testSeqDict.getSamSeqDict() != null) {
      proc.getGenomicCoordsHistory().setGenome(tokens);
    } else {
      System.err.println(
          Utils.padEndMultiLine("Cannot set genome from " + tokens, Utils.getTerminalWidth()));
      this.interactiveInputExitCode = ExitCode.ERROR;
    }
  }

  /**
   * Execute arbitrary system command and print its output
   *
   * @param cmdInput: String, in contrast to other coomands, process the raw string, not the
   *     tokenized version so you don't mess up with single quotes inside the system command.
   */
  private void execSysCmd(String cmdInput, int userWindowSize) {

    String rawSysCmd = cmdInput.trim().replaceAll("^sys +", ""); // Remove command name
    boolean isLiteral = false;
    if (rawSysCmd.trim().startsWith("-L ")) {
      // Check if -L option is present and remove it if yes.
      isLiteral = true;
      rawSysCmd = rawSysCmd.trim().replaceAll("-L +", "");
    }
    if (rawSysCmd.isEmpty()) {
      System.err.println(
          Utils.padEndMultiLine(
              "Please provide a system comand to execute. Use `sys -h` for help.", userWindowSize));
      return;
    }

    List<String> tokens = new ArrayList<String>();

    if (isLiteral) {
      tokens.add(rawSysCmd);
    } else { // w/o -L option execute command as bash string
      tokens.add("bash");
      tokens.add("-c");
      tokens.add(rawSysCmd);
    }

    try {
      Process p = new ProcessBuilder().inheritIO().command(tokens).start();
      BufferedReader reader = new BufferedReader(new InputStreamReader(p.getInputStream()));

      String line = "";
      while ((line = reader.readLine()) != null) {
        System.err.println(Utils.padEndMultiLine(line, userWindowSize));
      }
      p.waitFor();

    } catch (Exception e) {
      System.err.println(Utils.padEndMultiLine(e.getMessage(), userWindowSize));
    }
  }

  private void posHistory(List<String> cmdInput, List<GenomicCoords> history, int userWindowSize)
      throws InvalidCommandLineException {

    List<String> args = new ArrayList<String>(cmdInput);
    args.remove(0); // Remove cmd name.

    int n = 10;
    if (args.contains("-n")) {
      try {
        n = Integer.parseInt(args.get(args.indexOf("-n") + 1));
        args.remove(args.get(args.indexOf("n") + 1));
        args.remove("-n");
      } catch (Exception e) {
        System.err.println(
            Utils.padEndMultiLine("Argument to -n parameter must be an integer", userWindowSize));
        throw new InvalidCommandLineException();
      }
    }

    int start = 0; // Start listing positions from this index
    if (history.size() > n && n > 0) {
      start = history.size() - n;
    }

    for (int i = start; i < history.size(); i++) {
      GenomicCoords xgc = history.get(i);
      System.err.println(Utils.padEndMultiLine(xgc.toString(), userWindowSize));
    }
  }

  /**
   * Move to next feature using parameters in cmdInput. First arg in cmdInput is command name
   * itself. The side effect is to modify the TrackProcessor obj to update the position.
   */
  private void next(List<String> cmdInput, TrackProcessor proc)
      throws InvalidGenomicCoordsException, IOException, InvalidCommandLineException {

    List<String> args = new ArrayList<String>(cmdInput);
    args.remove(0); // Remove command name

    int zo = 5;
    if (args.contains("-zo")) {
      try {
        zo = Integer.parseInt(args.get(args.indexOf("-zo") + 1));
        args.remove(args.get(args.indexOf("-zo") + 1));
        args.remove("-zo");
        if (zo < 0) {
          zo = 0;
        }
      } catch (Exception e) {
        System.err.println(
            Utils.padEndMultiLine(
                "Argument to -zo parameter must be an integer", proc.getWindowSize()));
        throw new InvalidCommandLineException();
      }
    }

    GenomicCoords gc = (GenomicCoords) proc.getGenomicCoordsHistory().current().clone();
    boolean start = false;
    if (args.contains("-start")) {
      start = true;
      args.remove("-start");
    }
    boolean center = false;
    if (args.contains("-c")) {
      center = true;
      args.remove("-c");
    }
    boolean getPrevious = false;
    if (args.contains("-back")) {
      getPrevious = true;
      args.remove("-back");
    }
    String trackId = "";
    if (!args.isEmpty()) {
      trackId = args.get(0);
    }
    if (start) {
      proc.getGenomicCoordsHistory()
          .add(proc.getTrackSet().goToNextFeatureOnFile(trackId, gc, -1.0, getPrevious));
    } else if (center) {
      proc.getGenomicCoordsHistory()
          .add(proc.getTrackSet().goToNextFeatureOnFile(trackId, gc, 0, getPrevious));
    } else {
      proc.getGenomicCoordsHistory()
          .add(proc.getTrackSet().goToNextFeatureOnFile(trackId, gc, zo, getPrevious));
    }
  }

  /**
   * Edit visualization setting in TrackProcessor as appropriate.
   *
   * @return
   * @throws InvalidCommandLineException
   * @throws IOException
   * @throws InvalidGenomicCoordsException
   */
  private ExitCode show(List<String> cmdTokens, TrackProcessor proc)
      throws InvalidCommandLineException, InvalidGenomicCoordsException, IOException {
    List<String> args = new ArrayList<String>(cmdTokens);
    args.remove(0);
    if (args.size() == 0) {
      System.err.println("At least one argument is required.");
      throw new InvalidCommandLineException();
    }
    // With .startsWith() we allow partial matching of input to argument. I.e. "ge" will be enough
    // to
    // recognize "genome"."genome"
    if ("genome".startsWith(args.get(0))) {
      ExitCode exitCode;
      try {
        exitCode = this.showGenome(cmdTokens, proc);
      } catch (InvalidGenomicCoordsException e) {
        System.err.println(e.getMessage());
        exitCode = ExitCode.ERROR;
      }
      return exitCode;

    } else if ("trackInfo".startsWith(args.get(0))) {
      String info = Utils.padEndMultiLine(proc.getTrackSet().showTrackInfo(), proc.getWindowSize());
      System.err.println(info);
      return ExitCode.CLEAN_NO_FLUSH;

    } else if ("gruler".startsWith(args.get(0))) {
      proc.setShowGruler(!proc.isShowGruler());
      return ExitCode.CLEAN;

    } else if ("pctRuler".startsWith(args.get(0))) {
      proc.setShowCruler(!proc.isShowCruler());
      return ExitCode.CLEAN;

    } else {
      System.err.println("Unrecognized option: " + args.get(0));
      throw new InvalidCommandLineException();
    }
  }

  /**
   * Print to screen the genome file
   *
   * @throws InvalidCommandLineException
   */
  private ExitCode showGenome(List<String> cmdTokens, TrackProcessor proc)
      throws InvalidGenomicCoordsException, IOException, InvalidCommandLineException {
    int maxLines;
    try {
      String n = Utils.getArgForParam(cmdTokens, "-n", "50");
      maxLines = Integer.valueOf(n);
    } catch (NumberFormatException e) {
      System.err.println("Argument to -n must be an integer");
      return ExitCode.CLEAN_NO_FLUSH;
    }
    if (maxLines < 0) {
      maxLines = Integer.MAX_VALUE;
    }
    String genome =
        proc.getGenomicCoordsHistory()
            .current()
            .printSequenceDictionary(
                proc.getTrackSet().getKnownContigs(),
                -1,
                -1,
                ".*",
                ContigOrder.SIZE_DESC,
                30,
                maxLines,
                proc.isNoFormat());
    if (genome != null && !genome.isEmpty()) {
      System.err.println(Utils.padEndMultiLine(genome, proc.getWindowSize()));
    }
    return ExitCode.CLEAN_NO_FLUSH;
  }

  private String cmdHistoryToString(List<String> cmdInput) throws InvalidCommandLineException {

    List<String> args = new ArrayList<String>(cmdInput);
    args.remove(0);

    String re = ".*";
    int nmax = Integer.MAX_VALUE;
    if (args.contains("-grep")) {
      re = args.get(args.indexOf("-grep") + 1);
    }
    if (args.contains("-n")) {
      nmax = Integer.parseInt(args.get(args.indexOf("-n") + 1));
    }

    Pattern pattern = Pattern.compile(re); // .matcher(x).find();
    List<String> cmd = new ArrayList<String>();
    int i = 1;
    for (Entry x : console.getHistory()) {
      if (pattern.matcher(x.value().toString()).find()) {
        cmd.add(i + ": \t" + x.value().toString());
      }
      i++;
    }

    if (cmd.size() > nmax) { // Trim list to
      cmd = cmd.subList(cmd.size() - nmax, cmd.size());
    }

    List<String> cmdTab = Utils.tabulateList(cmd, -1, " ");
    String tab = "";
    for (String x : cmdTab) {
      tab += (x + "\n");
    }
    return tab.trim();
  }

  private void addTracks(List<String> cmdTokens, TrackProcessor proc)
      throws IOException, InvalidGenomicCoordsException, InvalidCommandLineException {
    List<String> globbed = Utils.globFiles(cmdTokens);
    if (globbed.isEmpty()) {
      globbed = this.openFilesFromIndexes(proc.getTrackSet().getOpenedFiles(), cmdTokens);
    }

    if (globbed.isEmpty()) {
      String msg = Utils.padEndMultiLine(cmdTokens + ": No file found.", proc.getWindowSize());
      System.err.println(msg);
      this.interactiveInputExitCode = ExitCode.ERROR;

    } else {

      for (String sourceName : globbed) {
        String msg = Utils.padEndMultiLine("Adding: " + sourceName, proc.getWindowSize());
        System.err.println(msg);
        try {
          proc.getTrackSet()
              .addTrackFromSource(sourceName, proc.getGenomicCoordsHistory().current(), null);
        } catch (Exception e) {
          try {
            // It may be that you are in position that doesn't exist in the sequence
            // dictionary that
            // came with this new file. To recover, find an existing position, move there and
            // try to reload the
            // file. This fixes issue#23
            String region = Main.initRegion(globbed, null, null, debug);

            GenomicCoords gc = (GenomicCoords) proc.getGenomicCoordsHistory().current().clone();
            this.repositionGenomicCoords(gc, region, Utils.getTerminalWidth());
            proc.getGenomicCoordsHistory().add(gc);
            proc.getTrackSet()
                .addTrackFromSource(sourceName, proc.getGenomicCoordsHistory().current(), null);
          } catch (Exception x) {
            x.printStackTrace();
            msg = Utils.padEndMultiLine("Failed to add: " + sourceName, proc.getWindowSize());
            System.err.println(msg);
            this.interactiveInputExitCode = ExitCode.CLEAN_NO_FLUSH;
          }
        }

        if (proc.getGenomicCoordsHistory().current().getSamSeqDict() == null
            || proc.getGenomicCoordsHistory().current().getSamSeqDict().isEmpty()) {
          GenomicCoords testSeqDict =
              new GenomicCoords("default", Utils.getTerminalWidth(), null, null);
          List<String> candidateSourceGenome = new ArrayList<String>();
          candidateSourceGenome.add(sourceName);
          testSeqDict.setGenome(candidateSourceGenome, false);
          if (testSeqDict.getSamSeqDict() != null) {
            candidateSourceGenome.add(0, "cmd");
            proc.getGenomicCoordsHistory().setGenome(candidateSourceGenome);
          }
        }
      }
    }
  }

  public ExitCode getInteractiveInputExitCode() {
    return interactiveInputExitCode;
  }

  public void setInteractiveInputExitCode(ExitCode exitCode) {
    this.interactiveInputExitCode = exitCode;
  }

  public String getCurrentSessionName() {
    return currentSessionName;
  }

  public void setCurrentSessionName(String currentSessionName) {
    this.currentSessionName = currentSessionName;
  }
}
