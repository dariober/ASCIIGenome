package samTextViewer;

import colouring.Config;
import com.google.common.base.Splitter;
import exceptions.InvalidCommandLineException;
import exceptions.InvalidConfigException;
import exceptions.InvalidGenomicCoordsException;
import exceptions.InvalidRecordException;
import exceptions.SessionException;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileDescriptor;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.List;
import jline.console.ConsoleReader;
import org.junit.Test;
import session.SessionHandler;
import tracks.Track;
import tracks.TrackPileup;
import tracks.TrackSet;

import static org.junit.Assert.*;

public class InteractiveInputTest {

  public static class ProcessInput {
    public String stderr;
    public String stdout;
  }

  public ProcessInput processInput(InteractiveInput ip, String cmd, TrackProcessor p)
      throws SQLException,
          InvalidGenomicCoordsException,
          InvalidCommandLineException,
          IOException,
          ClassNotFoundException,
          InvalidRecordException {
    ByteArrayOutputStream err = new ByteArrayOutputStream();
    System.setErr(new PrintStream(err));
    ip.processInput(cmd, p);
    String errStr = err.toString();
    System.setErr(new PrintStream(new FileOutputStream(FileDescriptor.err)));
    ProcessInput pi = new ProcessInput();
    pi.stderr = errStr;
    return pi;
  }

  public static TrackProcessor gimmeTrackProcessor(String region, int terminalWidth)
      throws InvalidGenomicCoordsException,
          IOException,
          ClassNotFoundException,
          InvalidRecordException,
          SQLException {
    GenomicCoords gc = new GenomicCoords(region, terminalWidth, null, null);
    List<String> gf = new ArrayList<String>();
    gf.add("test_data/ds051.actb.bam");
    gc.setGenome(gf, false);
    GenomicCoordsHistory gch = new GenomicCoordsHistory();
    gch.add(gc);
    TrackSet trackSet = new TrackSet(new ArrayList<String>(), gc);
    TrackProcessor proc = new TrackProcessor(trackSet, gch);
    return proc;
  }

  @Test
  public void canListSessions()
      throws IOException,
          SQLException,
          InvalidGenomicCoordsException,
          ClassNotFoundException,
          InvalidRecordException,
          InvalidCommandLineException {
    TrackProcessor proc = gimmeTrackProcessor("chr7:1001-1800", 80);
    InteractiveInput ip = new InteractiveInput(new ConsoleReader(), 1);
    ProcessInput pi = processInput(ip, "session list -f test_data/session.yaml -n 1", proc);
    assertTrue(pi.stderr.contains("sessionName: no-fastafile"));

    assertFalse(pi.stderr.replaceFirst("sessionName", "").contains("sessionName"));
    assertTrue(pi.stderr.contains("test_data/session.yaml"));

    pi = processInput(ip, "session list -f test_data/missing.yml", proc);
    assertTrue(pi.stderr.contains("does not exist or is not readable"));
    pi = processInput(ip, "session list -f test_data/broken.yml", proc);
    assertTrue(pi.stderr.contains("Failed to process"));
  }

  @Test
  public void canSaveCurrentSessionInPlace()
      throws IOException,
          InvalidConfigException,
          SQLException,
          InvalidGenomicCoordsException,
          ClassNotFoundException,
          InvalidRecordException,
          InvalidCommandLineException,
          SessionException {
    new Config(null);
    TrackProcessor proc = gimmeTrackProcessor("chr7:1001-1800", 80);
    InteractiveInput ip = new InteractiveInput(new ConsoleReader(), 1);
    new File("tmp.yml").delete();
    new File("tmp2.yml").delete();

    ip.processInput("session save -f tmp.yml foo1", proc);
    String yml = new String(Files.readAllBytes(Paths.get("tmp.yml")));
    assertTrue(yml.contains("foo1"));

    ip.processInput("session save -f tmp.yml foo2", proc);
    yml = new String(Files.readAllBytes(Paths.get("tmp.yml")));
    assertTrue(yml.contains("foo1") && yml.contains("foo2"));

    ip.processInput("session save foo3", proc);
    yml = new String(Files.readAllBytes(Paths.get("tmp.yml")));
    assertTrue(yml.contains("foo1") && yml.contains("foo2") && yml.contains("foo3"));

    ip.processInput("session save -f tmp2.yml foo4", proc);
    yml = new String(Files.readAllBytes(Paths.get("tmp2.yml")));
    assertTrue(yml.contains("foo4") && !yml.contains("foo3"));

    proc.getGenomicCoordsHistory().current().setTo(100000);
    ip.processInput("session save", proc);

    assertEquals(
        100000, (long) new SessionHandler(new File("tmp2.yml")).get("foo4").getGenome().to);

    new File("tmp.yml").delete();
    new File("tmp2.yml").delete();
  }

  @Test
  public void canReplaceGenome()
      throws SQLException,
          InvalidGenomicCoordsException,
          IOException,
          ClassNotFoundException,
          InvalidRecordException,
          InvalidCommandLineException,
          InvalidConfigException {
    new Config(null);
    TrackProcessor proc = gimmeTrackProcessor("chr7:1001-1800", 80);
    InteractiveInput ip = new InteractiveInput(new ConsoleReader(), 1);
    ip.processInput("session open -f test_data/session.yaml newSession", proc);
    String fasta = new String(proc.getGenomicCoordsHistory().current().getSequenceFromFasta());
    assertEquals("TTATT", fasta.substring(0, 5));
    ip.processInput("session open -f test_data/session.yaml fastafile-not-found", proc);
    assertEquals(null, proc.getGenomicCoordsHistory().current().getFastaFile());
    ip.processInput("setGenome test_data/chr7.fa", proc);
    fasta = new String(proc.getGenomicCoordsHistory().current().getSequenceFromFasta());
    assertEquals("ACACG", fasta.substring(0, 5));
  }

  @Test
  public void canOpenSessionWithFileMissing()
      throws IOException,
          InvalidConfigException,
          SQLException,
          InvalidGenomicCoordsException,
          ClassNotFoundException,
          InvalidRecordException,
          InvalidCommandLineException {
    new Config(null);
    TrackProcessor proc = gimmeTrackProcessor("chr7:1001-1800", 80);

    InteractiveInput ip = new InteractiveInput(new ConsoleReader(), 1);

    ByteArrayOutputStream err = new ByteArrayOutputStream();
    System.setErr(new PrintStream(err));
    ip.processInput("session open -f test_data/session.yaml file-not-found", proc);
    String errStr = err.toString();
    assertTrue(errStr.contains("Sequence dictionary"));
    assertTrue(errStr.contains("xs#2"));
    System.setErr(new PrintStream(new FileOutputStream(FileDescriptor.err)));

    assertEquals(5567000, (long) proc.getGenomicCoordsHistory().current().getFrom());
    assertEquals("xs#1", proc.getTrackSet().getTrackList().get(0).getTrackTag());
  }

  @Test
  public void canOpenSessionWithFastaFileMissing()
      throws IOException,
          InvalidConfigException,
          SQLException,
          InvalidGenomicCoordsException,
          ClassNotFoundException,
          InvalidRecordException,
          InvalidCommandLineException {
    new Config(null);
    TrackProcessor proc = gimmeTrackProcessor("chr7:1001-1800", 80);

    InteractiveInput ip = new InteractiveInput(new ConsoleReader(), 1);

    ByteArrayOutputStream err = new ByteArrayOutputStream();
    System.setErr(new PrintStream(err));
    ip.processInput("session open -f test_data/session.yaml fastafile-not-found", proc);
    String errStr = err.toString();
    assertTrue(errStr.contains("missing.fa"));
    System.setErr(new PrintStream(new FileOutputStream(FileDescriptor.err)));
    assertNull(proc.getGenomicCoordsHistory().current().getFastaFile());
    assertEquals("xy#1", proc.getTrackSet().getTrackList().get(0).getTrackTag());
  }

  @Test
  public void canHandleSessionFileNotFound()
      throws IOException,
          InvalidConfigException,
          SQLException,
          InvalidGenomicCoordsException,
          ClassNotFoundException,
          InvalidRecordException,
          InvalidCommandLineException {
    new Config(null);
    TrackProcessor proc = gimmeTrackProcessor("chr7:1001-1800", 80);
    InteractiveInput ip = new InteractiveInput(new ConsoleReader(), 0);
    ProcessInput pi = this.processInput(ip, "session open -f missing.yml foo", proc);
    assertTrue(pi.stderr.contains("does not exist"));

    pi = this.processInput(ip, "session save -f foo/bar/tmp.yaml foo", proc);
    assertTrue(pi.stderr.contains("Directory '"));

    pi = this.processInput(ip, "session save -f /tmp.yaml foo", proc);
    assertTrue(pi.stderr.contains("Cannot write "));
  }

  @Test
  public void canHandleSessionNameNotFound()
      throws IOException,
          InvalidConfigException,
          SQLException,
          InvalidGenomicCoordsException,
          ClassNotFoundException,
          InvalidRecordException,
          InvalidCommandLineException {
    new Config(null);
    TrackProcessor proc = gimmeTrackProcessor("chr7:1001-1800", 80);
    InteractiveInput ip = new InteractiveInput(new ConsoleReader(), 1);
    ProcessInput pi = this.processInput(ip, "session open -f test_data/session.yaml spam", proc);
    assertTrue(pi.stderr.contains("Cannot find session with name 'spam'"));
  }

  @Test
  public void canOpenSession()
      throws IOException,
          InvalidConfigException,
          SQLException,
          InvalidGenomicCoordsException,
          ClassNotFoundException,
          InvalidRecordException,
          InvalidCommandLineException {
    new Config(null);
    TrackProcessor proc = gimmeTrackProcessor("chr7:1001-1800", 80);
    InteractiveInput ip = new InteractiveInput(new ConsoleReader(), 1);
    ip.processInput("session", proc);
    assertEquals(ExitCode.ERROR, ip.getInteractiveInputExitCode());

    ip.processInput("session open foo bar", proc);
    assertEquals(ExitCode.ERROR, ip.getInteractiveInputExitCode());

    ip.processInput("session open -f test_data/session.yaml foobar", proc);
    assertEquals(ExitCode.ERROR, ip.getInteractiveInputExitCode());

    ip.processInput("session open -f test_data/session.yaml 1", proc);
    assertEquals(2, proc.getTrackSet().getTrackList().size());
  }

  @Test
  public void canOpenSessionAndShowGenome()
      throws IOException,
          InvalidConfigException,
          SQLException,
          InvalidGenomicCoordsException,
          ClassNotFoundException,
          InvalidRecordException,
          InvalidCommandLineException {
    new Config(null);
    TrackProcessor proc = gimmeTrackProcessor("chr7:1001-1800", 80);
    InteractiveInput ip = new InteractiveInput(new ConsoleReader(), 1);
    this.processInput(ip, "session open -f test_data/session.yaml 1", proc);
    ProcessInput pi = this.processInput(ip, "show genome", proc);
    assertTrue(pi.stderr.contains("159138663"));
  }

  @Test
  public void canSaveSession()
      throws IOException,
          InvalidConfigException,
          SQLException,
          InvalidGenomicCoordsException,
          ClassNotFoundException,
          InvalidRecordException,
          InvalidCommandLineException,
          SessionException {
    new Config(null);
    TrackProcessor proc = gimmeTrackProcessor("chr7:1001-1800", 80);
    GenomicCoords gc = proc.getGenomicCoordsHistory().current();
    Track tr = new TrackPileup("test_data/ds051.actb.bam", gc);
    proc.getTrackSet().addTrack(tr, "tr#1");
    InteractiveInput ip = new InteractiveInput(new ConsoleReader(), 1);
    ip.processInput("session save -f tmp.yml tr1", proc);
    SessionHandler sh = new SessionHandler(new File("tmp.yml"));
  }

  @Test
  public void canPrintGenome()
      throws IOException,
          InvalidConfigException,
          ClassNotFoundException,
          InvalidGenomicCoordsException,
          InvalidRecordException,
          SQLException,
          InvalidCommandLineException {
    new Config(null);
    TrackProcessor proc;
    InteractiveInput ip = new InteractiveInput(new ConsoleReader(), 1);
    proc = gimmeTrackProcessor("chr7:1001-1800", 80);
    ip.processInput("show ge -n -1", proc);
    ip.processInput("show ge -n 10", proc);
  }

  @Test
  public void canMovePositionByColumn()
      throws InvalidGenomicCoordsException,
          IOException,
          InvalidRecordException,
          ClassNotFoundException,
          SQLException,
          InvalidCommandLineException,
          InvalidConfigException {

    new Config(null);
    TrackProcessor proc;

    InteractiveInput ip = new InteractiveInput(new ConsoleReader(), 1);
    proc = gimmeTrackProcessor("chr7:1001-1800", 80);
    GenomicCoords gc2 = ip.processInput("[", proc).getGenomicCoordsHistory().current();
    assertEquals(1011, (int) gc2.getFrom());
    assertEquals(1810, (int) gc2.getTo());

    proc = gimmeTrackProcessor("chr7:1001-1800", 80);
    gc2 = ip.processInput("]", proc).getGenomicCoordsHistory().current();
    assertEquals(1001 - 10, (int) gc2.getFrom());
    assertEquals(1800 - 10, (int) gc2.getTo());

    proc = gimmeTrackProcessor("chr7:1001-1800", 80);
    gc2 = ip.processInput("[ 20", proc).getGenomicCoordsHistory().current();
    assertEquals(1001 + (20 * 10), (int) gc2.getFrom());
    assertEquals(1800 + (20 * 10), (int) gc2.getTo());

    proc = gimmeTrackProcessor("chr7:1001-1800", 80);
    gc2 = ip.processInput("[20", proc).getGenomicCoordsHistory().current();
    assertEquals(1001 + (20 * 10), (int) gc2.getFrom());
    assertEquals(1800 + (20 * 10), (int) gc2.getTo());

    proc = gimmeTrackProcessor("chr7:1001-1800", 80);
    gc2 = ip.processInput("]] 3", proc).getGenomicCoordsHistory().current();
    assertEquals(1001 - (6 * 10), (int) gc2.getFrom());
    assertEquals(1800 - (6 * 10), (int) gc2.getTo());

    proc = gimmeTrackProcessor("chr7:1001-1800", 80);
    gc2 = ip.processInput("] 0", proc).getGenomicCoordsHistory().current();
    assertEquals(1001, (int) gc2.getFrom());
    assertEquals(1800, (int) gc2.getTo());

    // Test left bound
    proc = gimmeTrackProcessor("chr7:1001-1800", 80);
    gc2 = ip.processInput("]] 30000", proc).getGenomicCoordsHistory().current();
    assertEquals(1, (int) gc2.getFrom());
    assertEquals(800, (int) gc2.getTo());

    // Test right bound
    proc = gimmeTrackProcessor("chr7:1001-1800", 80);
    gc2 = ip.processInput("[ 30000000", proc).getGenomicCoordsHistory().current();
    assertEquals(159138663, (int) gc2.getTo());
    assertEquals(159138663 - 800 + 1, (int) gc2.getFrom());

    // Invalid input
    ip.processInput("[]", proc).getGenomicCoordsHistory().current();
    assertEquals(ExitCode.ERROR, ip.getInteractiveInputExitCode());

    // Ensure this is fine
    ip.processInput("]", proc).getGenomicCoordsHistory().current();
    assertEquals(ExitCode.CLEAN, ip.getInteractiveInputExitCode());

    // Another invalid input
    ip.processInput("] foo", proc).getGenomicCoordsHistory().current();
    assertEquals(ExitCode.ERROR, ip.getInteractiveInputExitCode());
  }

  @Test
  public void canPrintHelp()
      throws InvalidGenomicCoordsException,
          IOException,
          ClassNotFoundException,
          InvalidRecordException,
          SQLException,
          InvalidCommandLineException {

    InteractiveInput ip = new InteractiveInput(null, 1);

    GenomicCoords gc = new GenomicCoords("chr7:1-100", 80, null, null);
    GenomicCoordsHistory gch = new GenomicCoordsHistory();
    gch.add(gc);
    TrackProcessor proc = new TrackProcessor(null, gch);

    // Send the output of the help from stderr to baos
    ByteArrayOutputStream baos = new ByteArrayOutputStream();
    System.setErr(new PrintStream(baos));

    // Various ways of getting general help
    ip.processInput("-h", proc);
    String H1 = baos.toString();
    baos.reset();
    assertTrue(H1.contains("show this help"));

    System.out.println(H1);

    ip.processInput("h", proc);
    String H2 = baos.toString();
    baos.reset();

    ip.processInput("help", proc);
    String H3 = baos.toString();
    baos.reset();

    ip.processInput("  ?", proc);
    String H4 = baos.toString();
    baos.reset();

    assertEquals(H1, H2);
    assertEquals(H1, H3);
    assertEquals(H1, H4);

    // Various way of getting command help
    ip.processInput("next -h", proc);
    String h1 = baos.toString();
    baos.reset();
    assertTrue(h1.contains("Move to the next"));

    ip.processInput("?next", proc);
    String h2 = baos.toString();
    baos.reset();

    ip.processInput("  help  next ", proc);
    String h3 = baos.toString();
    baos.reset();

    assertEquals(h1, h2);
    assertEquals(h1, h3);
  }

  @Test
  public void canGoToRegion() throws InvalidGenomicCoordsException, IOException {

    InteractiveInput ip = new InteractiveInput(null, 1);

    GenomicCoords gc = new GenomicCoords("chr1:1-1000", 80, null, null);
    String region = ip.gotoOnCurrentChrom(Splitter.on(" ").splitToList("1 100"), gc);
    assertEquals("chr1:1-100", region);

    region = ip.gotoOnCurrentChrom(Splitter.on(" ").splitToList("1 20 30 100"), gc);
    assertEquals("chr1:1-100", region);

    region = ip.gotoOnCurrentChrom(Splitter.on(" ").splitToList("2"), gc);
    assertEquals("chr1:2-81", region);
  }

  @Test
  public void canGoToRegionAndCenter() throws InvalidGenomicCoordsException, IOException {

    InteractiveInput ip = new InteractiveInput(null, 0);

    GenomicCoords gc = new GenomicCoords("chr1:1-20", 20, null, null);

    String region = ip.gotoOnCurrentChrom(Splitter.on(" ").splitToList("16 c"), gc);
    assertEquals("chr1:6-25", region);

    gc = new GenomicCoords("chr1:1-21", 20, null, null);

    region = ip.gotoOnCurrentChrom(Splitter.on(" ").splitToList("16 c"), gc);
    assertEquals("chr1:6-26", region);
  }
}
