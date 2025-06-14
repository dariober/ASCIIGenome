package samTextViewer;

import static org.junit.Assert.*;

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
import java.util.regex.Pattern;
import jline.console.ConsoleReader;
import org.junit.Test;
import session.SessionHandler;
import tracks.Track;
import tracks.TrackPileup;
import tracks.TrackSet;

public class InteractiveInputTest {

  public static class ProcessInput {
    public String stderr;
    public String stdout;
  }

  public ProcessInput processInput(InteractiveInput ip, String cmd, TrackProcessor p)
      throws InvalidGenomicCoordsException, IOException {
    ByteArrayOutputStream err = new ByteArrayOutputStream();
    System.setErr(new PrintStream(err));
    ByteArrayOutputStream out = new ByteArrayOutputStream();
    System.setOut(new PrintStream(out));
    ip.processInput(cmd, p);
    String errStr = err.toString();
    System.setErr(new PrintStream(new FileOutputStream(FileDescriptor.err)));
    String outStr = out.toString();
    System.setOut(new PrintStream(new FileOutputStream(FileDescriptor.out)));
    ProcessInput pi = new ProcessInput();
    pi.stderr = errStr;
    pi.stdout = outStr;
    return pi;
  }

  public static TrackProcessor gimmeTrackProcessor(String region, int terminalWidth, String genome)
      throws InvalidGenomicCoordsException,
          IOException,
          ClassNotFoundException,
          InvalidRecordException,
          SQLException {
    GenomicCoords gc = new GenomicCoords(region, terminalWidth, null, null);
    List<String> gf = new ArrayList<String>();
    gf.add(genome);
    gc.setGenome(gf, false);
    GenomicCoordsHistory gch = new GenomicCoordsHistory();
    gch.add(gc);
    TrackSet trackSet = new TrackSet(new ArrayList<String>(), gc);
    TrackProcessor proc = new TrackProcessor(trackSet, gch);
    return proc;
  }

//  public static TrackProcessor gimmeTrackProcessor(String region, int terminalWidth)
//      throws SQLException,
//          InvalidGenomicCoordsException,
//          IOException,
//          ClassNotFoundException,
//          InvalidRecordException {
//    return gimmeTrackProcessor(region, terminalWidth, "test_data/ds051.actb.bam");
//  }

  @Test
  public void canListSessions()
      throws IOException,
          SQLException,
          InvalidGenomicCoordsException,
          ClassNotFoundException,
          InvalidRecordException,
          InvalidCommandLineException {
    TrackProcessor proc = gimmeTrackProcessor("chr7:1001-1800", 80, "test_data/ds051.actb.bam");
    InteractiveInput ip = new InteractiveInput(new ConsoleReader(), 1);
    ProcessInput pi = processInput(ip, "sessionList -f test_data/session.yaml -n 1", proc);
    assertTrue(pi.stderr.contains("sessionName: no-fastafile"));

    assertFalse(pi.stderr.replaceFirst("sessionName", "").contains("sessionName"));
    assertTrue(pi.stderr.contains("test_data/session.yaml"));

    pi = processInput(ip, "sessionList -f test_data/missing.yml", proc);
    assertTrue(pi.stderr.contains("does not exist or is not readable"));
    pi = processInput(ip, "sessionList -f test_data/broken.yml", proc);
    assertTrue(pi.stderr.contains("Failed to process"));
  }

  @Test
  public void canListCurrentSessionNameAndFile()
      throws IOException,
          SQLException,
          InvalidGenomicCoordsException,
          ClassNotFoundException,
          InvalidRecordException,
          InvalidCommandLineException,
          InvalidConfigException {
    new Config(null);
    TrackProcessor proc = gimmeTrackProcessor("chr7:1001-1800", 80, "test_data/ds051.actb.bam");
    InteractiveInput ip = new InteractiveInput(new ConsoleReader(), 1);

    ProcessInput pi = processInput(ip, "sessionList", proc);
    assertTrue(
        Pattern.compile("Session file: /.*/.asciigenome/session.yml").matcher(pi.stderr).find());
    assertTrue(pi.stderr.contains("Current session: n/a"));

    ip.processInput("sessionOpen -f test_data/session.yaml no-fastafile", proc);
    pi = processInput(ip, "sessionList -f test_data/session.yaml no-fastafile", proc);
    assertTrue(
        Pattern.compile("Session file: /.*/test_data/session.yaml").matcher(pi.stderr).find());
    assertTrue(pi.stderr.contains("Current session: no-fastafile"));
  }

  @Test
  public void canDeleteSession()
      throws IOException,
          InvalidConfigException,
          SQLException,
          InvalidGenomicCoordsException,
          ClassNotFoundException,
          InvalidRecordException,
          InvalidCommandLineException {
    new Config(null);
    TrackProcessor proc = gimmeTrackProcessor("chr7:1001-1800", 80, "test_data/ds051.actb.bam");
    InteractiveInput ip = new InteractiveInput(new ConsoleReader(), 1);

    ProcessInput pi = processInput(ip, "sessionDelete", proc);
    pi = processInput(ip, "sessionDelete -f test_data/session.yaml", proc);
    assertEquals(ExitCode.ERROR, ip.getInteractiveInputExitCode());
    assertEquals("Please provide the name of the session to delete", pi.stderr.trim());

    pi = processInput(ip, "sessionDelete", proc);
    assertEquals(ExitCode.ERROR, ip.getInteractiveInputExitCode());
    assertEquals("Please provide the name of the session to delete", pi.stderr.trim());

    File f = new File("tmp.yml");
    f.delete();
    Files.copy(Paths.get("test_data/session.yaml"), Paths.get("tmp.yml"));

    processInput(ip, "sessionDelete -f tmp.yml newSession", proc);
    assertEquals(ExitCode.CLEAN_NO_FLUSH, ip.getInteractiveInputExitCode());
    String txt = new String(Files.readAllBytes(Paths.get("tmp.yml")));
    assertTrue(!txt.contains("newSession"));

    processInput(ip, "sessionDelete -f tmp.yml newSession", proc);
    assertEquals(ExitCode.ERROR, ip.getInteractiveInputExitCode());
    f.delete();
  }

  @Test
  public void canSaveCurrentSessionInPlace()
      throws IOException,
          InvalidConfigException,
          SQLException,
          InvalidGenomicCoordsException,
          ClassNotFoundException,
          InvalidRecordException,
          SessionException {
    new Config(null);
    TrackProcessor proc = gimmeTrackProcessor("chr7:1001-1800", 80, "test_data/ds051.actb.bam");
    InteractiveInput ip = new InteractiveInput(new ConsoleReader(), 1);
    new File("tmp.yml").delete();
    new File("tmp2.yml").delete();

    ip.processInput("sessionSave -f tmp.yml foo1", proc);
    String yml = new String(Files.readAllBytes(Paths.get("tmp.yml")));
    assertTrue(yml.contains("foo1"));

    ip.processInput("sessionSave -f tmp.yml foo2", proc);
    yml = new String(Files.readAllBytes(Paths.get("tmp.yml")));
    assertTrue(yml.contains("foo1") && yml.contains("foo2"));

    ip.processInput("sessionSave foo3", proc);
    yml = new String(Files.readAllBytes(Paths.get("tmp.yml")));
    assertTrue(yml.contains("foo1") && yml.contains("foo2") && yml.contains("foo3"));

    ip.processInput("sessionSave -f tmp2.yml foo4", proc);
    yml = new String(Files.readAllBytes(Paths.get("tmp2.yml")));
    assertTrue(yml.contains("foo4") && !yml.contains("foo3"));

    proc.getGenomicCoordsHistory().current().setTo(100000);
    ip.processInput("sessionSave", proc);

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
          InvalidConfigException {
    new Config(null);
    TrackProcessor proc = gimmeTrackProcessor("chr7:1001-1800", 80, "test_data/ds051.actb.bam");
    InteractiveInput ip = new InteractiveInput(new ConsoleReader(), 1);
    ip.processInput("sessionOpen -f test_data/session.yaml newSession", proc);
    String fasta = new String(proc.getGenomicCoordsHistory().current().getSequenceFromFasta());
    assertEquals("TTATT", fasta.substring(0, 5));
    ip.processInput("sessionOpen -f test_data/session.yaml fastafile-not-found", proc);
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
          InvalidRecordException {
    new Config(null);
    TrackProcessor proc = gimmeTrackProcessor("chr7:1001-1800", 80, "test_data/ds051.actb.bam");

    InteractiveInput ip = new InteractiveInput(new ConsoleReader(), 1);

    ByteArrayOutputStream err = new ByteArrayOutputStream();
    System.setErr(new PrintStream(err));
    ip.processInput("sessionOpen -f test_data/session.yaml file-not-found", proc);
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
          InvalidRecordException {
    new Config(null);
    TrackProcessor proc = gimmeTrackProcessor("chr7:1001-1800", 80, "test_data/ds051.actb.bam");

    InteractiveInput ip = new InteractiveInput(new ConsoleReader(), 1);

    ByteArrayOutputStream err = new ByteArrayOutputStream();
    System.setErr(new PrintStream(err));
    ip.processInput("sessionOpen -f test_data/session.yaml fastafile-not-found", proc);
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
    TrackProcessor proc = gimmeTrackProcessor("chr7:1001-1800", 80, "test_data/ds051.actb.bam");
    InteractiveInput ip = new InteractiveInput(new ConsoleReader(), 0);
    ProcessInput pi = this.processInput(ip, "sessionOpen -f missing.yml foo", proc);
    assertTrue(pi.stderr.contains("does not exist"));

    pi = this.processInput(ip, "sessionSave -f foo/bar/tmp.yaml foo", proc);
    assertTrue(pi.stderr.contains("Directory '"));

    pi = this.processInput(ip, "sessionSave -f /tmp.yaml foo", proc);
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
    TrackProcessor proc = gimmeTrackProcessor("chr7:1001-1800", 80, "test_data/ds051.actb.bam");
    InteractiveInput ip = new InteractiveInput(new ConsoleReader(), 1);
    ProcessInput pi = this.processInput(ip, "sessionOpen -f test_data/session.yaml spam", proc);
    assertTrue(pi.stderr.contains("Cannot find session with name 'spam'"));
  }

  @Test
  public void canOpenSession()
      throws IOException,
          InvalidConfigException,
          SQLException,
          InvalidGenomicCoordsException,
          ClassNotFoundException,
          InvalidRecordException {
    new Config(null);
    TrackProcessor proc = gimmeTrackProcessor("chr7:1001-1800", 80, "test_data/ds051.actb.bam");
    InteractiveInput ip = new InteractiveInput(new ConsoleReader(), 1);
    ip.processInput("session", proc);
    assertEquals(ExitCode.ERROR, ip.getInteractiveInputExitCode());

    ip.processInput("sessionOpen foo bar", proc);
    assertEquals(ExitCode.ERROR, ip.getInteractiveInputExitCode());

    ip.processInput("sessionOpen -f test_data/session.yaml foobar", proc);
    assertEquals(ExitCode.ERROR, ip.getInteractiveInputExitCode());

    ip.processInput("sessionOpen -f test_data/session.yaml 1", proc);
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
    TrackProcessor proc = gimmeTrackProcessor("chr7:1001-1800", 80, "test_data/ds051.actb.bam");
    InteractiveInput ip = new InteractiveInput(new ConsoleReader(), 1);
    Files.deleteIfExists(new File("test_data/tmp.yml").toPath());
    this.processInput(ip, "sessionSave -f test_data/tmp.yml test", proc);
    this.processInput(ip, "sessionOpen -f test_data/tmp.yml test", proc);
    ProcessInput pi = this.processInput(ip, "show genome", proc);
    assertTrue(pi.stderr.contains("159138663"));
    Files.deleteIfExists(new File("test_data/tmp.yml").toPath());
  }

  @Test
  public void canSaveSession()
      throws IOException,
          InvalidConfigException,
          SQLException,
          InvalidGenomicCoordsException,
          ClassNotFoundException,
          InvalidRecordException,
          SessionException {
    new Config(null);
    TrackProcessor proc = gimmeTrackProcessor("chr7:1001-1800", 80, "test_data/ds051.actb.bam");
    GenomicCoords gc = proc.getGenomicCoordsHistory().current();
    Track tr = new TrackPileup("test_data/ds051.actb.bam", gc);
    proc.getTrackSet().addTrack(tr, "tr#1");
    InteractiveInput ip = new InteractiveInput(new ConsoleReader(), 1);
    ip.processInput("sessionSave -f tmp.yml tr1", proc);
    SessionHandler sh = new SessionHandler(new File("tmp.yml"));
    Files.deleteIfExists(new File("tmp.yml").toPath());
  }

  @Test
  public void canPrintGenome()
      throws IOException,
          InvalidConfigException,
          ClassNotFoundException,
          InvalidGenomicCoordsException,
          InvalidRecordException,
          SQLException {
    new Config(null);
    TrackProcessor proc;
    InteractiveInput ip = new InteractiveInput(new ConsoleReader(), 1);
    proc = gimmeTrackProcessor("chr7:1001-1800", 80, "test_data/ds051.actb.bam");

    ProcessInput pi = this.processInput(ip, "show ge -n -1", proc);
    assertTrue(!pi.stderr.contains("omitted"));
    String[] out = pi.stderr.split("\n");
    assertEquals("Genome size: 3095693983; Number of contigs: 25", out[0].trim());
    assertEquals(26, out.length);

    pi = this.processInput(ip, "show ge -n 25", proc);
    assertTrue(!pi.stderr.contains("omitted"));
    out = pi.stderr.split("\n");
    assertEquals(26, out.length);

    pi = this.processInput(ip, "show ge -n 24", proc);
    out = pi.stderr.split("\n");
    assertEquals(26, out.length);
    assertTrue(pi.stderr.trim().endsWith("1 contigs omitted]"));

    pi = this.processInput(ip, "show ge -n 23", proc);
    out = pi.stderr.split("\n");
    assertEquals(25, out.length);
    assertTrue(pi.stderr.trim().endsWith("2 contigs omitted]"));
  }

  @Test
  public void canPrintTranslation()
      throws SQLException,
          InvalidGenomicCoordsException,
          IOException,
          ClassNotFoundException,
          InvalidRecordException,
          InvalidConfigException {
    new Config(null);
    TrackProcessor proc;
    proc = gimmeTrackProcessor("chr7:10001-10061", 80, "test_data/chr7.fa");
    proc.setNoFormat(true);
    InteractiveInput ip = new InteractiveInput(new ConsoleReader(), 1);
    ProcessInput pi = this.processInput(ip, "translate", proc);
    assertTrue(pi.stdout.contains("\nctaaccctaaccctaaccctaaccctaaccctaaccctaaccctaaccctaaccctaaccc\n"));
    assertTrue(pi.stdout.contains("\n3 *  P  *  P  *  P  *  P  *  P  *  P  *  P  *  P  *  P  *  P \n"));
    assertTrue(pi.stdout.contains("\n2L  T  L  T  L  T  L  T  L  T  L  T  L  T  L  T  L  T  L  T  \n"));
    assertTrue(pi.stdout.contains("\n1  N  P  N  P  N  P  N  P  N  P  N  P  N  P  N  P  N  P  N   \n"));
    assertTrue(pi.stdout.contains("\n   V  R  V  R  V  R  V  R  V  R  V  R  V  R  V  R  V  R  V  1\n"));
    assertTrue(pi.stdout.contains("\n  L  G  L  G  L  G  L  G  L  G  L  G  L  G  L  G  L  G  L  G2\n"));
    assertTrue(pi.stdout.contains("\n *  G  *  G  *  G  *  G  *  G  *  G  *  G  *  G  *  G  *  G 3\n"));

    pi = this.processInput(ip, "translate -geneticCode vertebrate_mitochondrial", proc);
    assertTrue(pi.stdout.contains("\n   V  *  V  *  V  *  V  *  V  *  V  *  V  *  V  *  V  *  V  1\n"));

    pi = this.processInput(ip, "+10000", proc);
    assertTrue(pi.stdout.contains("\n1 A  *  K  A  D  T  S  S  K  S  I  T  E  A  M  V  Q  P  K  L \n"));

    pi = this.processInput(ip, "translate -codon start_and_stop", proc);
    assertTrue(pi.stdout.contains("    *                                      M                \n"));

    pi = this.processInput(ip, "translate -codon start", proc);
    assertTrue(pi.stdout.contains("                                           M                \n"));

    pi = this.processInput(ip, "translate -codon stop", proc);
    assertTrue(pi.stdout.contains("    *                                                       \n"));

    pi = this.processInput(ip, "translate -codon all", proc);
    assertTrue(pi.stdout.contains("C  Q  E  S  *  H  I  I  K  I  H  Y  *  G  Y  S  S  A  K  A  \n"));

    pi = this.processInput(ip, "zo", proc);
    assertFalse(pi.stdout.contains("C  Q  E"));

    // MEMO: zo && zi does not return to *exactly* the same position
    pi = this.processInput(ip, "zi", proc);
    assertTrue(pi.stdout.contains("C  Q  E"));

    pi = this.processInput(ip, "translate -frame none", proc);
    assertFalse(pi.stdout.contains("C  Q  E"));
  }

  @Test
  public void canHandleTranslateWithoutReference() throws IOException, InvalidConfigException, SQLException, InvalidGenomicCoordsException, ClassNotFoundException, InvalidRecordException {
    new Config(null);
    TrackProcessor proc;
    proc = gimmeTrackProcessor("chr7:10001-10061", 80, "test_data/ds051.actb.bam");
    proc.setNoFormat(true);
    InteractiveInput ip = new InteractiveInput(new ConsoleReader(), 1);
    ProcessInput pi = this.processInput(ip, "translate", proc);
    assertEquals("Unable to translate without a reference sequence", pi.stderr.trim());
  }

  @Test
  public void canMoveToNextChrom()
      throws SQLException,
          InvalidGenomicCoordsException,
          IOException,
          ClassNotFoundException,
          InvalidRecordException {
    TrackProcessor proc = gimmeTrackProcessor("small", 150, "test_data/seq_cg.fa");
    proc.setNoFormat(true);
    InteractiveInput ip = new InteractiveInput(new ConsoleReader(), 1);
    ProcessInput pi = this.processInput(ip, "nextChrom", proc);
    assertTrue(pi.stdout.contains("seq:1-120"));
    assertTrue(pi.stdout.contains("ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA"));
  }

  @Test
  public void canOmitSequenceWhenZoomout()
      throws SQLException, InvalidGenomicCoordsException, IOException, ClassNotFoundException, InvalidRecordException {
    TrackProcessor proc = gimmeTrackProcessor("chr7:10000-10060", 80, "test_data/chr7.fa");
    proc.setNoFormat(true);
    InteractiveInput ip = new InteractiveInput(new ConsoleReader(), 1);
    ProcessInput pi = this.processInput(ip, "goto chr7:200000-200060", proc);
    assertTrue(pi.stdout.contains("\nTTCTTGACACTGATTGATCTGCCAAAAGGGGAAGAATGAGTCCAGCTAGAATCCAGGACTA\n"));
    pi = this.processInput(ip, "zo", proc);
    assertFalse(pi.stdout.contains("TTCTTGACACTGATTGATCTGCCAAAAGGGGAAGAATGAGTCCAGCTAGAATCCAGGACTA"));
    pi = this.processInput(ip, "zi", proc);
    assertTrue(pi.stdout.contains("TTCTTGACACTGATTGATCTGCCAAAAGGGGAAGAATGAG"));
  }

  @Test
  public void canUpdateSequenceWhenMovingCoords()
      throws SQLException, InvalidGenomicCoordsException, IOException, ClassNotFoundException, InvalidRecordException {
    TrackProcessor proc = gimmeTrackProcessor("chr7:10000-10060", 80, "test_data/chr7.fa");
    proc.setNoFormat(true);
    InteractiveInput ip = new InteractiveInput(new ConsoleReader(), 1);
    ProcessInput pi = this.processInput(ip, "+1", proc);
    assertTrue(pi.stdout.contains("\nctaaccctaaccctaaccctaaccctaaccctaaccctaaccctaaccctaaccctaaccc\n"));

    pi = this.processInput(ip, "+1", proc);
    assertTrue(pi.stdout.contains("\ntaaccctaaccctaaccctaaccctaaccctaaccctaaccctaaccctaaccctaaccct\n"));

    pi = this.processInput(ip, "[", proc);
    assertTrue(pi.stdout.contains("\naaccctaaccctaaccctaaccctaaccctaaccctaaccctaaccctaaccctaacccta\n"));

    pi = this.processInput(ip, "100000", proc);
    assertTrue(pi.stdout.contains("\ncagaaggaaaacgggaaacttcacaattagtgaatatttaaaaacagactcttaagaaacc\n"));

    pi = this.processInput(ip, "goto chr7:100010-100070", proc);
    assertTrue(pi.stdout.contains("\nacgggaaacttcacaattagtgaatatttaaaaacagactcttaagaaaccaaaggatcaa\n"));

    pi = this.processInput(ip, "0.17 0.33", proc);
    assertTrue(pi.stdout.contains("\ntcacaattagt\n"));

    pi = this.processInput(ip, "zo", proc);
    assertTrue(pi.stdout.contains("\naaacttcacaattagtgaata\n"));

    this.processInput(ip, "goto chr7:100010-100070", proc);
    this.processInput(ip, "f", proc);
    pi = this.processInput(ip, "ff", proc);
    assertTrue(pi.stdout.contains("\ngactcttaagaaaccaaaggatcaaggaagataccacagggaaaaatagagaatatctcaa\n"));

    this.processInput(ip, "goto chr7:100010-100070", proc);
    this.processInput(ip, "bb", proc);
    pi = this.processInput(ip, "b", proc);
    assertTrue(pi.stdout.contains("\naaaggaatgaaactagaaatcaacagcagaaggaaaacgggaaacttcacaattagtgaat\n"));

    this.processInput(ip, "goto chr7:100010-100020", proc);
    pi = this.processInput(ip, "extend 10", proc);
    assertTrue(pi.stdout.contains("\ncagaaggaaaacgggaaacttcacaattagt\n"));

    this.processInput(ip, "goto chr7:200000-200060", proc);
    this.processInput(ip, "goto chr7:200010-200070", proc);
    pi = this.processInput(ip, "p", proc);
    assertTrue(pi.stdout.contains("\nTTCTTGACACTGATTGATCTGCCAAAAGGGGAAGAATGAGTCCAGCTAGAATCCAGGACTA\n"));

    this.processInput(ip, "goto chr7:200000-200060", proc);
    this.processInput(ip, "goto chr7:200010-200070", proc);
    this.processInput(ip, "p", proc);
    pi = this.processInput(ip, "n", proc);
    assertTrue(pi.stdout.contains("\nTGATTGATCTGCCAAAAGGGGAAGAATGAGTCCAGCTAGAATCCAGGACTAACCAGCGGGT\n"));

    this.processInput(ip, "goto chr7:200000-200060", proc);
    this.processInput(ip, "open test_data/hg19_genes.gtf.gz", proc);
    pi = this.processInput(ip, "next -c", proc);
    assertTrue(pi.stdout.contains("\nTGACCCTGTTTCTCTCCCTCCTTCCTGCAGCCATGAAGTCGGGGGGCACGCAGCTGAAGCT\n"));
    pi = this.processInput(ip, "next -c", proc);
    assertTrue(pi.stdout.contains("\nACCGTTTCTGTTTCTGTCTTGTTTTCTCAGACAAACGAGGGAGCAGGAGACACCCCCTGAC\n"));

    pi = this.processInput(ip, "zo", proc);
    assertFalse(pi.stdout.contains("ACCGTTTCTGTTTCTGTCTTGTTTTCTCAGACAAACGAGGGAGCAGGAGACACCCCCTGAC"));
  }

  @Test
  public void canMovePositionByColumn()
      throws InvalidGenomicCoordsException,
          IOException,
          InvalidRecordException,
          ClassNotFoundException,
          SQLException,
          InvalidConfigException {

    new Config(null);
    TrackProcessor proc;

    InteractiveInput ip = new InteractiveInput(new ConsoleReader(), 1);
    proc = gimmeTrackProcessor("chr7:1001-1800", 80, "test_data/ds051.actb.bam");
    GenomicCoords gc2 = ip.processInput("[", proc).getGenomicCoordsHistory().current();
    assertEquals(1011, (int) gc2.getFrom());
    assertEquals(1810, (int) gc2.getTo());

    proc = gimmeTrackProcessor("chr7:1001-1800", 80, "test_data/ds051.actb.bam");
    gc2 = ip.processInput("]", proc).getGenomicCoordsHistory().current();
    assertEquals(1001 - 10, (int) gc2.getFrom());
    assertEquals(1800 - 10, (int) gc2.getTo());

    proc = gimmeTrackProcessor("chr7:1001-1800", 80, "test_data/ds051.actb.bam");
    gc2 = ip.processInput("[ 20", proc).getGenomicCoordsHistory().current();
    assertEquals(1001 + (20 * 10), (int) gc2.getFrom());
    assertEquals(1800 + (20 * 10), (int) gc2.getTo());

    proc = gimmeTrackProcessor("chr7:1001-1800", 80, "test_data/ds051.actb.bam");
    gc2 = ip.processInput("[20", proc).getGenomicCoordsHistory().current();
    assertEquals(1001 + (20 * 10), (int) gc2.getFrom());
    assertEquals(1800 + (20 * 10), (int) gc2.getTo());

    proc = gimmeTrackProcessor("chr7:1001-1800", 80, "test_data/ds051.actb.bam");
    gc2 = ip.processInput("]] 3", proc).getGenomicCoordsHistory().current();
    assertEquals(1001 - (6 * 10), (int) gc2.getFrom());
    assertEquals(1800 - (6 * 10), (int) gc2.getTo());

    proc = gimmeTrackProcessor("chr7:1001-1800", 80, "test_data/ds051.actb.bam");
    gc2 = ip.processInput("] 0", proc).getGenomicCoordsHistory().current();
    assertEquals(1001, (int) gc2.getFrom());
    assertEquals(1800, (int) gc2.getTo());

    // Test left bound
    proc = gimmeTrackProcessor("chr7:1001-1800", 80, "test_data/ds051.actb.bam");
    gc2 = ip.processInput("]] 30000", proc).getGenomicCoordsHistory().current();
    assertEquals(1, (int) gc2.getFrom());
    assertEquals(800, (int) gc2.getTo());

    // Test right bound
    proc = gimmeTrackProcessor("chr7:1001-1800", 80, "test_data/ds051.actb.bam");
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
      throws InvalidGenomicCoordsException, IOException, InvalidRecordException {

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
  public void canGoToInteractive()
      throws SQLException, InvalidGenomicCoordsException, IOException, ClassNotFoundException, InvalidRecordException {
    TrackProcessor proc = gimmeTrackProcessor("chr7:10000-10060", 80, "test_data/chr7.fa");
    proc.setNoFormat(true);
    InteractiveInput ip = new InteractiveInput(new ConsoleReader(), 1);
    ProcessInput pi = this.processInput(ip, "goto chr7:200000-200060", proc);
    assertTrue(pi.stdout.contains("chr7:200000-200060"));
    assertTrue(pi.stdout.contains("\n200000 "));

    pi = this.processInput(ip, "goto chr7:100-10", proc);
    assertTrue(pi.stderr.contains("Start coordinate (100) is larger than end coordinate (10)"));

    pi = this.processInput(ip, "goto chr7:0-10", proc);
    assertTrue(pi.stderr.contains("Start position must be greater or equal to 1"));

    pi = this.processInput(ip, "goto chr7:10000", proc);
    assertTrue(pi.stdout.contains("chr7:10000-10080;"));

    pi = this.processInput(ip, "goto chr7:159138663", proc);
    assertTrue(pi.stdout.contains("-159138663;")); // Something like chr7:159138583-159138663
    System.out.println(pi.stdout);
    System.err.println(pi.stderr);
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
