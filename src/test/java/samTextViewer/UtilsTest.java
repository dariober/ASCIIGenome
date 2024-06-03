package samTextViewer;

import static org.hamcrest.CoreMatchers.is;
import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertThat;
import static org.junit.Assert.assertTrue;

import colouring.Config;
import com.google.common.base.Joiner;
import com.google.common.base.Stopwatch;
import exceptions.InvalidColourException;
import exceptions.InvalidCommandLineException;
import exceptions.InvalidConfigException;
import exceptions.InvalidGenomicCoordsException;
import exceptions.InvalidRecordException;
import faidx.UnindexableFastaFileException;
import filter.FirstOfPairFilter;
import filter.FlagToFilter;
import filter.ReadNegativeStrandFilter;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.filter.AlignedFilter;
import htsjdk.samtools.filter.MappingQualityFilter;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintStream;
import java.net.MalformedURLException;
import java.nio.charset.StandardCharsets;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.regex.Pattern;
import jline.console.ConsoleReader;
import jline.console.history.History;
import jline.console.history.MemoryHistory;
import org.apache.commons.io.FileUtils;
import org.junit.BeforeClass;
import org.junit.Test;
import tracks.IntervalFeature;
import tracks.TrackFormat;
import tracks.TrackReads;
import utils.Tokenizer;

public class UtilsTest {

  @BeforeClass
  public static void init() throws IOException, InvalidConfigException {
    new Config(null);
  }

  static SamReaderFactory srf = SamReaderFactory.make();
  static SamReader samReader = srf.open(new File("test_data/ds051.short.bam"));
  public static SAMSequenceDictionary samSeqDict =
      samReader.getFileHeader().getSequenceDictionary();

  public static String fastaFile = "test_data/chr7.fa";

  @Test
  public void test() throws IOException, InterruptedException {
    // Input list is [foo0, foo1, ..., foo99999]
    List<String> inputList = new ArrayList<String>();
    for (int i = 0; i < 10; i++) {
      inputList.add("foo" + i);
    }
    for (int i = 0; i < 10; i++) {
      inputList.add("bar" + i);
    }

    List<String> cmd = new ArrayList<String>();
    cmd.add("grep");
    cmd.add("foo");

    for (int i = 0; i < 100; i++) {
      ArrayList<String> results =
          Utils.execSystemCommand(inputList.toArray(new String[inputList.size()]), cmd);
      assertEquals(10, results.size());
      // System.out.println("First element:" + results.get(0));
      // System.out.println("Last element:" + results.get(results.size() - 1));
    }
  }

  @Test
  public void canGetSamReaderFromBam() throws MalformedURLException {
    SamReader sr = Utils.getSamReader("test_data/ds051.actb.bam", null);
    SAMRecordIterator iter = sr.query("chr7", 1, 5566791, false);
    assertTrue(iter.hasNext());
  }

  @Test
  public void canGetSamReaderFromCram() throws MalformedURLException {
    SamReader sr = Utils.getSamReader("test_data/ds051.actb.cram", "test_data/chr7.fa");
    SAMRecordIterator iter = sr.query("chr7", 1, 5566791, false);
    assertTrue(iter.hasNext());
  }

  @Test
  public void canGuessTrackType() throws Exception {

    // MEMO:
    // * Check URL
    // * Check empty

    List<TrackFormat> fmt = new ArrayList<TrackFormat>();
    for (TrackFormat x : TrackFormat.values()) {
      fmt.add(x);
    }

    assertEquals(TrackFormat.VCF, Utils.guessTrackType("test_data/malformed.vcf.gz", null));
    assertEquals(TrackFormat.VCF, Utils.guessTrackType("test_data/invalid.vcf", null));
    assertEquals(
        TrackFormat.VCF, Utils.guessTrackType("test_data/CEU.exon.2010_06.genotypes.vcf", null));
    assertEquals(
        TrackFormat.VCF, Utils.guessTrackType("test_data/CEU.exon.2010_06.genotypes.vcf.gz", null));

    assertEquals(TrackFormat.BAM, Utils.guessTrackType("test_data/ds051.actb.bam", null));
    assertEquals(TrackFormat.BAM, Utils.guessTrackType("test_data/ds051.noindex.sam", null));
    fmt.remove(TrackFormat.BAM);

    assertEquals(
        TrackFormat.BAM, Utils.guessTrackType("test_data/ds051.actb.cram", "test_data/chr7.fa"));
    fmt.remove(TrackFormat.BAM);

    assertEquals(
        TrackFormat.GFF,
        Utils.guessTrackType("test_data/Homo_sapiens.GRCh38.86.ENST00000331789.gff3", null));
    assertEquals(TrackFormat.GFF, Utils.guessTrackType("test_data/ovl.gff", null));
    fmt.remove(TrackFormat.GFF);

    // assertEquals(TrackFormat.GTF, Utils.guessTrackType("test_data/hg19_genes_head.gtf"));

    // assertEquals(TrackFormat.BED, Utils.guessTrackType("test_data/refSeq.hg19.short.bed"));

    if (!fmt.isEmpty()) {
      // throw new Exception("These file types have not been tested: " + fmt);
    }
    ;
  }

  @Test
  public void canInterpretStringAsRegexOrAwkAndFilter() throws IOException {
    String[] rawrecs = {"foo\t20", "bar\t30", "baz\t40"};

    boolean[] matched = Utils.matchByAwkOrRegex(rawrecs, "ba");
    assertArrayEquals(new boolean[] {false, true, true}, matched);

    matched = Utils.matchByAwkOrRegex(rawrecs, "'$2 > 25'");
    assertArrayEquals(new boolean[] {false, true, true}, matched);

    matched = Utils.matchByAwkOrRegex(rawrecs, "'$2 > 25 && $1 == \"bar\"'");
    assertArrayEquals(new boolean[] {false, true, false}, matched);

    matched = Utils.matchByAwkOrRegex(rawrecs, "foo|40$");
    assertArrayEquals(new boolean[] {true, false, true}, matched);
  }

  @Test
  public void canTestForOverlappingSegments() {
    // Intersect/contained
    assertTrue(Utils.isOverlapping(1, 10, 8, 20));
    assertTrue(Utils.isOverlapping(1, 10, 2, 5));
    assertTrue(Utils.isOverlapping(1, 10, -2, 5));

    // One) bp overlap
    assertTrue(Utils.isOverlapping(1, 10, 10, 15));
    assertTrue(Utils.isOverlapping(1, 10, -10, 1));

    assertTrue(!Utils.isOverlapping(1, 10, 11, 20));
    assertTrue(!Utils.isOverlapping(1, 10, -10, 0));

    boolean pass = false;
    try {
      Utils.isOverlapping(10, 1, 11, 20);
    } catch (ArithmeticException e) {
      pass = true;
    }
    assertTrue(pass);
  }

  @Test
  public void canGetBoolean() {

    assertTrue(Utils.asBoolean("true"));
    assertTrue(Utils.asBoolean("T"));
    assertTrue(Utils.asBoolean("Y"));
    assertTrue(Utils.asBoolean("ye"));
    assertTrue(Utils.asBoolean("On"));

    assertTrue(!Utils.asBoolean("FALSE"));
    assertTrue(!Utils.asBoolean("N"));
    assertTrue(!Utils.asBoolean("OFF"));

    boolean pass = false;
    try {
      Utils.asBoolean("FOO");
    } catch (IllegalArgumentException e) {
      pass = true;
    }
    assertTrue(pass);

    pass = false;
    try {
      Utils.asBoolean("");
    } catch (IllegalArgumentException e) {
      pass = true;
    }
    assertTrue(pass);
  }

  @Test
  public void canWinsoriseData() {

    int[] ints = new int[] {-50, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 100};
    List<Float> x = new ArrayList<Float>();
    for (int z : ints) {
      x.add((float) z);
    }
    List<Float> w = Utils.winsor2(x, 3.0);
    assertEquals(-7.8434, w.get(0), 0.0001);
    assertEquals(1, w.get(1), 0.0001);
    assertEquals(18.8434, w.get(w.size() - 1), 0.0001);

    w = Utils.winsor2(x, 22);
    assertEquals(x, w);

    w = Utils.winsor2(x, 1);
    assertEquals(1.0522, w.get(0), 0.0001);
    assertEquals(1.0522, w.get(1), 0.0001);
    assertEquals(2, w.get(2), 0.0001);
    assertEquals(9.9478, w.get(w.size() - 1), 0.0001);

    boolean pass = false;
    try {
      Utils.winsor2(x, 0);
    } catch (Exception e) {
      pass = true;
    }
    assertTrue(pass);
  }

  @Test
  public void canGetVCFHeaderAsString() {
    VCFFileReader reader =
        new VCFFileReader(new File("test_data/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz"));
    VCFHeader hdr = reader.getFileHeader();
    reader.close();
    List<String> str = Utils.vcfHeaderToStrings(hdr);
    assertEquals(71, str.size());
    assertTrue(str.get(str.size() - 1).startsWith("#CHROM"));

    reader = new VCFFileReader(new File("test_data/CEU.exon.2010_06.genotypes.vcf.gz"));
    hdr = reader.getFileHeader();
    reader.close();
    str = Utils.vcfHeaderToStrings(hdr);
    assertEquals(12, str.size());
    System.err.println(str);
  }

  @Test
  public void canRoundNumbersInString() {
    String x =
        "chr9.1234x\t9.1234\t2.2 003.1234 DP=4.467; XX=9.987,8.987; \"6.6\" AC=\"8.8\" AC=9.9876="
            + " x5.5 7.7x 1-1.2345 x\"9.99\";\tBAR";
    String fmt = Utils.roundNumbers(x, 2, TrackFormat.BED);
    assertTrue(fmt.startsWith("chr9.1234x\t"));
    assertTrue(fmt.contains("\t9.12\t"));
    assertTrue(fmt.contains(" 003.1234 ")); // More than one leading zeros make NaN
    assertTrue(fmt.contains("DP=4.47;"));
    assertTrue(fmt.contains("XX=9.99,8.99;"));
    assertTrue(fmt.endsWith("\tBAR"));
    assertTrue(fmt.contains(" AC=9.9876= x5.5 7.7x 1-1.2345 ")); // Not rounded

    x = "chrx\t10\t20";
    fmt = Utils.roundNumbers(x, 2, TrackFormat.BED);
    assertEquals(x, fmt);

    x = "chrx\t10\t20\t9.9911";
    fmt = Utils.roundNumbers(x, 2, TrackFormat.BED);
    assertEquals("chrx\t10\t20\t9.99", fmt);

    x = "chrx\t10\t20\t9.9910\tFOO";
    fmt = Utils.roundNumbers(x, 2, TrackFormat.BED);
    assertEquals("chrx\t10\t20\t9.99\tFOO", fmt);

    x = "chrx\t10\t20\t9.9910";
    fmt = Utils.roundNumbers(x, 0, TrackFormat.BED);
    assertEquals("chrx\t10\t20\t10", fmt);

    x = "chrx\t10\t20\t9.9910";
    fmt = Utils.roundNumbers(x, -1, TrackFormat.BED);
    assertEquals(x, fmt);

    x = "chrx 10 20 \"9.9910\"";
    fmt = Utils.roundNumbers(x, 2, TrackFormat.BED);
    assertEquals(x, fmt);

    x = "chrx 10 20 \"9.9910\""; // Round inside quotes
    fmt = Utils.roundNumbers(x, 2, TrackFormat.GTF);
    assertEquals("chrx 10 20 \"9.99\"", fmt);
  }

  @Test
  public void canTestForEqualReadNames() {
    Stopwatch sw = Stopwatch.createStarted();
    for (int i = 0; i < 10000000; i++) {
      Utils.equalReadNames(
          "HSQ9103:403:C6F0HANXX:5:2302:20709:5219", "HSQ9103:403:C6F0HANXX:5:2302:20709:5219");
    }
    System.err.println(sw);
    assertTrue(Utils.equalReadNames("foo", "foo"));
    assertTrue(Utils.equalReadNames("foo index1", "foo index2"));
    assertTrue(Utils.equalReadNames("foo/1", "foo/2"));
    assertTrue(!Utils.equalReadNames("foo", "bar"));
    assertTrue(!Utils.equalReadNames("foo/1foo", "foo/2foo"));
    assertTrue(!Utils.equalReadNames("/1", "/2"));
  }

  @Test
  public void canSortAndIndexCram() throws IOException {
    Utils.sortAndIndexSamOrBam(
        "test_data/ds051.noindex.cram", "sorted.bam", true, "test_data/chr7.fa");
    assertTrue(new File("sorted.bai").length() > 1000);
    assertTrue(new File("sorted.bam").length() > 1000);

    boolean pass = false;
    try {
      Utils.sortAndIndexSamOrBam(
          "test_data/ds051.noindex.cram", "sorted.cram", true, "test_data/chr7.fa");
    } catch (IOException e) {
      assertTrue(e.getMessage().contains("CRAM output"));
      pass = true;
    }
    assertTrue(pass);
  }

  @Test
  public void canSortAndIndexSamOrBam() throws IOException {

    Utils.sortAndIndexSamOrBam("test_data/ds051.noindex.bam", "sorted.bam", true, null);
    assertTrue(new File("sorted.bai").length() > 1000);
    assertTrue(new File("sorted.bam").length() > 1000);

    // With SAM input
    Utils.sortAndIndexSamOrBam("test_data/ds051.noindex.sam", "sorted1.bam", true, null);
    assertTrue(new File("sorted1.bai").length() > 1000);
    assertTrue(new File("sorted1.bam").length() > 1000);

    // Works also with URL
    File sorted2 = new File("sorted2.bam");

    Utils.sortAndIndexSamOrBam(
        "https://raw.githubusercontent.com/dariober/ASCIIGenome/master/test_data/ds051.noindex.bam",
        sorted2.getAbsolutePath(),
        true,
        null);
    assertTrue(new File("sorted2.bai").length() > 1000);
    assertTrue(new File("sorted2.bam").length() > 1000);
  }

  @Test
  public void roundNumber() {
    assertEquals(10.12, Utils.round(10.123, 2), 0.000001);
    assertEquals(10.13, Utils.round(10.129, 2), 0.000001);
    assertEquals(10.50, Utils.round(10.505, 2), 0.000001);
    assertEquals(10.52, Utils.round(10.515, 2), 0.000001);
  }

  @Test
  public void testParallel() {
    List<List<String>> list = new ArrayList<List<String>>();
    List<String> inList = new ArrayList<String>();
    inList.add("foo");
    inList.add("bar");
    inList.add("baz");
    List<String> inList2 = new ArrayList<String>();
    inList2.add("foo2");
    inList2.add("bar2");
    inList2.add("baz2");
    List<String> inList3 = new ArrayList<String>();
    inList3.add("foo3");
    inList3.add("bar3");
    inList3.add("baz3");

    list.add(inList);
    list.add(inList2);
    list.add(inList3);
    System.err.println(list);

    // final List<String> outList= new ArrayList<String>();
    ExecutorService exec = Executors.newFixedThreadPool(2);
    try {
      for (final List<String> o : list) {
        exec.submit(
            new Runnable() {
              @Override
              public void run() {
                o.add("X");
              }
            });
      }
    } finally {
      exec.shutdown();
    }
    System.err.println(list);
  }

  @Test
  public void canGetIndexOfCharsOnFormattedLine() {
    // Unformatted: "  10mFOOBAR[0mfoobar]"
    // String fline= "  \033[38;5;10m10m\033[0mFOOBAR\033[48;5;10m[0mfoobar]\033[0m";
    String fline = "foo";

    // Reconstruct the string without codes
    String chars = "";
    for (int idx : Utils.indexOfCharsOnFormattedLine(fline)) {
      chars += fline.charAt(idx);
    }
    assertEquals("foo", chars);

    fline = "\033[0mFOO\033[38;5;10m  BAR";
    chars = "";
    for (int idx : Utils.indexOfCharsOnFormattedLine(fline)) {
      chars += fline.charAt(idx);
    }
    assertEquals("FOO  BAR", chars);

    // No printable chars:
    fline = "\033[0m";
    assertEquals(0, Utils.indexOfCharsOnFormattedLine(fline).size());
    fline = "";
    assertEquals(0, Utils.indexOfCharsOnFormattedLine(fline).size());
  }

  //    @Test
  //    public void canHighlightCenterColumn(){
  //        String profile= "AAANAAA\n"
  //                      + "AAA AAA\n"
  //                      + "TTTTT";
  //        String hProfile= Utils.highlightCenterColumn(profile);
  //        String exp= "AAA\033[27m;7mNAAA\n";
  //    }

  @Test
  public void canParseStringToCoords() throws InvalidGenomicCoordsException {

    assertEquals(
        Arrays.asList(new String[] {"chr1", "1", "10"}),
        Utils.parseStringCoordsToList("chr1:1-10"));

    assertEquals(
        Arrays.asList(new String[] {"chr1", "1", "536870912"}),
        Utils.parseStringCoordsToList("chr1"));

    assertEquals(
        Arrays.asList(new String[] {"chr1", "10", "10"}), Utils.parseStringCoordsToList("chr1:10"));

    // Can handle : in chrom name
    assertEquals(
        Arrays.asList(new String[] {"foo:bar", "10", "100"}),
        Utils.parseStringCoordsToList("foo:bar:10-100"));

    assertEquals(
        Arrays.asList(new String[] {"foo:bar", "10", "10"}),
        Utils.parseStringCoordsToList("foo:bar:10"));

    // Can handle spaces and comma separators
    assertEquals(
        Arrays.asList(new String[] {"chr1", "1000", "2000"}),
        Utils.parseStringCoordsToList("chr1 : 1,000 - 2,000"));

    // Chrom name with : and without from-to part
    assertEquals(
        Arrays.asList(new String[] {"foo:bar", "1", "536870912"}),
        Utils.parseStringCoordsToList("foo:bar"));

    // Chrom name with : and without from-to part
    assertEquals(
        Arrays.asList(new String[] {"foo:bar", "1", "1"}),
        Utils.parseStringCoordsToList("foo:bar:1"));

    // Chrom name with ':'
    assertEquals(
        Arrays.asList(new String[] {"foo:bar:1", "10", "10"}),
        Utils.parseStringCoordsToList("foo:bar:1:10"));

    assertEquals(
        Arrays.asList(new String[] {"chr1", "55681590", "55681890"}),
        Utils.parseStringCoordsToList("chr1:55681590-55681890"));

    // Invalid strings:
    boolean pass = false;
    try {
      Utils.parseStringCoordsToList("chr1:0-10");
    } catch (InvalidGenomicCoordsException e) {
      pass = true;
    }
    assertTrue(pass);

    pass = false;
    try {
      Utils.parseStringCoordsToList("chr1:10-5");
    } catch (InvalidGenomicCoordsException e) {
      pass = true;
    }
    assertTrue(pass);

    pass = false;
    try {
      Utils.parseStringCoordsToList("chr1:1-536870913"); // To large span: > 2^29 (536870912)
    } catch (InvalidGenomicCoordsException e) {
      pass = true;
    }
    assertTrue(pass);
    // Max span:
    assertEquals(
        Arrays.asList(new String[] {"chr1", "1", "536870912"}),
        Utils.parseStringCoordsToList("chr1:1-536870912"));
  }

  @Test
  public void canGetCommandFlag() {
    String[] cmd = {"-r", "foo", "bar"};
    List<String> argList = new LinkedList<String>(Arrays.asList(cmd));
    assertTrue(Utils.argListContainsFlag(argList, "foo"));
    assertTrue(!argList.contains("foo")); // Arg has been removed

    // Flag is not present
    argList = new LinkedList<String>(Arrays.asList(cmd));
    assertTrue(!Utils.argListContainsFlag(argList, "spam"));
    assertEquals(cmd.length, argList.size()); // No change
  }

  @Test
  public void canGetCommandWithMultipleArgs() throws InvalidCommandLineException {
    String[] cmd = {"-baz", "-r", "foo", "bar"};
    List<String> argList = new LinkedList<String>(Arrays.asList(cmd));

    List<String> args = new ArrayList<String>();
    args.add("foo");
    args.add("bar");

    assertEquals(args, Utils.getNArgsForParam(argList, "-r", 2));
    assertEquals(1, argList.size()); // Item left in original list

    // Missing arg
    argList = new LinkedList<String>(Arrays.asList(cmd));
    assertEquals(null, Utils.getNArgsForParam(argList, "-z", 2));
    assertEquals(cmd.length, argList.size());

    // Asked for too many args
    argList = new LinkedList<String>(Arrays.asList(cmd));
    boolean passed = false;
    try {
      Utils.getNArgsForParam(argList, "-r", 3);
    } catch (InvalidCommandLineException e) {
      passed = true;
    }
    assertTrue(passed);
    assertEquals(cmd.length, argList.size()); // Input list left unchanged
  }

  @Test
  public void canGetCommandArg() throws InvalidCommandLineException {
    String[] cmd = {"-r", "foo", "-x", "-bar"};
    List<String> argList = new LinkedList<String>(Arrays.asList(cmd));

    assertEquals("foo", Utils.getArgForParam(argList, "-r", null));
    assertTrue(!argList.contains("-r")); // Arg has been removed
    assertTrue(!argList.contains("foo")); // Arg has been removed

    // Param is not present
    argList = new LinkedList<String>(Arrays.asList(cmd));
    assertNull(Utils.getArgForParam(argList, "-z", null));
    assertEquals(cmd.length, argList.size()); // No change

    // Default arg:
    argList = new LinkedList<String>(Arrays.asList(cmd));
    assertEquals("default", Utils.getArgForParam(argList, "-z", "default"));
    assertEquals(cmd.length, argList.size()); // No change

    // Miss-specified arg throws exception:
    argList = new LinkedList<String>(Arrays.asList(cmd));
    boolean pass = false;
    try {
      Utils.getArgForParam(argList, "-bar", null);
    } catch (InvalidCommandLineException e) {
      pass = true;
    }
    assertTrue(pass);

    // Param given more than once: Return first found:
    String[] cmd2 = {"-r", "first", "-X", "-r", "second"};
    argList = new LinkedList<String>(Arrays.asList(cmd2));
    assertEquals("first", Utils.getArgForParam(argList, "-r", null));
    assertEquals("second", Utils.getArgForParam(argList, "-r", null));
    assertEquals("-X", argList.get(0));
  }

  @Test
  public void canPadMultilineString() {

    assertEquals("foo  ", Utils.padEndMultiLine("foo", 5));

    // Empty string is returned as is.
    assertEquals("", Utils.padEndMultiLine("", 3));

    // One newline is expanded to TWO strings
    assertEquals("   \n   ", Utils.padEndMultiLine("\n", 3));

    // Note starting with emtpy line, which is going to be padded.
    String x = "\nfoo\n1234567890";
    String padded = Utils.padEndMultiLine(x, 5);

    String[] p = padded.split("\n");
    assertEquals("     ", p[0]);
    assertEquals("foo  ", p[1]);
    assertEquals("1234567890", p[2]);
  }

  @Test
  public void canFilterArrayUsingAwk() throws IOException {

    String[] in3 = {"chr1\t1\t100", "chr1\t10\t100", "chr1\t2\t100"};
    boolean[] results = Utils.passAwkFilter(in3, "-v VAR=5 '$2 > VAR && $1'");
    assertTrue(!results[0]);
    assertTrue(results[1]);
    assertTrue(!results[2]);

    String[] in = {"chr1\t1\t100", "chr1\t1\t100", "chr1\t2\t100"};
    results = Utils.passAwkFilter(in, "-v VAR=5 '$2 > VAR && $1'");
    assertTrue(!results[0]);
    assertTrue(!results[1]);
    assertTrue(!results[2]);

    String[] in4 = {"chr1\t10\t100", "chr1\t10\t100", "chr1\t2\t100"};
    results = Utils.passAwkFilter(in4, "-v VAR=5 '$2 > VAR && $1'");
    assertTrue(results[0]);
    assertTrue(results[1]);
    assertTrue(!results[2]);

    // Zero length
    String[] in2 = {};
    results = Utils.passAwkFilter(in2, "-v VAR=5 '$2 > VAR && $1'");
    assertEquals(0, results.length);
  }

  @Test
  public void canFilterUsingAwk() throws IOException {

    String[] in = {"chr1\t1\t100", "chr1\t1\t100", "chr1\t2\t100"};
    boolean[] results = Utils.passAwkFilter(in, "-v VAR=5 '$2 > VAR && $1'");
    System.err.println(results[0]);
    System.err.println(results[1]);
    System.err.println(results[2]);

    // Note single quotes around the awk script
    assertTrue(Utils.passAwkFilter(new String[] {"chr1\t10\t100"}, "-v VAR=5 '$2 > VAR && $1'")[0]);
    assertTrue(!Utils.passAwkFilter(new String[] {"'chr1\t10\t100"}, "-v VAR=50 '$2 > VAR'")[0]);
    assertTrue(Utils.passAwkFilter(new String[] {"'chr1\t10\t100"}, "'($3 - $2) > 50'")[0]);
    assertTrue(!Utils.passAwkFilter(new String[] {"'chr1\t10\t100"}, "'($3 - $2) > 500'")[0]);

    assertTrue(
        Utils.passAwkFilter(new String[] {"'chr1\t10\t100'"}, "")[
            0]); // Empty script equals to no filter.
    assertTrue(Utils.passAwkFilter(new String[] {"'chr1\t10\t100'"}, "  ")[0]);

    // Valid awk script but output is not empty and not equal to input.
    assertTrue(!Utils.passAwkFilter(new String[] {"'chr1\t10\t100'"}, "'{print $1}'")[0]);
  }

  @Test
  public void canHandleQuotesInAwkScript() throws IOException {
    // NB: A single \ needs to be \\ in Java strings

    // We use triple quotes to wrap the script to avoid an escaping hell
    assertTrue(Utils.passAwkFilter(new String[] {"chr'1"}, "'''$1 == \"chr'1\";'''")[0]);

    // Double quote. Note three '\': Two to represent a single \ and one to escape
    // the double quote in Java string
    assertTrue(Utils.passAwkFilter(new String[] {"chr\"1"}, "'''$1 == \"chr\\\"1\";'''")[0]);

    // In v < 1.16 the awk script was passed as command line string instead of using -f <file>
    // and the mess below was required to match single quote.
    // See https://www.gnu.org/software/gawk/manual/html_node/Quoting.html
    // and see
    // http://stackoverflow.com/questions/9899001/how-to-escape-single-quote-in-awk-inside-printf
    assertTrue(Utils.passAwkFilter(new String[] {"chr'1"}, "'''$1 ~ \"chr\\x27\"'''")[0]);
  }

  @Test
  public void canQuoteString() {
    assertEquals(2, new Tokenizer("foo bar").tokenize().size()); // Just ensure tokenozer works!
    assertEquals(1, new Tokenizer(Utils.quote("foo bar")).tokenize().size());
    assertEquals(1, new Tokenizer(Utils.quote("'foo bar'")).tokenize().size());
    assertEquals(1, new Tokenizer(Utils.quote("\"foo bar\"")).tokenize().size());
    assertEquals(1, new Tokenizer(Utils.quote("foo'bar")).tokenize().size());
  }

  @Test
  public void canHandleEscapeInAwkScript() throws IOException {
    // NB: A single \ needs to be \\ in Java strings
    assertTrue(Utils.passAwkFilter(new String[] {"chr\\1"}, "'$1 == \"chr\\\\1\"'")[0]);
  }

  @Test
  public void canDetectBrokenAwkScript() throws IOException, InterruptedException {

    for (int i = 0; i < 20; i++) {
      // Broken awk script:
      boolean pass = false;
      try {
        assertTrue(!Utils.passAwkFilter(new String[] {"'chr1\t10\t100'"}, "'print {'")[0]);
      } catch (IOException e) {
        pass = true;
      }
      assertTrue(pass);
      Thread.sleep(3000);
    }
  }

  @Test
  public void canFilterGtfAwkFunc() throws IOException {
    // NB: You need to set the field separator to \\t.

    String[] gtf = {
      "GL873520\tchr1\tstop_codon\t8064\t8066\t0.000000\t-\t.\tgene_id 100; trax_id \"ACTB\";"
    };
    assertTrue(Utils.passAwkFilter(gtf, "-F '\\t' 'getGtfTag(\"gene_id\") == 100'")[0]);
    assertTrue(Utils.passAwkFilter(gtf, "-F '\\t' 'getGtfTag(\"trax_id\") == \"ACTB\"'")[0]);

    // Empty string if key not found
    assertTrue(Utils.passAwkFilter(gtf, "-F '\\t' 'getGtfTag(\"SPAM\") == \"\"'")[0]);

    // No attributes at all: Empty string returned
    gtf = new String[] {"GL873520\tchr1\tstop_codon\t8064\t8066\t0.000000\t-\t.\t."};
    assertTrue(Utils.passAwkFilter(gtf, "-F '\\t' 'getGtfTag(\"gene_id\") == \"\"'")[0]);

    // No attribute column at all (this would be an invalid GTF, by the way)
    gtf = new String[] {"GL873520\tchr1\tstop_codon\t8064\t8066\t0.000000\t-\t."};
    assertTrue(Utils.passAwkFilter(gtf, "-F '\\t' 'getGtfTag(\"gene_id\") == \"\"'")[0]);
  }

  @Test
  public void canFilterGffAwkFunc() throws IOException {
    // NB: You need to set the field separator to \\t.

    String[] x = {
      ".|.|.|.|.|.|.|.|Tag=100; ID = foo : bar ; Alias=spam,bar;".replaceAll("\\|", "\t")
    };
    assertTrue(Utils.passAwkFilter(x, "-F '\\t' 'getGffTag(\"Tag\") == 100'")[0]);
    assertTrue(Utils.passAwkFilter(x, "-F '\\t' 'getGffTag(\"ID\") == \"foo : bar\"'")[0]);
    assertTrue(Utils.passAwkFilter(x, "-F '\\t' 'getGffTag(\"Alias\") == \"spam,bar\"'")[0]);
    assertTrue(Utils.passAwkFilter(x, "-F '\\t' 'getGffTag(\"Alias\", 1) == \"spam\"'")[0]);
    assertTrue(Utils.passAwkFilter(x, "-F '\\t' 'getGffTag(\"Alias\", 2) == \"bar\"'")[0]);
    assertTrue(Utils.passAwkFilter(x, "-F '\\t' 'getGffTag(\"Alias\", 99) == \"\"'")[0]);
    assertTrue(Utils.passAwkFilter(x, "-F '\\t' 'getGffTag(\"SPAM\") == \"\"'")[0]);

    // NB: Missing tag i.e., empty string, evaluates to 0!!
    assertTrue(Utils.passAwkFilter(x, "-F '\\t' 'getGffTag(\"SPAM\") == 0'")[0]);

    x = new String[] {".|.|.|.|.|.|.|.|.".replaceAll("\\|", "\t")};
    assertTrue(Utils.passAwkFilter(x, "-F '\\t' 'getGffTag(\"Alias\") == \"\"'")[0]);

    x = new String[] {".|.|.|.|.|.|.|.".replaceAll("\\|", "\t")};
    assertTrue(Utils.passAwkFilter(x, "-F '\\t' 'getGffTag(\"Alias\") == \"\"'")[0]);

    // The value of the gff tag includes double-quotes. So it is "X" not X
    x =
        new String[] {
          ".|.|.|.|.|.|.|.|Tag=\"X\"".replaceAll("\\|", "\t")
        }; // Double quotes are not stripped
    assertTrue(Utils.passAwkFilter(x, "-F '\\t' 'getGffTag(\"Tag\") == \"\\\"X\\\"\"'")[0]);
  }

  @Test
  public void canFilterInfoVcfWithAwkFunc() throws IOException {
    String[] vcf = {
      "chr1 75888 . A T . . IMPRECISE;SVTYPE=DEL;DP=20,30;SVLEN=-32945;FOLD_CHANGE=0.723657;FOLD_CHANGE_LOG=-0.466623;PROBES=21"
          .replaceAll(" ", "\t")
    };
    assertTrue(Utils.passAwkFilter(vcf, "'getInfoTag(\"IMPRECISE\") == 1'")[0]);
    assertTrue(Utils.passAwkFilter(vcf, "'getInfoTag(\"ABSENT_TAG\") == 0'")[0]);
    assertTrue(Utils.passAwkFilter(vcf, "'getInfoTag(\"FOLD_CHANGE\") <= 0.723657'")[0]);
    assertTrue(!Utils.passAwkFilter(vcf, "'getInfoTag(\"FOLD_CHANGE\") > 0.723657'")[0]);

    // Split list of values
    assertTrue(Utils.passAwkFilter(vcf, "'getInfoTag(\"DP\") == \"20,30\"'")[0]);
    assertTrue(Utils.passAwkFilter(vcf, "'getInfoTag(\"DP\", 1) == 20'")[0]);
    assertTrue(Utils.passAwkFilter(vcf, "'getInfoTag(\"DP\", 2) == 30'")[0]);
    // Out of range
    assertTrue(Utils.passAwkFilter(vcf, "'getInfoTag(\"DP\", 3) == \"\"'")[0]);

    // Missing INFO
    vcf = new String[] {"chr1 75888 . A T . . .".replaceAll(" ", "\t")};
    assertTrue(Utils.passAwkFilter(vcf, "'getInfoTag(\"FOO\") == 0'")[0]);

    // INFO column not there at all
    vcf = new String[] {"chr1 75888 . A T . .".replaceAll(" ", "\t")};
    assertTrue(Utils.passAwkFilter(vcf, "'getInfoTag(\"FOO\") == 0'")[0]);
  }

  @Test
  public void canFilterFormatVcfWithAwkFunc() throws IOException {
    String[] vcf = {"chr1 75888 . A T . . . GT:GQ 11:21,10 22:99,100".replaceAll(" ", "\t")};
    assertTrue(Utils.passAwkFilter(vcf, "'getFmtTag(\"GT\") == 11'")[0]); // Default sample_idx= 1
    assertTrue(Utils.passAwkFilter(vcf, "'getFmtTag(\"GT\", 1) == 11'")[0]);
    assertTrue(Utils.passAwkFilter(vcf, "'getFmtTag(\"GT\", 2) == 22'")[0]);
    assertTrue(!Utils.passAwkFilter(vcf, "'getFmtTag(\"GT\", 2) < 22'")[0]);

    // Get value from list
    assertTrue(Utils.passAwkFilter(vcf, "'getFmtTag(\"GQ\", 2) == \"99,100\"'")[0]);
    assertTrue(Utils.passAwkFilter(vcf, "'getFmtTag(\"GQ\", 2, 1) == \"99\"'")[0]);
    assertTrue(Utils.passAwkFilter(vcf, "'getFmtTag(\"GQ\", 2, 2) == \"100\"'")[0]);
    // Out range
    assertTrue(Utils.passAwkFilter(vcf, "'getFmtTag(\"GQ\", 2, 3) == \"\"'")[0]);

    // Tag not found
    assertTrue(Utils.passAwkFilter(vcf, "'getFmtTag(\"ABSENT\", 1) == \"\"'")[0]);

    // Invalid indexes
    assertTrue(Utils.passAwkFilter(vcf, "'getFmtTag(\"GT\", 99) == \"\"'")[0]);
    assertTrue(Utils.passAwkFilter(vcf, "'getFmtTag(\"GT\", \"foobar\") == \"\"'")[0]);
  }

  //    @Test
  //    public void canFilterSamStructVarWithAwk() throws IOException{
  //        String rec= "r1 2113 chr1 951  60 50M50H = 1801 851  * *
  // SA:Z:chr1,1951,-,50M50S,60,0".replaceAll(" ", "\t");
  //        rec= "r1 81 chr1 951 60 50M50H = 1801 851 * * SA:Z:chr1,1951,-,50M50S,60,0".replaceAll("
  // ", "\t");
  //
  //        assertTrue(Utils.passAwkFilter(rec, "'isSV() > 0'"));
  //
  //        long t0= System.currentTimeMillis();
  //        int i= 0;
  //        while(i < 1000){
  //            Utils.passAwkFilter(rec, "'getSamTag(\"NH\") > 0'");
  //            i++;
  //        }
  //        long t1= System.currentTimeMillis();
  //        assertTrue((t1-t0) < 3000); // It can filter reasonably fast (?)
  //    }

  @Test
  public void canFilterSamTagWithAwk() throws IOException {
    String[] rec = {
      "read\t0\tchr7\t5566778\t50\t5M\t*\t0\t0\tCTCAT\tIIIII\tMD:Z:75\tRG:Z:1\tXG:i:0\tNH:i:1"
          + "\tNM:i:0\tXM:i:0\tXN:i:0\tXO:i:0\tAS:i:0\tYT:Z:UU"
    };

    // Filter for NH tag value
    assertTrue(Utils.passAwkFilter(rec, "'getSamTag(\"NH\") > 0'")[0]);
    assertFalse(Utils.passAwkFilter(rec, "'getSamTag(\"NH\") > 10'")[0]);

    // Missing tag
    assertFalse(Utils.passAwkFilter(rec, "'getSamTag(\"ZZ\") > 0'")[0]);

    // Missing tag searched but not used
    assertTrue(Utils.passAwkFilter(rec, "'{getSamTag(\"ZZ\"); print $0}'")[0]);
    assertFalse(Utils.passAwkFilter(rec, "'getSamTag(\"NM\") > 0'")[0]);

    // Tags missing altogether returns empty string
    rec = new String[] {"read\t0\tchr7\t5566778\t50\t5M\t*\t0\t0\tCTCAT\tIIIII"};
    assertTrue(Utils.passAwkFilter(rec, "'getSamTag(\"NM\") == \"\"'")[0]);

    long t0 = System.currentTimeMillis();
    int i = 0;
    while (i < 1000) {
      Utils.passAwkFilter(rec, "'getSamTag(\"NH\") > 0'");
      i++;
    }
    long t1 = System.currentTimeMillis();
    assertTrue((t1 - t0) < 10000); // It can filter reasonably fast (?)
  }

  /*
  @Test
  public void canFilterSamAlnEnd() throws IOException{
      String[] rec;

      // Missing CIGAR string
      rec = new String[] {"read 0 chr7 1 255 *".replaceAll(" ", "\t")};
      assertTrue(Utils.passAwkFilter(rec, "'getAlnEnd() == -1'")[0]);

      rec= new String[] {"read 0 chr7 1 255 5M".replaceAll(" ", "\t")};
      assertTrue(Utils.passAwkFilter(rec, "'getAlnEnd() == 5'")[0]);

      rec= new String[] {"read 0 chr7 10 255 5M".replaceAll(" ", "\t")};
      assertTrue(Utils.passAwkFilter(rec, "'getAlnEnd() == 14'")[0]);

      rec= new String[] {"read 0 chr7 1 255 50M".replaceAll(" ", "\t")};
      assertTrue(Utils.passAwkFilter(rec, "'getAlnEnd() == 50'")[0]);

      rec= new String[] {"read 0 chr7 1 255 9H50M10I5X3S".replaceAll(" ", "\t")};
      assertTrue(Utils.passAwkFilter(rec, "'getAlnEnd() == 55'")[0]);

      rec= new String[] {"read 0 chr7 1 255 1S2H3M4I5D6N7P8=9X99S33H".replaceAll(" ", "\t")};
      assertTrue(Utils.passAwkFilter(rec, "'getAlnEnd() == 31'")[0]);
  } */

  @Test
  public void canAwkFilterSamAln()
      throws IOException,
          InvalidGenomicCoordsException,
          ClassNotFoundException,
          InvalidRecordException,
          SQLException,
          InvalidColourException,
          InvalidCommandLineException {

    GenomicCoords gc = new GenomicCoords("chr7:1-100", 100, null, null);
    TrackReads tr = new TrackReads("test_data/pairs.sam", gc);
    tr.setNoFormat(true);
    tr.setAwk("'getAlnEnd() == 5'");
    assertEquals("NANAN\nTTTTT", tr.printToScreen().replaceAll(" ", ""));

    tr.setAwk("'getAlnEnd() == 24'");
    assertEquals("ntntn", tr.printToScreen().replaceAll(" ", ""));

    tr.setAwk("'getAlnLen() == 3'");
    assertEquals("NNN", tr.printToScreen().replaceAll(" ", ""));

    /*
    String[] rec;

    String cigar = "read 0 chr7 1 255 ";
    for(int i=0; i < 1000; i++) {
        cigar = cigar.concat("1M");
    }

    rec= new String[] {cigar.replaceAll(" ", "\t")};
    long t0;
    long t1;

    t0 = System.currentTimeMillis();
    for(int i = 0; i < 100; i++) {
        Utils.passAwkFilter(rec, "'getAlnEnd() > 0'");
    }
    t1 = System.currentTimeMillis();
    System.out.println("getAlnLen: " + (t1 - t0));


    t0 = System.currentTimeMillis();
    for(int i = 0; i < 100; i++) {
        Utils.passAwkFilter(rec, "'$1 != $2'");
    }
    t1 = System.currentTimeMillis();
    System.out.println("SIMPLE: " + (t1 - t0));
    */
  }

  /*
  @Test
  public void canFilterSamAlnLen() throws IOException{
      String[] rec;

      // Missing CIGAR string
      rec = new String[] {"read 0 chr7 1 255 *".replaceAll(" ", "\t")};
      assertTrue(Utils.passAwkFilter(rec, "'getAlnLen() == -1'")[0]);

      rec = new String[] {"read 0 chr7 1 255 1M".replaceAll(" ", "\t")};
      assertTrue(Utils.passAwkFilter(rec, "'getAlnLen() == 1'")[0]);

      rec = new String[] {"read 0 chr7 1 255 10M".replaceAll(" ", "\t")};
      assertTrue(Utils.passAwkFilter(rec, "'getAlnLen() == 10'")[0]);

      rec = new String[] {"read 0 chr7 5 255 10M".replaceAll(" ", "\t")};
      assertTrue(Utils.passAwkFilter(rec, "'getAlnLen() == 10'")[0]);

      rec = new String[] {"read 0 chr7 5 255 10M10S".replaceAll(" ", "\t")};
      assertTrue(Utils.passAwkFilter(rec, "'getAlnLen() == 10'")[0]);
  }
  */

  @Test
  public void canFilterSamBitflagWithAwk() throws IOException {

    assertTrue(Utils.passAwkFilter(new String[] {"read\t16"}, "'bitset(16) == 1'")[0]);
    assertTrue(Utils.passAwkFilter(new String[] {"read\t17"}, "'bitset(16) == 1'")[0]);
    assertTrue(Utils.passAwkFilter(new String[] {"read\t1"}, "'bitset(16) == 0'")[0]);
    assertTrue(Utils.passAwkFilter(new String[] {"read\t3585"}, "'bitset(1025) == 1'")[0]);
    assertTrue(Utils.passAwkFilter(new String[] {"read\t3584"}, "'bitset(1025) == 0'")[0]);
    assertTrue(Utils.passAwkFilter(new String[] {"read\t3840"}, "'bitset(3840) == 1'")[0]);
    assertTrue(Utils.passAwkFilter(new String[] {"read\t3840"}, "'bitset(3841) == 0'")[0]);

    assertTrue(Utils.passAwkFilter(new String[] {"read\t3840"}, "'bitset(0) == 1'")[0]);
    assertTrue(Utils.passAwkFilter(new String[] {"read\t0"}, "'bitset(0) == 1'")[0]);

    // Invalid input: Always return 0 (false)
    assertTrue(Utils.passAwkFilter(new String[] {"read\tfoo"}, "'bitset(0) == 0'")[0]);
    assertTrue(Utils.passAwkFilter(new String[] {"read\t16"}, "'bitset(\"foo\") == 0'")[0]);
    assertTrue(Utils.passAwkFilter(new String[] {"read\t16"}, "'bitset(-16) == 0'")[0]);
    assertTrue(Utils.passAwkFilter(new String[] {"read\t-16"}, "'bitset(0) == 0'")[0]);
    assertTrue(Utils.passAwkFilter(new String[] {"read\t16"}, "'bitset(17) == 0'")[0]);
  }

  @Test
  public void canCheckForUpdates() throws IOException {

    try {
      List<String> up = Utils.checkUpdates(50000);

      System.err.println(up);

      assertEquals(2, up.size());
      assertTrue(Character.isDigit(up.get(0).charAt(0)));
      assertTrue(Character.isDigit(up.get(1).charAt(0)));

      assertEquals(0, Utils.versionCompare("1.0.0", "1.0.0"));
      assertEquals(1, Utils.versionCompare("1.0.0", "0.0.9")); // Running version ahead of repo

      assertEquals(-1, Utils.versionCompare("1.0.0", "1.0.1")); // Running version out of date
      assertEquals(-1, Utils.versionCompare("1.0.0", "1.0.0.1"));

      // This should timeput and throw a warning
      up = Utils.checkUpdates(1);
    } catch (IOException e) {
      if (e.getMessage().equals("Server returned HTTP response code: 403 for URL")) {
        // This happen when travis runs from github. It would be nice to fix it.
      }
    }
  }

  @Test
  public void testAwk() throws Exception {

    String fin = "test_data/hg19.gencode_genes_v19.gtf.gz";
    InputStream is = new FileInputStream(fin);

    String[] args = {"-F", "\t", "NR <= 100000"};
    System.err.println(Arrays.toString(args));

    PrintStream stdout = System.out;
    ByteArrayOutputStream baos = new ByteArrayOutputStream();
    try {
      PrintStream os = new PrintStream(baos);
      System.err.println("START");
      long t0 = System.currentTimeMillis();
      new org.jawk.Main(args, is, os, System.err);
      long t1 = System.currentTimeMillis();
      System.err.println(t1 - t0);

    } finally {
      System.setOut(stdout);
      is.close();
    }
    System.err.println("DONE");
    long t0 = System.currentTimeMillis();
    String content = new String(baos.toByteArray(), StandardCharsets.UTF_8);
    // assertEquals(10, content.split("\n").length);

    List<String> full = FileUtils.readLines(new File(fin), "UTF-8");
    Set<String> keepme = new HashSet<String>(Arrays.asList(content.split("\n")));
    for (String line : full) {
      if (keepme.contains(line)) {
        // System.err.println(line);
      }
    }
    long t1 = System.currentTimeMillis();
    System.err.println(t1 - t0);
  }

  @Test
  public void testHistory() throws IOException {

    ConsoleReader console = new ConsoleReader();
    // History history= new History(new File(System.getProperty("user.home") + File.separator +
    // ".asciigenome_history"));
    History history = new MemoryHistory();
    history.add("foobar");
    history.add("baz");
    console.setHistory(history);
    System.out.println(console.getHistory());
    console.close();
  }

  @Test
  public void canGlobFiles() throws IOException {

    ArrayList<String> cmdInput = Utils.tokenize("test_data/ear*{bam,tdf} README.*", " ");
    List<String> globbed = Utils.globFiles(cmdInput);
    assertEquals(3, globbed.size());

    cmdInput = Utils.tokenize("test_data", " ");
    globbed = Utils.globFiles(cmdInput);
    assertTrue(globbed.size() > 10);

    cmdInput = Utils.tokenize("test_data/*", " ");
    globbed = Utils.globFiles(cmdInput);
    assertTrue(globbed.size() > 10);

    cmdInput = Utils.tokenize("test_data//*", " ");
    globbed = Utils.globFiles(cmdInput);
    assertTrue(globbed.size() > 10);

    cmdInput = Utils.tokenize("test_data/../test_data", " ");
    globbed = Utils.globFiles(cmdInput);
    assertTrue(globbed.size() > 10);

    cmdInput = Utils.tokenize("test_data/*.gtf.*", " ");
    globbed = Utils.globFiles(cmdInput);
  }

  @Test
  public void canGetRangeOfListOfValues() {
    List<Float> y = new ArrayList<Float>();
    y.add((float) 1.0);
    y.add((float) 10.0);
    y.add((float) 3.0);
    y.add(Float.NaN);
    assertEquals(1.0, Utils.range(y)[0], 0.0001); // Min
    assertEquals(10.0, Utils.range(y)[1], 0.0001); // Max

    // Only NaN
    List<Float> nan = new ArrayList<Float>();
    nan.add(Float.NaN);
    nan.add(Float.NaN);
    nan.add(Float.NaN);
    assertTrue(Utils.range(nan)[0].isNaN());
    assertTrue(Utils.range(nan)[1].isNaN());

    // Length of one
    List<Float> y1 = new ArrayList<Float>();
    y1.add((float) 1.0);
    assertEquals(1.0, Utils.range(y1)[0], 0.0001);
    assertEquals(1.0, Utils.range(y1)[1], 0.0001);

    // Zero length
    List<Float> y0 = new ArrayList<Float>();
    assertTrue(Utils.range(y0)[0].isNaN());
    assertTrue(Utils.range(y0)[1].isNaN());
  }

  @Test
  public void testSortByValueReverse() {

    Map<Character, Integer> baseCount = new LinkedHashMap<Character, Integer>();
    baseCount.put('A', 0);
    baseCount.put('C', 0);
    baseCount.put('G', 0);
    baseCount.put('T', 0);
    baseCount.put('N', 0);

    for (int i = 0; i < 10000; i++) {
      int count = baseCount.get('G') + 1;
      baseCount.put('G', count);
    }
    assertEquals('G', (char) Utils.sortByValue(baseCount).keySet().iterator().next());
  }

  @Test
  public void createFastaIndex() throws IOException, UnindexableFastaFileException {
    String fastaFile = "test_data/noindex.fa";
    if (new File("test_data/noindex.fa.fai").isFile()) {
      new File("test_data/noindex.fa.fai").delete();
    }

    Utils.checkFasta(fastaFile, 0);

    assertTrue((new File("test_data/noindex.fa.fai")).isFile());
    assertTrue((new File("test_data/noindex.fa.fai")).length() > 10);

    boolean pass = false;
    try {
      Utils.checkFasta("foo.bar", 2);
    } catch (IOException e) {
      pass = true;
    }
    assertTrue(pass);
  }

  @Test
  public void canMergeIntervals() throws InvalidGenomicCoordsException, InvalidColourException {

    // Zero len list
    List<IntervalFeature> intv = new ArrayList<IntervalFeature>();
    assertEquals(0, Utils.mergeIntervalFeatures(intv, false).size());

    // Fully contained feature
    intv.clear();
    intv.add(
        new IntervalFeature(
            "chr1 . . 100 1000 . . .".replaceAll(" ", "\t"), TrackFormat.GTF, null, -1));
    intv.add(
        new IntervalFeature(
            "chr1 . . 200 300 . . .".replaceAll(" ", "\t"), TrackFormat.GTF, null, -1));
    assertEquals(1, Utils.mergeIntervalFeatures(intv, false).size());
    assertEquals(100, Utils.mergeIntervalFeatures(intv, false).get(0).getFrom());
    assertEquals(1000, Utils.mergeIntervalFeatures(intv, false).get(0).getTo());

    // Partial overlap contained feature
    intv.clear();
    intv.add(
        new IntervalFeature(
            "chr1 . . 100 1000 . . .".replaceAll(" ", "\t"), TrackFormat.GTF, null, -1));
    intv.add(
        new IntervalFeature(
            "chr1 . . 200 300 . . .".replaceAll(" ", "\t"), TrackFormat.GTF, null, -1));
    intv.add(
        new IntervalFeature(
            "chr1 . . 500 5000 . . .".replaceAll(" ", "\t"), TrackFormat.GTF, null, -1));
    assertEquals(1, Utils.mergeIntervalFeatures(intv, false).size());
    assertEquals(100, Utils.mergeIntervalFeatures(intv, false).get(0).getFrom());
    assertEquals(5000, Utils.mergeIntervalFeatures(intv, false).get(0).getTo());

    /* MEMO: Start of bed features must be augmented by 1 */
    // One feature
    intv.clear();
    intv.add(new IntervalFeature("chr1 0 10 x1".replaceAll(" ", "\t"), TrackFormat.BED, null, -1));
    assertEquals(1, Utils.mergeIntervalFeatures(intv, false).get(0).getFrom());
    // Test the name is taken from the original feature since only one interval is merged (i.e. no
    // merging at all)
    assertEquals(intv.get(0).getName(), Utils.mergeIntervalFeatures(intv, false).get(0).getName());

    // One feature overalapping
    intv.add(new IntervalFeature("chr1 5 10".replaceAll(" ", "\t"), TrackFormat.BED, null, -1));
    IntervalFeature expected =
        new IntervalFeature("chr1 0 10".replaceAll(" ", "\t"), TrackFormat.BED, null, -1);

    assertEquals(expected.getFrom(), Utils.mergeIntervalFeatures(intv, false).get(0).getFrom());
    assertTrue(expected.equals(Utils.mergeIntervalFeatures(intv, false).get(0)));

    intv.add(new IntervalFeature("chr1 20 100".replaceAll(" ", "\t"), TrackFormat.BED, null, -1));
    assertEquals(2, Utils.mergeIntervalFeatures(intv, false).size());
    assertEquals(21, Utils.mergeIntervalFeatures(intv, false).get(1).getFrom());
    assertEquals(100, Utils.mergeIntervalFeatures(intv, false).get(1).getTo());

    intv.add(new IntervalFeature("chr1 30 110".replaceAll(" ", "\t"), TrackFormat.BED, null, -1));
    intv.add(new IntervalFeature("chr1 50 110".replaceAll(" ", "\t"), TrackFormat.BED, null, -1));
    assertEquals(2, Utils.mergeIntervalFeatures(intv, false).size());
    assertEquals(21, Utils.mergeIntervalFeatures(intv, false).get(1).getFrom());
    assertEquals(110, Utils.mergeIntervalFeatures(intv, false).get(1).getTo());

    // Touching features get merged into a single one
    intv.clear();
    intv.add(new IntervalFeature("chr1 0 10".replaceAll(" ", "\t"), TrackFormat.BED, null, -1));
    intv.add(new IntervalFeature("chr1 10 20".replaceAll(" ", "\t"), TrackFormat.BED, null, -1));
    assertEquals(1, Utils.mergeIntervalFeatures(intv, false).size());
    assertEquals(1, Utils.mergeIntervalFeatures(intv, false).get(0).getFrom());
    assertEquals(20, Utils.mergeIntervalFeatures(intv, false).get(0).getTo());

    // Touching GFF feature
    intv.clear();
    intv.add(
        new IntervalFeature(
            "chr1 . . 1 10 . . .".replaceAll(" ", "\t"), TrackFormat.GTF, null, -1));
    intv.add(
        new IntervalFeature(
            "chr1 . . 11 20 . . .".replaceAll(" ", "\t"), TrackFormat.GTF, null, -1));
    assertEquals(1, Utils.mergeIntervalFeatures(intv, false).size());
    assertEquals(1, Utils.mergeIntervalFeatures(intv, false).get(0).getFrom());
    assertEquals(20, Utils.mergeIntervalFeatures(intv, false).get(0).getTo());

    // Nothing to merge
    intv.clear();
    intv.add(
        new IntervalFeature(
            "chr1 . . 1 10 . . .".replaceAll(" ", "\t"), TrackFormat.GTF, null, -1));
    intv.add(
        new IntervalFeature(
            "chr1 . . 20 30 . . .".replaceAll(" ", "\t"), TrackFormat.GTF, null, -1));
    intv.add(
        new IntervalFeature(
            "chr1 . . 40 50 . . .".replaceAll(" ", "\t"), TrackFormat.GTF, null, -1));
    assertEquals(3, Utils.mergeIntervalFeatures(intv, false).size());

    intv.add(
        new IntervalFeature(
            "chr1 . . 40 50 . . .".replaceAll(" ", "\t"), TrackFormat.GTF, null, -1));
    assertEquals(3, Utils.mergeIntervalFeatures(intv, false).size());
  }

  @Test
  public void testStringContainsRegex() {
    String x = "foobarbaz";
    String regex = "b.r";
    System.out.println("PATTERN:" + Pattern.compile(regex).matcher(x).find());
  }

  @Test
  public void canParseGoToRegion() {
    assertEquals("1-1000", Utils.parseGoToRegion("1-1000"));
    assertEquals("1-1000", Utils.parseGoToRegion("1 - 1,000  "));
    assertEquals("1000", Utils.parseGoToRegion("1,000  "));
    assertEquals("1000-1400", Utils.parseGoToRegion("1000 1200 1300 1400"));
    assertEquals("1000-1400", Utils.parseGoToRegion("1k 1200 1300 1.4k"));
    assertEquals("1000-1400000", Utils.parseGoToRegion(" **1k*******1.4M** "));
    assertEquals("1000-5000000", Utils.parseGoToRegion("--**1k*******1.4M**------5.0M--"));
  }

  @Test
  public void canParseZoomArg() throws InvalidCommandLineException {
    assertEquals(2, Utils.parseZoom("zo 2", 1));
    assertEquals(3, Utils.parseZoom("zo 3 foo", 1));
    assertEquals(4, Utils.parseZoom("zo   4  ", 1));
    assertEquals(0, Utils.parseZoom("zo   0", 1));
    assertEquals(1, Utils.parseZoom("zo", 1));

    // Invalid input:
    boolean pass = false;
    try {
      Utils.parseZoom("zo -3", 1);
    } catch (InvalidCommandLineException e) {
      pass = true;
    }

    assertTrue(pass);
    pass = false;
    try {
      Utils.parseZoom("zo foo", 1);
    } catch (NumberFormatException e) {
      pass = true;
    }

    assertTrue(pass);
    pass = false;
    try {
      Utils.parseZoom("zo 3.3", 1);
    } catch (NumberFormatException e) {
      pass = true;
    }
    assertTrue(pass);
  }

  @Test
  public void canTabulateListOfFeatures() {
    List<String> rawList = new ArrayList<String>();
    rawList.add("1\tgenedb\tgene\t2964\t45090");
    rawList.add("chr1\tgenedb_long\tgene\t2964\t45090");
    rawList.add("1\tfoo\tna\t2"); // Missing last field
    rawList.add("1\tfoo\tna\t2\t10");
    rawList.add("1\tfoo\t\t2\t10"); // Empty cell

    List<String> expList = new ArrayList<String>();
    expList.add("1    genedb      gene 2964 45090");
    expList.add("chr1 genedb_long gene 2964 45090");
    expList.add("1    foo         na   2");
    expList.add("1    foo         na   2    10");
    expList.add("1    foo              2    10");

    List<String> obsList = Utils.tabulateList(rawList, -1, " ");
    assertThat(obsList, is(obsList));

    // Flush left "foo" rows as it has too much white space
    obsList =
        Utils.tabulateList(
            rawList, 5 * 4, " "); // *4 has to be adjusted according to what you have in the code.
    System.err.println(Joiner.on("\n").join(obsList));
    assertTrue(Joiner.on("\n").join(obsList).contains("1   "));
    assertTrue(Joiner.on("\n").join(obsList).contains("foo na"));

    obsList = Utils.tabulateList(rawList, 5 * 4, " | ");
    assertTrue(Joiner.on("\n").join(obsList).contains("1    | "));
    assertTrue(Joiner.on("\n").join(obsList).contains("foo | na"));

    // Handling region with no features to print
    rawList = new ArrayList<String>();
    obsList = Utils.tabulateList(rawList, -1, " ");
    expList = new ArrayList<String>();
    assertThat(expList, is(obsList));
  }

  @Test
  public void canGetBamReadCount() throws IOException {
    assertEquals(15098, Utils.getAlignedReadCount("test_data/ds051.actb.bam"));
    // Painfully slow!
    // assertEquals(6337212,
    // Utils.getAlignedReadCount("http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeCaltechRnaSeq/wgEncodeCaltechRnaSeqGm12878R2x75Il400SplicesRep2V2.bam"));
  }

  @Test
  public void canCountReadsInWindow2() throws InvalidGenomicCoordsException, IOException {
    GenomicCoords gc = new GenomicCoords("chr7:5524838-5611878", 80, samSeqDict, fastaFile);
    List<SamRecordFilter> filters = new ArrayList<SamRecordFilter>();

    assertEquals(100377, Utils.countReadsInWindow("test_data/ear045.oxBS.actb.bam", gc, filters));

    filters.add(new AlignedFilter(true));
    assertEquals(100265, Utils.countReadsInWindow("test_data/ear045.oxBS.actb.bam", gc, filters));

    filters = new ArrayList<SamRecordFilter>();
    filters.add(new AlignedFilter(true));
    filters.add(new ReadNegativeStrandFilter(false));
    assertEquals(50157, Utils.countReadsInWindow("test_data/ear045.oxBS.actb.bam", gc, filters));

    filters = new ArrayList<SamRecordFilter>();
    filters.add(new AlignedFilter(true));
    filters.add(new ReadNegativeStrandFilter(true));
    assertEquals(50108, Utils.countReadsInWindow("test_data/ear045.oxBS.actb.bam", gc, filters));

    filters = FlagToFilter.flagToFilterList(80, 1026);
    assertEquals(2729, Utils.countReadsInWindow("test_data/ear045.oxBS.actb.bam", gc, filters));

    filters = FlagToFilter.flagToFilterList(80, 1026);
    filters.add(new MappingQualityFilter(30));
    assertEquals(1592, Utils.countReadsInWindow("test_data/ear045.oxBS.actb.bam", gc, filters));
  }

  @Test
  public void canCountReadsInWindow() throws InvalidGenomicCoordsException, IOException {
    GenomicCoords gc = new GenomicCoords("chr7:5522436-5613572", 80, samSeqDict, fastaFile);
    List<SamRecordFilter> filters = new ArrayList<SamRecordFilter>();

    filters.add(new MappingQualityFilter(30)); // Same as
    filters.add(new FirstOfPairFilter(true)); // samtools view -q 30 -f 64

    long t0 = System.currentTimeMillis();
    for (int i = 0; i < 10; i++) {
      assertEquals(42770, Utils.countReadsInWindow("test_data/ear045.oxBS.actb.bam", gc, filters));
    }
    long t1 = System.currentTimeMillis();
    System.out.println("TIME TO FILTER: " + (t1 - t0));

    gc = new GenomicCoords("chr7:5524838-5611878", 80, samSeqDict, fastaFile);
  }

  @Test
  public void canTestForTabixIndex() throws IOException {
    assertTrue(Utils.hasTabixIndex("test_data/test.bedGraph.gz"));
    assertTrue(!Utils.hasTabixIndex("test_data/test.bedGraph"));

    // This ftp file has index but see https://github.com/samtools/htsjdk/issues/797
    assertTrue(
        !Utils.hasTabixIndex(
            "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/analysis_results/input_call_sets/ALL.wex.union_illumina_wcmc_bcm_bc_bi.20110521.snps.exome.sites.vcf.gz"));

    // HTTP is ok.
    // If this file does not exist, put any valid tabix file and its index on Dropbox/Public and use
    // the dropbox link here.
    assertTrue(
        Utils.hasTabixIndex("http://genome.ucsc.edu/goldenPath/help/examples/vcfExample.vcf.gz"));

    // NB: Uncompressed files give a OutOfMemoryError: Java heap space
    assertTrue(
        !Utils.hasTabixIndex(
            "ftp://ftp.solgenomics.net/tomato_genome/annotation/ITAG3.0_release/ITAG3.0_RepeatModeler_repeats_light.gff"));
  }

  @Test
  public void canGetFileTypeFromName() {

    assertEquals(TrackFormat.VCF, Utils.getFileTypeFromName("test/gz.vcf.bgz"));

    assertEquals(TrackFormat.BIGWIG, Utils.getFileTypeFromName("http://foo/bar/wgEncode.bigWig"));

    assertEquals(TrackFormat.BAM, Utils.getFileTypeFromName("test/foo.bam"));

    assertEquals(TrackFormat.BAM, Utils.getFileTypeFromName("test/foo.SAM.GZ"));

    assertEquals(TrackFormat.BAM, Utils.getFileTypeFromName("test/foo.CRAM"));
  }

  @Test
  public void canInitRegionFromCram()
      throws IOException,
          InvalidGenomicCoordsException,
          ClassNotFoundException,
          InvalidCommandLineException,
          InvalidRecordException,
          SQLException {
    assertEquals(
        "chr7:5566778", Utils.initRegionFromFile("test_data/ds051.actb.cram", "test_data/chr7.fa"));
  }

  @Test
  public void canInitRegion()
      throws IOException,
          InvalidGenomicCoordsException,
          ClassNotFoundException,
          InvalidCommandLineException,
          InvalidRecordException,
          SQLException {

    // Files with no records
    assertEquals("chr1", Utils.initRegionFromFile("test_data/norecords.vcf"));
    assertEquals("Undefined_contig", Utils.initRegionFromFile("test_data/norecords_nodict.vcf"));
    assertEquals("chr7", Utils.initRegionFromFile("test_data/norecords.sam"));
    assertEquals("Undefined_contig", Utils.initRegionFromFile("test_data/empty.bedGraph"));
    assertEquals("Undefined_contig", Utils.initRegionFromFile("test_data/empty2.bedGraph"));
    // Note: empty bigBed and bigWig seems to fail for reasons independent of ASCIIGenome

    assertEquals("chr7:5566778", Utils.initRegionFromFile("test_data/ds051.short.bam"));
    assertEquals(
        "chr7:5566778",
        Utils.initRegionFromFile(
            "https://raw.githubusercontent.com/dariober/ASCIIGenome/master/test_data/ds051.short.bam"));
    assertEquals("chr9", Utils.initRegionFromFile("test_data/hg18_var_sample.wig.v2.1.30.tdf"));
    assertEquals(
        "chr9:15077",
        Utils.initRegionFromFile(
            "test_data/wgEncodeCaltechRnaSeqGm12878R2x75Il400SigRep2V2.sample.bigWig"));
    assertEquals("chr1:67208779", Utils.initRegionFromFile("test_data/refSeq.hg19.short.bed"));
    assertEquals(
        "chr1:8404074", Utils.initRegionFromFile("test_data/refSeq.hg19.short.sort.bed.gz"));
    assertEquals("chr1:11874", Utils.initRegionFromFile("test_data/hg19_genes_head.gtf.gz"));
    assertEquals(
        "chr9:41066",
        Utils.initRegionFromFile("test_data/wgEncodeDukeDnase8988T.fdr01peaks.hg19.bb"));
    assertEquals(
        "1:113054374", Utils.initRegionFromFile("test_data/CEU.exon.2010_06.genotypes.vcf"));

    boolean pass = false;
    try {
      Utils.initRegionFromFile(
          "http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeCaltechRnaSeq/wgEncodeCaltechRnaSeqGm12878R2x75Il400SigRep2V2.bigWig");
    } catch (InvalidGenomicCoordsException e) {
      pass = true;
    }
    assertTrue(pass);
  }

  @Test
  public void canInitRegionWithMissingContig()
      throws IOException,
          InvalidGenomicCoordsException,
          ClassNotFoundException,
          InvalidCommandLineException,
          InvalidRecordException,
          SQLException {
    // Reproduce issue#86
    assertEquals("chr7", Utils.initRegionFromFile("test_data/missing_contig.sam"));
  }

  @Test
  public void testBamHasIndex() throws IOException {
    assertTrue(Utils.bamHasIndex("test_data/ds051.short.bam"));
    assertTrue(!Utils.bamHasIndex("test_data/ds051.noindex.bam"));
  }

  @Test
  public void canGetClosestIndex() {

    List<Double> seq = Utils.seqFromToLenOut(10, 50, 5);

    System.err.println(seq);
    assertEquals(2, Utils.getIndexOfclosestValue(30, seq));
    assertEquals(
        3,
        Utils.getIndexOfclosestValue(
            35,
            seq)); // Value is right in the middle of the interval (but take care you are working
    // with floating points)
    assertEquals(2, Utils.getIndexOfclosestValue(29, seq));
    assertEquals(2, Utils.getIndexOfclosestValue(31, seq));
    assertEquals(4, Utils.getIndexOfclosestValue(50, seq));
    assertEquals(0, Utils.getIndexOfclosestValue(3, seq)); // Genome pos is to the left of window
    // This should not happen though.

    int windowSize = 150;
    List<Double> mapping = Utils.seqFromToLenOut(1, 1000000, windowSize);
    for (int i = 0; i < windowSize; i++) {
      System.out.println("Index: " + i + " position: " + mapping.get(i));
    }

    assertEquals(
        windowSize - 1,
        Arrays.binarySearch(mapping.toArray(new Double[mapping.size()]), 1000000.0));

    assertEquals(windowSize - 1, Utils.getIndexOfclosestValue(1000000, mapping));
    assertEquals((windowSize / 2) - 1, Utils.getIndexOfclosestValue(1000000 / 2, mapping));
    assertEquals(0, Utils.getIndexOfclosestValue(1, mapping));
    assertEquals(0, Utils.getIndexOfclosestValue((6712 / 2.0) - 1, mapping));
    assertEquals(1, Utils.getIndexOfclosestValue((6712 / 2.0) + 1, mapping));
    assertEquals(0, Utils.getIndexOfclosestValue((6712 / 2.0), mapping));
    assertEquals(102, Utils.getIndexOfclosestValue(684564, mapping));
    assertEquals(112, Utils.getIndexOfclosestValue(750000, mapping));
  }

  @Test
  public void canGenerateSequence() {
    Utils.seqFromToLenOut(15, 17, 13);
    Utils.seqFromToLenOut(17, 15, 13);
    Utils.seqFromToLenOut(15, 26, 13).size();
    Utils.seqFromToLenOut(15, 15, 13);

    // Length of 1 returns "from" like R seq(0, 10, length.out= 1) -> 0
    assertEquals(1, Utils.seqFromToLenOut(0, 10, 1).size());
    assertEquals(0, Utils.seqFromToLenOut(0, 10, 1).get(0), 0.00001);
    assertEquals(0, Utils.seqFromToLenOut(0, 10, 0).size()); // Zero-length sequence

    assertEquals((Double) Double.NaN, Utils.seqFromToLenOut(Double.NaN, Double.NaN, 10).get(0));
  }

  @Test
  public void canTestForAllNaN() {
    ArrayList<Double> x = new ArrayList<Double>();
    x.add(Double.NaN);
    x.add(Double.NaN);
    x.add(Double.NaN);
    assertTrue(Utils.allIsNaN(x));

    x.add((Double) 1.0);
    assertFalse(Utils.allIsNaN(x));
  }

  @Test
  public void canParseInputAndUpdateGenomicCoords()
      throws InvalidGenomicCoordsException, IOException, InvalidCommandLineException {
    GenomicCoords gc = new GenomicCoords("chr7:100-200", 80, samSeqDict, fastaFile);
    List<String> tokens = new ArrayList<String>();
    tokens.add("+10");
    String region = Utils.parseConsoleInput(tokens, gc);
    assertEquals("chr7:110-210", region);

    tokens.set(0, "-100000");
    region = Utils.parseConsoleInput(tokens, gc);
    assertTrue(region.startsWith("chr7:1-"));

    tokens.set(0, "f");
    tokens.add(1, "1");
    String fregion = Utils.parseConsoleInput(tokens, gc);
    assertTrue(!region.equals(fregion));

    boolean pass = false;
    try {
      tokens.clear();
      tokens.add("f");
      tokens.add("foobar");
      Utils.parseConsoleInput(tokens, gc);
    } catch (InvalidCommandLineException e) {
      pass = true;
    }
    assertTrue(pass);

    pass = false;
    try {
      tokens.clear();
      tokens.add("b");
      tokens.add("foobar");
      Utils.parseConsoleInput(tokens, gc);
    } catch (InvalidCommandLineException e) {
      pass = true;
    }
    assertTrue(pass);
  }

  @Test
  public void canRoundNumbersToSignificantDigits() {

    Arrays.asList(Utils.roundToSignificantDigits(85477601.0, 85657825.0, 2));

    double x = 1000.123456789;
    double y = 1001.123456789;
    int nSignif = 3;
    String[] rounded = Utils.roundToSignificantDigits(x, y, nSignif);
    assertEquals("1000.123", rounded[0]); // Regular rounding
    assertEquals("1001.123", rounded[1]);

    x = 1000.00012345;
    y = 1000.0012345;
    nSignif = 2;
    rounded = Utils.roundToSignificantDigits(x, y, nSignif);
    assertEquals("1000.00012", rounded[0]);
    assertEquals("1000.00123", rounded[1]);

    x = 1000.0009876;
    y = 1000.00987654;
    nSignif = 3;
    rounded = Utils.roundToSignificantDigits(x, y, nSignif);
    assertEquals("1000.000988", rounded[0]);
    assertEquals("1000.009877", rounded[1]);
  }

  @Test
  public void canRoundToSiginficantDigits() {
    assertEquals(0.001234, Utils.roundToSignificantFigures(0.001234, 5), 1e-16);
    assertEquals(1200, (int) Utils.roundToSignificantFigures(1234, 2));
  }

  @Test
  public void canAddMetricSuffixToInt() {
    assertEquals("123M", Utils.parseIntToMetricSuffix(123000000));
    assertEquals("123k", Utils.parseIntToMetricSuffix(123000));
    assertEquals("123", Utils.parseIntToMetricSuffix(123));
  }

  @Test
  public void canAddTracksToList()
      throws IOException, InvalidGenomicCoordsException, InvalidCommandLineException {
    List<String> inputFileList = new ArrayList<String>();
    inputFileList.add("foo");
    inputFileList.add("bar");
    List<String> newFileNames = new ArrayList<String>();
    newFileNames.add("test_data/ds051.actb.bam");
    newFileNames.add("nonsense");

    boolean pass = false;
    try {
      Utils.addSourceName(inputFileList, newFileNames, 0);
    } catch (InvalidCommandLineException e) {
      pass = true;
    }
    assertTrue(pass);
    assertEquals(2, inputFileList.size());
  }

  //    @Test
  //    public void canAddTrackFromStdin() throws IOException, InvalidGenomicCoordsException,
  // InvalidCommandLineException{
  //        List<String> inputFileList= new ArrayList<String>();
  //        List<String> newFileNames= new ArrayList<String>();
  //        newFileNames.add("-");
  //
  //        this.simulateStdin("test_data/hg19_genes_head.gtf");
  //        Utils.addSourceName(inputFileList, newFileNames, 0);
  //        assertEquals(1, inputFileList.size());
  //        assertTrue(new File(inputFileList.get(0)).exists());
  //        assertTrue(inputFileList.get(0).endsWith(".gtf"));
  //    }

  //    @Test
  //    public void canSplitCommandLineInTokens(){
  //
  //        List<String> x=
  // Lists.newArrayList(Splitter.on(Pattern.compile("&&(?=([^']*'[^']*')*[^']*$)"))
  //            .trimResults()
  //            .split("awk '$4 > 910000 && $4 < 1000000' && grep -i foo"));
  //
  //        System.err.println(x);
  //
  //        List<String> cmdInput= Utils.tokenize("awk '$4 > 910000 && $4 < 1000000'", "&&");
  //        System.err.println(cmdInput);
  //        assertEquals(1, cmdInput.size());
  //
  //        // Note single quotes
  //        cmdInput= Utils.tokenize("foo && bar && 'baz && baz'", "&&");
  //        assertEquals(3, cmdInput.size());
  //        assertEquals("baz && baz", cmdInput.get(2));
  //    }

  @Test
  public void canSplitStringInTokens() {

    String str = "foo    bar    baz   ";
    assertTrue(Utils.tokenize(str, " ").contains("bar"));
    assertTrue(Utils.tokenize(str, " ").contains("baz"));
    // assertEquals("bar", Utils.tokenize(str, " ").get(1));

    // Note use of quotes
    ArrayList<String> xx =
        Utils.tokenize(
            "'foo && bar' " + "&& bar" + "&&baz " + "&& 'foo && biz'" + "&& 'foo && \' biz'", "&&");
    for (String token : xx) {
      System.out.println(token);
    }

    assertEquals("foo && bar", xx.get(0));
    assertEquals("bar", xx.get(1));
    assertEquals("baz", xx.get(2));
    assertEquals("foo && biz", xx.get(3));

    xx = Utils.tokenize("gene \"ACTB\"", "&&");
    assertEquals("gene \"ACTB\"", xx.get(0));

    // Unusual input:
    assertEquals(null, Utils.tokenize(null, " "));
    assertTrue(Utils.tokenize("", " ").size() == 0);
    //        assertTrue(Utils.tokenize("   ", " ").size() == 0);

    // Empty token:
    assertTrue(Utils.tokenize("foo -f '' -bar", " ").get(2).isEmpty());

    // Reverse token
    String cmdInput = "goto chr1 0 100";
    assertEquals(cmdInput, Joiner.on(" ").join(Utils.tokenize(cmdInput, " ")));
    assertEquals(
        "goto chr1 0 100", Joiner.on(" ").join(Utils.tokenize("   goto   chr1   0   100  ", " ")));
  }

  @Test
  public void canGetWritableFileOrNull() throws InvalidGenomicCoordsException, IOException {
    GenomicCoords gc = new GenomicCoords("chr7:1-200", 80, samSeqDict, fastaFile);
    String x = Utils.parseCmdinputToGetSnapshotFile("save", gc);
    assertEquals("chr7_1_200.txt", x);
    x = Utils.parseCmdinputToGetSnapshotFile("save /tmp/foo.txt", gc);
    assertEquals("/tmp/foo.txt", x);
  }

  @Test
  public void canConvertCoordsToString() {
    assertEquals("chr1:1-100", Utils.coordinatesToString("chr1", 1, 100));
    assertEquals("chr1:1", Utils.coordinatesToString("chr1", 1, null));
    assertEquals("chr1:1", Utils.coordinatesToString("chr1", 1, null));
    assertEquals("chr1:1", Utils.coordinatesToString("chr1", null, null));
    assertEquals("chr1:1", Utils.coordinatesToString("chr1", null, -1));
    assertEquals("chr1:10", Utils.coordinatesToString("chr1", 10, 9));
  }

  @Test
  public void canExpandTildeToHomeDir() {
    // No change
    assertEquals("/foo/bar/baz", Utils.tildeToHomeDir("/foo/bar/baz"));

    // Expand to home dir
    assertTrue(Utils.tildeToHomeDir("~/foo/bar/baz").startsWith(File.separator));

    // No change
    assertEquals("/~foo/~bar/baz", Utils.tildeToHomeDir("/~foo/~bar/baz"));

    // No change
    assertEquals("~foo", Utils.tildeToHomeDir("~foo"));

    // Expanded to /Users/berald01/
    assertEquals(System.getProperty("user.home") + File.separator, Utils.tildeToHomeDir("~/"));

    // This will fail most likelly fail on systems other than *nix:
    assertTrue(
        Utils.tildeToHomeDir("~/foo/bar/baz").startsWith("/Users/")
            || Utils.tildeToHomeDir("~/foo/bar/baz").startsWith("/home/"));
  }

  @Test
  public void canGetTemplateNameFromReadName() {
    assertEquals("keepme", Utils.templateNameFromSamReadName("keepme"));
    assertEquals("keepme", Utils.templateNameFromSamReadName("keepme "));
    assertEquals("keepme", Utils.templateNameFromSamReadName("keepme   foo bar"));
    assertEquals("keepme", Utils.templateNameFromSamReadName("keepme/1"));
    assertEquals("keepme", Utils.templateNameFromSamReadName("keepme/2"));
    assertEquals("keepme", Utils.templateNameFromSamReadName("keepme /2"));
    Stopwatch sw = Stopwatch.createStarted();
    String x = "HSQ9103:404:C6F0VANXX:1:2208:4363:50381 foo bar /1";
    for (int i = 0; i < 1000000; i++) {
      Utils.templateNameFromSamReadName(x);
    }
    System.err.println(sw.stop());
  }
}
