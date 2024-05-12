package tracks;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import coloring.Config;
import coloring.ConfigKey;
import coloring.Xterm256;
import exceptions.InvalidColourException;
import exceptions.InvalidCommandLineException;
import exceptions.InvalidConfigException;
import exceptions.InvalidGenomicCoordsException;
import exceptions.InvalidRecordException;
import filter.FlagToFilter;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.filter.SamRecordFilter;
import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Pattern;
import org.junit.BeforeClass;
import org.junit.Test;
import samTextViewer.GenomicCoords;

public class TrackReadsTest {

  @BeforeClass
  public static void init() throws IOException, InvalidConfigException {
    new Config(null);
  }

  static SamReaderFactory srf = SamReaderFactory.make();
  static SamReader samReader = srf.open(new File("test_data/ds051.actb.bam"));
  public static SAMSequenceDictionary samSeqDict =
      samReader.getFileHeader().getSequenceDictionary();
  public static String fastaFile = "test_data/chr7.fa";

  @Test
  public void canReturnChromosomeNames()
      throws ClassNotFoundException,
          IOException,
          InvalidGenomicCoordsException,
          InvalidRecordException,
          SQLException {
    GenomicCoords gc = new GenomicCoords("chr7", 80, samSeqDict, null);
    TrackReads tr = new TrackReads("test_data/ds051.actb.bam", gc);
    assertTrue(tr.getChromosomeNames().contains("chr1"));
    assertTrue(tr.getChromosomeNames().contains("chr7"));
  }

  @Test
  public void canReloadTrack()
      throws ClassNotFoundException,
          IOException,
          InvalidGenomicCoordsException,
          InvalidRecordException,
          SQLException {
    GenomicCoords gc = new GenomicCoords("chr7:1-200", 80, null, null);
    TrackReads tr = new TrackReads("test_data/pairs.sam", gc);
    tr.reload();
  }

  @Test
  public void canFilterTopStrand()
      throws ClassNotFoundException,
          IOException,
          InvalidGenomicCoordsException,
          InvalidRecordException,
          SQLException,
          InvalidColourException {
    GenomicCoords gc = new GenomicCoords("chr7:1-80", 80, null, null);
    TrackReads tr = new TrackReads("test_data/testTopBottomStrand.sam", gc);
    tr.setNoFormat(true);

    // Top strand
    int f = 4096;
    int F = 0;
    tr.set_f_flag(f);

    List<SamRecordFilter> filters = new ArrayList<SamRecordFilter>();
    filters.addAll(FlagToFilter.flagToFilterList(f, F));
    tr.setSamRecordFilter(filters);
    for (String x : tr.getRecordsAsStrings()) {
      assertTrue(x.startsWith("top"));
    }

    // Bottom strand
    f = 0;
    F = 4096;
    tr.set_f_flag(f);

    filters = new ArrayList<SamRecordFilter>();
    filters.addAll(FlagToFilter.flagToFilterList(f, F));
    tr.setSamRecordFilter(filters);
    for (String x : tr.getRecordsAsStrings()) {
      assertTrue(x.startsWith("bottom"));
    }
  }

  @Test
  public void canShadeLowBaseQuality()
      throws InvalidGenomicCoordsException,
          InvalidColourException,
          ClassNotFoundException,
          IOException,
          InvalidRecordException,
          SQLException,
          InvalidCommandLineException,
          InvalidConfigException {

    String shade = Config.get(ConfigKey.shade_low_mapq);
    String xshade = Integer.toString(Xterm256.colorNameToXterm256(shade));

    GenomicCoords gc = new GenomicCoords("chr7:999-1041", 80, null, null);
    TrackReads tr = new TrackReads("test_data/missingReadSeq.bam", gc);
    System.err.println(tr.printToScreen());
    assertTrue(tr.printToScreen().trim().startsWith("[48;5;" + xshade));

    // Read with soft clipped bases
    gc = new GenomicCoords("chr7:9999-10050", 80, null, null);
    tr = new TrackReads("test_data/missingReadSeq.bam", gc);
    assertTrue(tr.printToScreen().contains(xshade));

    // Read with deletions and skipped bases
    gc = new GenomicCoords("chr7:19999-20050", 80, null, null);
    tr = new TrackReads("test_data/missingReadSeq.bam", gc);
    tr.printToScreen();
  }

  // @Test // STUB
  public void canChangeReadColourOnRegex()
      throws InvalidGenomicCoordsException,
          IOException,
          ClassNotFoundException,
          InvalidRecordException,
          SQLException,
          InvalidColourException {
    GenomicCoords gc = new GenomicCoords("chr7:5566778-5566943", 80, null, null);
    TrackReads tr = new TrackReads("test_data/ds051.short.bam", gc);
    List<Argument> list = new ArrayList<Argument>();
    Argument re = new Argument("NCNNNCCC", "red1", false);
    list.add(re);
    tr.changeFeatureColor(list);
    assertTrue(tr.printToScreen().contains("196;"));
  }

  @Test
  public void canReadReadsWithMissingSequence()
      throws ClassNotFoundException,
          IOException,
          InvalidGenomicCoordsException,
          InvalidRecordException,
          SQLException,
          InvalidColourException {

    // Window size larger then 1 bp per column
    GenomicCoords gc = new GenomicCoords("chr7:1-1000", 80, null, null);
    TrackReads tr = new TrackReads("test_data/missingReadSeq.bam", gc);
    tr.setNoFormat(true);
    assertTrue(tr.printToScreen().trim().startsWith(">"));
    assertTrue(tr.printToScreen().trim().endsWith(">"));
    assertTrue(tr.printToScreen().trim().length() > 1);

    // Window size < 1 bp per column
    gc = new GenomicCoords("chr7:100-120", 80, null, null);
    tr = new TrackReads("test_data/missingReadSeq.bam", gc);
    tr.setNoFormat(true);
    assertTrue(tr.printToScreen().trim().startsWith("N"));
    assertTrue(tr.printToScreen().trim().endsWith("N"));
    assertTrue(tr.printToScreen().trim().length() > 1);
  }

  public void canExplainSamFlag()
      throws ClassNotFoundException,
          IOException,
          InvalidGenomicCoordsException,
          InvalidRecordException,
          SQLException,
          InvalidColourException,
          InvalidCommandLineException {

    GenomicCoords gc = new GenomicCoords("chr7:5529094-5529267", 80, null, null);
    TrackReads tr = new TrackReads("test_data/ear045.oxBS.actb.bam", gc);
    tr.setNoFormat(true);
    tr.setPrintMode(PrintRawLine.FULL);
    tr.setExplainSamFlag(true);
    String printable = tr.printLines();
    System.err.println(printable);
    assertTrue(printable.contains("163|+|pair|prop-p|mate-rev|2nd"));
  }

  @Test
  public void canShowReadsAsPairs() throws Exception {
    GenomicCoords gc = new GenomicCoords("chr7:1-80", 80, null, null);
    TrackReads tr = new TrackReads("test_data/pairs.sam", gc);
    tr.setNoFormat(true);
    // Still unpaired
    assertTrue(tr.printToScreen().startsWith("NANAN  gcgcgcgcgc  ntntn "));
    assertTrue(tr.printToScreen().contains("ATATATATAT  "));

    tr.setReadsAsPairs(true); // Switch on pairing

    // Properly paired
    assertTrue(tr.printToScreen().contains(" GGGGG~~~~~~~~~~~~~~~ggggg "));

    // Overlapping pair
    assertTrue(tr.printToScreen().contains("ATATATAgcgcgcgcgc "));
    System.err.println(tr.printToScreen());

    // Pair where one read is fully contained in the other
    assertTrue(tr.printToScreen().contains(" CCaaaaaaCC "));

    // Properly paired but mate is not in the window
    assertTrue(tr.printToScreen().contains("TTTTT  "));

    // Properly paired flag is set but mateAlignmentStart is not set. Treat as unpaired:
    assertTrue(tr.printToScreen().contains("NANAN  "));
    assertTrue(tr.printToScreen().contains(" ntntn "));
  }

  @Test
  public void canReturnReadsAsRawStrings()
      throws ClassNotFoundException,
          IOException,
          InvalidGenomicCoordsException,
          InvalidRecordException,
          SQLException {
    GenomicCoords gc = new GenomicCoords("chr7:5566000-5567000", 80, null, null);
    TrackReads tr = new TrackReads("test_data/ds051.short.bam", gc);

    assertEquals(22, tr.getRecordsAsStrings().size());

    // No read in interval
    gc = new GenomicCoords("chr1:1-1000", 80, null, null);
    tr = new TrackReads("test_data/ds051.short.bam", gc);
    assertEquals(0, tr.getRecordsAsStrings().size());
  }

  @Test
  public void canPrintReads()
      throws InvalidGenomicCoordsException,
          IOException,
          ClassNotFoundException,
          InvalidRecordException,
          SQLException,
          InvalidColourException,
          InvalidConfigException,
          InvalidCommandLineException {

    GenomicCoords gc = new GenomicCoords("chr7:5566000-5567000", 80, null, null);
    TrackReads tr = new TrackReads("test_data/ds051.short.bam", gc);
    tr.setNoFormat(true);

    tr.setPrintMode(PrintRawLine.OFF);
    String printable = tr.printLines();
    assertEquals("", printable); // Empty string if printing is off.

    tr.setPrintMode(PrintRawLine.CLIP);

    tr.setPrintRawLineCount(5);
    printable = tr.printLines();
    assertTrue(printable.length() > 100);
    assertEquals(
        5 + 1, printable.split("\n").length); // Expect 6 lines: 5 for reads and 1 for info header.
  }

  @Test
  public void canShowReadsInWindow() throws Exception {
    GenomicCoords gc = new GenomicCoords("chr7:5566000-5567000", 80, samSeqDict, null);
    TrackReads tr = new TrackReads("test_data/ds051.short.bam", gc);
    tr.setNoFormat(true);
    tr.setyMaxLines(1000);
    assertTrue(tr.getTitle().contains("22")); // N. reads stacked in this interval before filtering
  }

  @Test
  public void canFilterReadsContainingVariants()
      throws InvalidGenomicCoordsException,
          IOException,
          ClassNotFoundException,
          InvalidRecordException,
          SQLException,
          InvalidColourException {

    GenomicCoords gc = new GenomicCoords("chr7:1000001-1000081", 100, samSeqDict, fastaFile);
    TrackReads tr = new TrackReads("test_data/variant_reads.sam", gc);
    tr.setNoFormat(true);

    tr.setyMaxLines(1000);
    assertEquals(
        11,
        tr.printToScreen()
            .split("\n")
            .length); // N. reads stacked in this interval before filtering

    tr.setVariantReadInInterval("chr7", 1000001, 1000001, true);
    assertEquals(2, tr.printToScreen().split("\n").length);
    assertTrue(tr.printToScreen().startsWith("A"));

    tr.setVariantReadInInterval("chr7", 1000015, 1000015, true);
    assertEquals(1, tr.printToScreen().split("\n").length);
    assertTrue(tr.printToScreen().trim().endsWith("T"));

    tr.setVariantReadInInterval("chr7", 1000003, 1000003, true); // Test deletion
    assertEquals(1, tr.printToScreen().split("\n").length);
    assertTrue(tr.printToScreen().contains("-"));

    tr.setVariantReadInInterval("chr7", 1000005, 1000005, true); // Test deletion
    assertEquals(1, tr.printToScreen().split("\n").length);

    tr.setVariantReadInInterval("chr7", 1000002, 1000002, true); // Test insertion
    assertEquals(2, tr.printToScreen().split("\n").length);

    tr.setVariantReadInInterval("chr7", 1000014, 1000014, true); // No variants at this site
    assertEquals("", tr.printToScreen().trim());

    tr.setVariantReadInInterval("chr7", 1000031, 1000040, true); // No read
    assertEquals("", tr.printToScreen().trim());

    tr.setVariantReadInInterval("chr7", 2000000, 2000040, true); // No read
    assertEquals("", tr.printToScreen().trim());

    tr.setVariantReadInInterval("chr7", 1000003, 1000010, true); // Test range
    assertEquals(3, tr.printToScreen().split("\n").length);
    assertTrue(tr.printToScreen().contains("-"));
    assertTrue(tr.printToScreen().contains("AG"));

    tr.setVariantReadInInterval(
        Filter.DEFAULT_VARIANT_CHROM.getValue(), -1, -1, true); // Remove filter
    assertEquals(
        11,
        tr.printToScreen()
            .split("\n")
            .length); // N. reads stacked in this interval before filtering

    // Genomic window does not overlap the variant range, but n reads in this window do overlap.
    gc = new GenomicCoords("chr7:1000002-1000081", 100, samSeqDict, fastaFile);
    tr = new TrackReads("test_data/variant_reads.sam", gc);
    tr.setNoFormat(true);
    tr.setyMaxLines(1000);
    tr.setVariantReadInInterval("chr7", 1000001, 1000001, true);
    assertEquals(2, tr.printToScreen().split("\n").length);
  }

  @Test
  public void canFilterReadsInIntervalKeepAll()
      throws InvalidGenomicCoordsException,
          IOException,
          ClassNotFoundException,
          InvalidRecordException,
          SQLException,
          InvalidColourException {

    GenomicCoords gc = new GenomicCoords("chr7:1000001-1000081", 100, samSeqDict, fastaFile);
    TrackReads tr = new TrackReads("test_data/variant_reads.sam", gc);
    tr.setNoFormat(true);
    tr.setVariantReadInInterval("chr7", 1000028, 1000028, true);
    assertEquals(2, tr.printToScreen().split("\n").length); // Default: only variants

    tr.setVariantReadInInterval("chr7", 1000028, 1000028, false);
    System.err.println(tr.printToScreen());
    assertEquals(4, tr.printToScreen().split("\n").length);
  }

  @Test
  public void canFilterReadsWithAwk()
      throws InvalidGenomicCoordsException,
          IOException,
          ClassNotFoundException,
          InvalidRecordException,
          SQLException,
          InvalidColourException {
    GenomicCoords gc = new GenomicCoords("chr7:5566000-5567000", 80, samSeqDict, null);
    TrackReads tr = new TrackReads("test_data/ds051.short.bam", gc);
    tr.setNoFormat(true);
    tr.setyMaxLines(1000);
    assertEquals(
        22,
        tr.printToScreen()
            .split("\n")
            .length); // N. reads stacked in this interval before filtering
    tr.setAwk("'$1 ~ \"NCNNNCCC\"'");
    assertEquals(6, tr.printToScreen().split("\n").length);
    assertTrue(tr.getTitle().contains("awk"));
    System.err.println(tr.getTitle());
  }

  @Test
  public void canFilterReadsWithGrep()
      throws InvalidGenomicCoordsException,
          IOException,
          ClassNotFoundException,
          InvalidRecordException,
          SQLException,
          InvalidColourException {
    GenomicCoords gc = new GenomicCoords("chr7:5566000-5567000", 80, samSeqDict, null);
    TrackReads tr = new TrackReads("test_data/ds051.short.bam", gc);
    tr.setNoFormat(true);
    tr.setyMaxLines(1000);
    assertEquals(
        22,
        tr.printToScreen()
            .split("\n")
            .length); // N. reads stacked in this interval before filtering
    tr.setShowHideRegex(Pattern.compile("NCNNNCCC"), Pattern.compile("\\t5566779\\t"));
    assertEquals(4, tr.printToScreen().split("\n").length);
  }

  @Test
  public void canFilterReadsWithGrepAndAwk()
      throws InvalidGenomicCoordsException,
          IOException,
          ClassNotFoundException,
          InvalidRecordException,
          SQLException,
          InvalidColourException {
    GenomicCoords gc = new GenomicCoords("chr7:5566000-5567000", 80, samSeqDict, null);
    TrackReads tr = new TrackReads("test_data/ds051.short.bam", gc);
    tr.setNoFormat(true);
    tr.setyMaxLines(1000);
    assertEquals(
        22,
        tr.printToScreen()
            .split("\n")
            .length); // N. reads stacked in this interval before filtering
    tr.setShowHideRegex(
        Pattern.compile("NCNNNCCC"), Pattern.compile(Filter.DEFAULT_HIDE_REGEX.getValue()));
    tr.setAwk("'$4 != 5566779'");
    assertEquals(4, tr.printToScreen().split("\n").length);
  }

  @Test
  public void canShowReadCount() throws Exception {
    GenomicCoords gc = new GenomicCoords("chr7:5565600-5567600", 80, null, null);
    TrackReads tr = new TrackReads("test_data/ear045.oxBS.actb.bam", gc);

    // Same as: samtools view -F 4 -c ear045.oxBS.actb.bam chr7:5565600-5567600
    // AND excluding also reads fully soft-clipped
    assertTrue(tr.getTitle().contains("2436"));
  }

  @Test
  public void canGetTitle()
      throws InvalidGenomicCoordsException,
          InvalidColourException,
          IOException,
          ClassNotFoundException,
          InvalidRecordException,
          SQLException {
    String bam = "test_data/adjacent.bam";
    GenomicCoords gc = new GenomicCoords("chr7:1-100", 80, samSeqDict, null);
    TrackReads tr = new TrackReads(bam, gc);

    tr.setNoFormat(true);
    tr.setTrackTag("aln.bam#1");

    assertTrue(tr.getTitle().trim().startsWith("aln.bam#1"));
  }

  @Test
  public void canShowActiveFiltersInTitle()
      throws ClassNotFoundException,
          IOException,
          InvalidGenomicCoordsException,
          InvalidRecordException,
          SQLException,
          InvalidColourException {
    String bam = "test_data/adjacent.bam";
    GenomicCoords gc = new GenomicCoords("chr7:1-100", 80, null, fastaFile);
    TrackReads tr = new TrackReads(bam, gc);
    tr.setNoFormat(true);
    assertTrue(!tr.getTitle().contains("filters"));

    tr.setMapq(Integer.valueOf(Filter.DEFAULT_MAPQ.getValue()));
    assertTrue(!tr.getTitle().contains("filters"));

    tr.setMapq(10);
    assertTrue(tr.getTitle().contains("mapq"));

    tr.set_F_flag(1024);
    assertTrue(tr.getTitle().contains("flag"));

    tr.setShowHideRegex(Pattern.compile(".*"), Pattern.compile("foo"));
    assertTrue(tr.getTitle().contains("grep"));

    tr.setVariantReadInInterval("chr7", 10, 100, true);
    assertTrue(tr.getTitle().contains("var-read"));

    tr.setAwk("'2 > 1'");
    assertTrue(tr.getTitle().contains("awk"));
  }

  @Test
  public void testOneLineStack()
      throws InvalidGenomicCoordsException,
          IOException,
          ClassNotFoundException,
          InvalidRecordException,
          SQLException,
          InvalidColourException {
    String bam = "test_data/adjacent.bam";
    int yMaxLines = 10;
    boolean bs = false;
    boolean noFormat = true;
    GenomicCoords gc = new GenomicCoords("chr7:1-50", 80, samSeqDict, null);

    TrackReads tr = new TrackReads(bam, gc);
    tr.setBisulf(bs);
    tr.setNoFormat(noFormat);

    // NB: The success of this test depends on the screen width of eclipse
    //        String exp=
    //        "AAAAAAAAAA           GGGGGGGGGG TTTTTTTTTT \n"+
    //        "          CCCCCCCCCC                       ";
    tr.setyMaxLines(yMaxLines);
    assertEquals(2, tr.printToScreen().split("\n").length); // Two lines

    System.out.println(tr.printToScreen());

    gc = new GenomicCoords("chr7:1-100", 80, samSeqDict, null);
    tr = new TrackReads(bam, gc);
    tr.setBisulf(bs);
    tr.setNoFormat(noFormat);
    tr.setyMaxLines(yMaxLines);
    assertTrue(tr.printToScreen().split("\n").length == 1);
  }

  @Test
  public void testNoReadsInRegion()
      throws InvalidGenomicCoordsException,
          IOException,
          ClassNotFoundException,
          InvalidRecordException,
          SQLException,
          InvalidColourException {
    String bam = "test_data/ds051.actb.bam";
    int yMaxLines = 50;
    boolean bs = false;
    boolean noFormat = true;
    GenomicCoords gc = new GenomicCoords("chr7:1-101", 80, samSeqDict, fastaFile);

    TrackReads tr = new TrackReads(bam, gc);
    tr.setyMaxLines(yMaxLines);
    tr.setBisulf(bs);
    tr.setNoFormat(noFormat);
    assertEquals("", tr.printToScreen());
  }

  @Test
  public void canResetToZeroLargeWindow()
      throws IOException,
          InvalidGenomicCoordsException,
          ClassNotFoundException,
          InvalidRecordException,
          SQLException,
          InvalidColourException {
    // If the genomic window is too large do not process the bam file and return zero height track.
    GenomicCoords gc = new GenomicCoords("chr7:1-100000000", 80, samSeqDict, fastaFile);
    TrackReads tr = new TrackReads("test_data/ds051.actb.bam", gc);
    assertEquals("", tr.printToScreen());
  }

  @Test
  public void canConstructFromUnsortedInput()
      throws IOException,
          InvalidGenomicCoordsException,
          ClassNotFoundException,
          InvalidRecordException,
          SQLException,
          InvalidColourException {
    // If the genomic window is too large do not process the bam file and return zero height track.
    GenomicCoords gc = new GenomicCoords("chr7:1-100000000", 80, samSeqDict, fastaFile);
    TrackReads tr = new TrackReads("test_data/ds051.noindex.sam", gc);
    assertEquals("", tr.printToScreen());
  }

  @Test
  public void canConstructFromUnsortedCram()
      throws ClassNotFoundException,
          IOException,
          InvalidGenomicCoordsException,
          InvalidRecordException,
          SQLException {
    GenomicCoords gc = new GenomicCoords("chr1:1-1000", 80, null, "test_data/chr7.fa");
    new TrackReads("test_data/ds051.noindex.cram", gc);
  }
}
