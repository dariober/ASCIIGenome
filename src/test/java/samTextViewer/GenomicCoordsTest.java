package samTextViewer;

import static org.junit.Assert.*;

import colouring.Config;
import com.google.common.base.Splitter;
import exceptions.InvalidColourException;
import exceptions.InvalidCommandLineException;
import exceptions.InvalidConfigException;
import exceptions.InvalidGenomicCoordsException;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.regex.Pattern;
import org.junit.Before;
import org.junit.Test;

public class GenomicCoordsTest {

  static SamReaderFactory srf = SamReaderFactory.make();
  static SamReader samReader = srf.open(new File("test_data/ds051.short.bam"));
  public static SAMSequenceDictionary samSeqDict =
      samReader.getFileHeader().getSequenceDictionary();

  public static String fastaFile = "test_data/chr7.fa";

  @Before
  public void initConfig() throws IOException, InvalidConfigException {
    new Config(null);
  }

  @Test
  public void canPrintSequenceDictionaryOnlyContigs()
      throws InvalidGenomicCoordsException, IOException {
    ArrayList<String> knownContigs = new ArrayList<String>();
    knownContigs.add("chr3");
    knownContigs.add("chr1");
    knownContigs.add("chr2");

    GenomicCoords gc = new GenomicCoords("chr3:1000-100000", 80, null, null);

    String out =
        gc.printSequenceDictionary(
            knownContigs, -1, -1, ".", ContigOrder.ALPHANUMERIC_ASC, 30, -1, true);
    List<String> outList = Splitter.on("\n").splitToList(out);
    assertTrue(outList.get(0).equals("Genome size: n/a; Number of contigs: 3"));
    assertTrue(outList.get(1).startsWith("chr1"));
  }

  @Test
  public void canPrintSequenceDictionaryBoldCurrentChrom()
      throws InvalidGenomicCoordsException, IOException {
    ArrayList<String> knownContigs = new ArrayList<String>();
    knownContigs.add("chr3");
    knownContigs.add("chr1");
    knownContigs.add("chr2");

    GenomicCoords gc = new GenomicCoords("chr3:1000-100000", 80, null, null);

    String out =
        gc.printSequenceDictionary(
            knownContigs, -1, -1, ".", ContigOrder.ALPHANUMERIC_ASC, 30, -1, false);
    assertTrue(out.contains("\033[1m" + "chr3"));

    out =
        gc.printSequenceDictionary(
            knownContigs, -1, -1, ".", ContigOrder.ALPHANUMERIC_ASC, 30, -1, true);
    assertTrue(Pattern.compile(".*chr3 <==.*").matcher(out).find());
  }

  @Test
  public void canPrintSequenceDictionary() throws InvalidGenomicCoordsException, IOException {
    GenomicCoords gc = new GenomicCoords("chrM:1000-2000", 80, samSeqDict, null);

    String out = gc.printSequenceDictionary(null, -1, -1, ".", ContigOrder.SIZE_DESC, 30, -1, true);
    List<String> outList = Splitter.on("\n").splitToList(out);

    assertTrue(outList.get(0).startsWith("Genome size: 3095693983; Number of contigs: 25"));
    assertTrue(outList.get(1).startsWith("chr1"));

    out =
        gc.printSequenceDictionary(null, -1, -1, ".", ContigOrder.ALPHANUMERIC_DESC, 30, -1, true);
    outList = Splitter.on("\n").splitToList(out);
    assertTrue(outList.get(1).startsWith("chrY"));

    out = gc.printSequenceDictionary(null, -1, -1, ".", ContigOrder.SIZE_DESC, 30, 3, false);
    outList = Splitter.on("\n").splitToList(out);
    assertTrue(outList.get(0).startsWith("Genome size: 3095693983; Number of contigs: 25"));
    assertEquals(5, outList.size());

    out = gc.printSequenceDictionary(null, -1, -1, ".", ContigOrder.ALPHANUMERIC_ASC, 30, -1, true);
    outList = Splitter.on("\n").splitToList(out);
    assertTrue(outList.get(0).startsWith("Genome size: 3095693983; Number of contigs: 25"));
    assertEquals(outList.size(), 26);
  }

  @Test
  public void canReturnNextChromNoDict() throws InvalidGenomicCoordsException, IOException {

    ArrayList<String> knownContigs = new ArrayList<String>();
    knownContigs.add("chr3");
    knownContigs.add("chr1");
    knownContigs.add("chr2");

    GenomicCoords gc = new GenomicCoords("chr3:1000-100000", 80, null, null);

    gc.nextChrom(knownContigs, -1, -1, ".*", ContigOrder.UNSORTED);
    assertEquals("chr1", gc.getChrom());

    // No effect of size limits
    gc = new GenomicCoords("chr1:1000-100000", 80, null, null);
    gc.nextChrom(knownContigs, 10, 100, ".*", ContigOrder.UNSORTED);
    assertEquals("chr2", gc.getChrom());
  }

  @Test
  public void canReturnNextChromSize() throws InvalidGenomicCoordsException, IOException {
    GenomicCoords gc = new GenomicCoords("chrM:1000-2000", 80, samSeqDict, null);

    // Next by size is...:
    gc.nextChrom(null, -1, -1, ".*", ContigOrder.SIZE_ASC);
    assertEquals("chr21", gc.getChrom());

    // Smaller by size is...:
    gc.nextChrom(null, -1, -1, ".*", ContigOrder.SIZE_DESC);
    assertEquals("chrM", gc.getChrom());

    // Next in dict is:...
    gc.nextChrom(null, -1, -1, ".*", ContigOrder.UNSORTED);
    assertEquals("chr1", gc.getChrom());

    // Wrap around: if you are on the largest chrom, go to the smallest:
    gc = new GenomicCoords("chr1:1000-2000", 80, samSeqDict, null);
    gc.nextChrom(null, -1, -1, ".*", ContigOrder.SIZE_ASC);
    assertEquals("chrM", gc.getChrom());
  }

  @Test
  public void canReturnNextChromAlphasort() throws InvalidGenomicCoordsException, IOException {
    GenomicCoords gc = new GenomicCoords("chrM:1000-2000", 80, samSeqDict, null);

    // Next by in alphanumeric order is...:
    gc.nextChrom(null, -1, -1, ".*", ContigOrder.ALPHANUMERIC_ASC);
    assertEquals("chrX", gc.getChrom());

    // and back
    gc.nextChrom(null, -1, -1, ".*", ContigOrder.ALPHANUMERIC_DESC);
    assertEquals("chrM", gc.getChrom());
  }

  @Test
  public void canReturnNextChromRegex() throws InvalidGenomicCoordsException, IOException {

    GenomicCoords gc = new GenomicCoords("chr1:1000-100000", 80, samSeqDict, null);
    gc.nextChrom(null, -1, -1, ".*", ContigOrder.UNSORTED);
    assertEquals("chr2", gc.getChrom());
    assertEquals((Integer) 1, gc.getFrom());
    assertEquals((Integer) 80, gc.getTo());

    gc.nextChrom(null, -1, -1, "Y$", ContigOrder.UNSORTED);
    assertEquals("chrY", gc.getChrom());

    boolean pass = false;
    try {
      gc.nextChrom(null, -1, -1, "chrY", ContigOrder.UNSORTED);
    } catch (InvalidGenomicCoordsException e) {
      assertEquals("chrY", gc.getChrom());
      pass = true;
    }
    assertTrue(pass);

    // Wrap around
    gc.nextChrom(null, -1, -1, ".", ContigOrder.UNSORTED);
    assertEquals("chrM", gc.getChrom());

    // Min and Max
    gc.nextChrom(null, 249000000, -1, ".", ContigOrder.UNSORTED);
    assertEquals("chr1", gc.getChrom());

    gc.nextChrom(null, -1, 20000, ".", ContigOrder.UNSORTED);
    assertEquals("chrM", gc.getChrom());
  }

  @Test
  public void canPrintPerecentRuler()
      throws InvalidGenomicCoordsException,
          IOException,
          InvalidConfigException,
          InvalidColourException {

    GenomicCoords gc = new GenomicCoords("chr1:101-2000", 80, samSeqDict, null);
    assertTrue(gc.printablePercentRuler(10, true).length() > 10);
    assertTrue((gc.printablePercentRuler(10, false).contains("[")));

    // Single digit % is handled correctly. I.e. 6 -> .06
    gc = new GenomicCoords("chr1:101-2000", 180, samSeqDict, null);
    assertTrue(gc.printablePercentRuler(10, true).startsWith("0         .06       .11"));
  }

  @Test
  public void canExtendCoordinates()
      throws InvalidGenomicCoordsException, IOException, InvalidCommandLineException {

    GenomicCoords gc = new GenomicCoords("chr7:1001-1009", 80, null, null);
    List<String> cmdInput = new ArrayList<String>();
    cmdInput.add("extend");
    cmdInput.add("20");
    cmdInput.add("30");
    cmdInput.add("mid");
    gc.cmdInputExtend(cmdInput);
    assertEquals(1005 - 20, (int) gc.getFrom());
    assertEquals(1005 + 30 - 1, (int) gc.getTo());

    // Check left doesn't go negative
    gc = new GenomicCoords("chr7:1-9", 80, samSeqDict, null);
    gc.cmdInputExtend(cmdInput);
    assertEquals(1, (int) gc.getFrom());
    assertEquals(5 + 30 - 1, (int) gc.getTo());

    // Check right doesn't go beyond fasta ref
    gc = new GenomicCoords("seq:101-109", 80, null, "test_data/seq_cg.fa");
    cmdInput.set(1, "50");
    gc.cmdInputExtend(cmdInput);
    assertEquals(120, (int) gc.getTo());

    // Check window mode
    gc = new GenomicCoords("seq:100-110", 80, null, null);
    cmdInput.set(1, "50");
    cmdInput.set(2, "60");
    cmdInput.set(3, "window");
    gc.cmdInputExtend(cmdInput);
    assertEquals(100 - 50, (int) gc.getFrom());
    assertEquals(110 + 60, (int) gc.getTo());

    // Check odd input: From coords after to swaps left with right.
    gc = new GenomicCoords("seq:100-200", 80, null, null);
    cmdInput.set(1, "-200");
    cmdInput.set(2, "0");
    cmdInput.set(3, "window");
    gc.cmdInputExtend(cmdInput);
    assertEquals(200, (int) gc.getFrom());
    assertEquals(300, (int) gc.getTo());
  }

  @Test
  public void canGetStringRegion() throws InvalidGenomicCoordsException, IOException {
    assertEquals("chr7:1-100", new GenomicCoords("chr7:1-100", 80, null, null).toStringRegion());

    GenomicCoords gc = new GenomicCoords("chr7:1", 80, null, null);
    assertTrue(
        gc.toStringRegion().matches("^chr7:1-\\d+$")); // String region is something like chr7:1-79

    gc = new GenomicCoords("chr7", 80, null, null);
    assertTrue(gc.toStringRegion().matches("^chr7:1-\\d+$")); // chr7:1-79
  }

  @Test
  public void canPrintChromMap()
      throws InvalidGenomicCoordsException,
          IOException,
          InvalidColourException,
          InvalidConfigException {

    GenomicCoords gc = new GenomicCoords("chr7:1-100", 80, samSeqDict, null);

    gc.getChromIdeogram(10, false);

    String chromMap = gc.getChromIdeogram(10, true);
    assertTrue(chromMap.startsWith("*---------"));

    gc = new GenomicCoords("chr7:1-100", 80, samSeqDict, null);
    chromMap = gc.getChromIdeogram(10, true);
    assertEquals(80, chromMap.length()); // 79 Should be the size of eclipse's terminal

    gc = new GenomicCoords("chr7:1-1000000000", 80, samSeqDict, null);
    chromMap = gc.getChromIdeogram(10, true);
    assertTrue(chromMap.startsWith("**********"));

    gc = new GenomicCoords("chr7:200000000-300000000", 80, samSeqDict, null);
    chromMap = gc.getChromIdeogram(10, true);
    assertTrue(chromMap.startsWith("1--------"));
  }

  @Test
  public void printRefSeq()
      throws InvalidGenomicCoordsException, IOException, InvalidColourException {
    GenomicCoords gc = new GenomicCoords("chr7:5540580-5540590", 80, null, "test_data/chr7.fa");
    assertEquals("ggccggctggg\n", gc.printableRefSeq(true));
  }

  @Test
  public void printRefSeqBgzip()
          throws InvalidGenomicCoordsException, IOException, InvalidColourException {
    GenomicCoords gc = new GenomicCoords("chr7:5540580-5540590", 80, null, "test_data/chr7.fa.gz");
    assertEquals("ggccggctggg\n", gc.printableRefSeq(true));
  }

  @Test
  public void canTestForEqualCoords() throws InvalidGenomicCoordsException, IOException {
    GenomicCoords gc = new GenomicCoords("chr1:1-10", 80, null, null);
    GenomicCoords other = new GenomicCoords("chr1:1-10", 80, null, null);
    assertTrue(gc.equalCoords(other));

    GenomicCoords other2 = new GenomicCoords("chr2:1-10", 80, null, null);
    assertTrue(!gc.equalCoords(other2));
    other2 = new GenomicCoords("chr1:2-10", 80, null, null);
    assertTrue(!gc.equalCoords(other2));
    other2 = new GenomicCoords("chr1:1-100", 80, null, null);
    assertTrue(!gc.equalCoords(other2));

    GenomicCoords gc2 = (GenomicCoords) gc.clone();
    System.out.println(gc);
    System.out.println(gc2);
    gc.zoomOut();
    System.out.println(gc);
    System.out.println(gc2);
    // assertTrue(gc2.equalCoords(gc));
  }

  @Test
  public void canInitializeSamSeqDictFromVCF() throws IOException, InvalidGenomicCoordsException {
    GenomicCoords gc = new GenomicCoords("1", 80, null, null);
    gc.setGenome(
        Arrays.asList(new String[] {"test_data/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz"}),
        true);
    assertEquals(25, gc.getSamSeqDict().size());

    gc = new GenomicCoords("1", 80, null, null);
    gc.setGenome(
        Arrays.asList(
            new String[] {
              "https://raw.githubusercontent.com/dariober/ASCIIGenome/master/test_data/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz"
            }),
        true);
    assertEquals(25, gc.getSamSeqDict().size());

    // Not indexed
    gc = new GenomicCoords("1", 80, null, null);
    gc.setGenome(
        Arrays.asList(new String[] {"test_data/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf"}),
        true);
    assertEquals(25, gc.getSamSeqDict().size());

    // VCF w/o sequence dict
    gc = new GenomicCoords("1", 80, null, null);
    gc.setGenome(Arrays.asList(new String[] {"CEU.exon.2010_06.genotypes.vcf.gz"}), true);
    assertNull(gc.getSamSeqDict());
  }

  @Test
  public void canInitializeSamSeqDictFromGenomeFile()
      throws IOException, InvalidGenomicCoordsException {

    List<String> insam = new ArrayList<String>();
    insam.add("hg19");

    // From resource:
    GenomicCoords gc = new GenomicCoords("chr1", 80, null, null);
    gc.setGenome(insam, true);
    assertEquals(93, gc.getSamSeqDict().size());

    // From bam header:
    insam.set(0, "test_data/ds051.short.bam");
    gc.setGenome(insam, true);
    assertEquals(25, gc.getSamSeqDict().size());

    // Check we get the full path to source file.
    assertTrue(gc.getSamSeqDictSource().length() > "test_data/ds051.short.bam".length());

    // From single item as string
    gc = new GenomicCoords("chr1", 80, null, null);
    gc.setGenome(Arrays.asList(new String[] {"test_data/ds051.short.bam"}), true);
    assertEquals(25, gc.getSamSeqDict().size());
  }

  @Test
  public void canInitializeSamSeqDictFromFasta() throws InvalidGenomicCoordsException, IOException {

    // From fasta
    GenomicCoords gc = new GenomicCoords("chr1", 80, null, null);
    gc.setGenome(Arrays.asList(new String[] {"test_data/chr7.fa"}), true);
    assertEquals(1, gc.getSamSeqDict().size());

    // From fasta without index
    gc = new GenomicCoords("chr1", 80, null, null);
    gc.setGenome(Arrays.asList(new String[] {"test_data/noindex.fa"}), false);
    assertEquals(1, gc.getSamSeqDict().size());
  }

  @Test
  public void initializeFromInvalidInput() throws InvalidGenomicCoordsException, IOException {
    // What if invalid input?
    GenomicCoords gc = new GenomicCoords("chr1", 80, null, null);
    assertEquals(null, gc.getSamSeqDict());

    // Non existent file or genome tag
    gc.setGenome(Arrays.asList(new String[] {"test_data/foo.fa"}), true);
    assertEquals(null, gc.getSamSeqDict());

    // Invalid after having set a valid one: No change:
    gc.setGenome(Arrays.asList(new String[] {"hg19"}), true);
    assertTrue(gc.getSamSeqDict() != null);
  }

  @Test
  public void canPrintRefSeq()
      throws InvalidGenomicCoordsException, IOException, InvalidColourException {
    GenomicCoords gc = new GenomicCoords("chr7:5566770-5566790", 80, samSeqDict, fastaFile);
    assertEquals("CACTTGGCCTCATTTTTAAGG\n", gc.printableRefSeq(true));
    // with format
    gc = new GenomicCoords("chr7:5566770-5566772", 80, samSeqDict, fastaFile);
    assertTrue(gc.printableRefSeq(false).contains("["));
  }

  @Test
  public void canConstructGenomicCoordsFromSpaceSepString()
      throws InvalidGenomicCoordsException, IOException {
    // Simple cases
    GenomicCoords gc = new GenomicCoords("chr7 1 100", 80, samSeqDict, fastaFile);
    assertEquals("chr7", gc.getChrom());
    assertEquals(1, (int) gc.getFrom());
    assertEquals(100, (int) gc.getTo());

    gc = new GenomicCoords("chr7 1,000 2,000", 80, samSeqDict, fastaFile);
    assertEquals("chr7", gc.getChrom());
    assertEquals(1000, (int) gc.getFrom());
    assertEquals(2000, (int) gc.getTo());

    gc = new GenomicCoords("chr7 11", 80, samSeqDict, fastaFile);
    assertEquals(11, (int) gc.getFrom());
    assertEquals(11 + 79, (int) gc.getTo());

    gc = new GenomicCoords("chr7 11", 80, samSeqDict, fastaFile);
    assertEquals(11, (int) gc.getFrom());
    assertEquals(11 + 79, (int) gc.getTo());

    gc = new GenomicCoords("chr7 11 1000 15000 2000", 80, samSeqDict, fastaFile);
    assertEquals(11, (int) gc.getFrom());
    assertEquals(2000, (int) gc.getTo());

    // Some mixing of space and - separator
    gc = new GenomicCoords("chr7 11-20", 80, samSeqDict, fastaFile);
    assertEquals(11, (int) gc.getFrom());
    assertEquals(20, (int) gc.getTo());

    gc = new GenomicCoords("chr7 11 -   20", 80, samSeqDict, fastaFile);
    assertEquals(11, (int) gc.getFrom());
    assertEquals(20, (int) gc.getTo());

    gc = new GenomicCoords("chr7 11- 20", 80, samSeqDict, fastaFile);
    assertEquals(11, (int) gc.getFrom());
    assertEquals(20, (int) gc.getTo());
  }

  @Test
  public void canConstructGenomicCoords() throws InvalidGenomicCoordsException, IOException {

    GenomicCoords gc = new GenomicCoords("chr7:1-100", 80, samSeqDict, fastaFile);
    assertEquals("chr7", gc.getChrom());
    assertEquals(1, (int) gc.getFrom());
    assertEquals(100, (int) gc.getTo());

    gc = new GenomicCoords("chr7:11", 80, samSeqDict, fastaFile);
    assertEquals(11 + 79, (int) gc.getTo());

    gc = new GenomicCoords("chr7", 80, samSeqDict, fastaFile);
    assertEquals(1, (int) gc.getFrom());
    assertEquals(80, (int) gc.getTo());

    gc =
        new GenomicCoords(
            "chr7:1000000000-1000000000", 80, samSeqDict, fastaFile); // Reset to size of chrom
    assertEquals(159138663, (int) gc.getFrom());

    gc = new GenomicCoords("chr7:100", 80, samSeqDict, fastaFile);
    assertEquals(100, (int) gc.getFrom());
    assertEquals(100 + 79, (int) gc.getTo());

    gc = new GenomicCoords("chr7:1,000-2,000", 80, samSeqDict, fastaFile);
    assertEquals(1000, (int) gc.getFrom());
    assertEquals(2000, (int) gc.getTo());

    gc = new GenomicCoords("chr7:1,000,000,000", 80, samSeqDict, null);
    assertEquals(159138663 - 79, (int) gc.getFrom()); // Reset to chrom size
    assertEquals(159138663, (int) gc.getTo());

    gc = new GenomicCoords("chr7:1,000,000,000", 80, null, null);
    assertEquals(1000000000, (int) gc.getFrom()); // Fine, no dict to check against.
    assertEquals(1000000000 + 79, (int) gc.getTo());

    gc = new GenomicCoords("chr99:1:1000-2000", 80, null, null);
    assertEquals("chr99:1", gc.getChrom());
    assertEquals(1000, (int) gc.getFrom());
    assertEquals(2000, (int) gc.getTo());

    gc = new GenomicCoords("chr99:1:1000", 80, null, null);
    assertEquals("chr99:1", gc.getChrom());
    assertEquals(1000, (int) gc.getFrom());
    assertEquals(1000 + 79, (int) gc.getTo());
  }

  @Test
  public void canGetRefSeq() throws InvalidGenomicCoordsException, IOException {
    GenomicCoords gc = new GenomicCoords("chr7:5566770-5566790", 80, samSeqDict, fastaFile);
    assertEquals("CACTTGGCCTCATTTTTAAGG", new String(gc.getRefSeq()));
    gc = new GenomicCoords("chr7:1-80", 79, samSeqDict, fastaFile);
    assertEquals(null, gc.getRefSeq());

    // Return seq even if len(seq) > windowSize
    // assertEquals("CACTTGGCCTCATTTTTAAGG", new String(gc.getRefSeq()));
  }

  @Test(expected = InvalidGenomicCoordsException.class)
  public void canThrowChromNotInDict() throws InvalidGenomicCoordsException, IOException {
    new GenomicCoords("nonsense:1-100", 80, samSeqDict, fastaFile);
  }

  @Test
  public void canZoom() throws InvalidGenomicCoordsException, IOException {
    GenomicCoords gc = new GenomicCoords("chr1:101-105", 80, samSeqDict, null);
    gc.zoomOut();
    assertEquals(99, (int) gc.getFrom()); // exp 95,111 if zoom fact is x2
    assertEquals(107, (int) gc.getTo());
    for (int i = 0; i < 40; i++) {
      gc.zoomOut();
    }
    assertEquals(1, (int) gc.getFrom());
    assertEquals(
        samSeqDict.getSequence("chr1").getSequenceLength(),
        (int) gc.getTo()); // Doesn't extend beyond chrom

    // Zoom out from a single base window at start of chrom
    gc = new GenomicCoords("chrM:1-1", 80, samSeqDict, null);
    gc.zoomOut();
    assertEquals(1, (int) gc.getFrom());
    assertTrue((int) gc.getTo() > 1);
    // 1 bp at the end of chrom
    gc = new GenomicCoords("chrM:16571-16571", 80, samSeqDict, null);
    gc.zoomOut();
    assertTrue((int) gc.getFrom() < 16571);
    assertEquals(16571, (int) gc.getTo());

    gc = new GenomicCoords("chr1:101-1000", 80, samSeqDict, null);
    gc.zoomIn();
    assertEquals(326, (int) gc.getFrom());
    assertEquals(776, (int) gc.getTo());

    // Zoom-in in small interval has no effect
    gc = new GenomicCoords("chr1:101-130", 80, samSeqDict, null);
    gc.zoomIn();
    assertEquals(101, (int) gc.getFrom());
    assertEquals(130, (int) gc.getTo());

    gc = new GenomicCoords("chr1:1-50", 80, samSeqDict, null);
    gc.zoomIn();
    assertEquals(1, (int) gc.getFrom());
    assertEquals(50, (int) gc.getTo());

    gc = new GenomicCoords("chrM:16561-16571", 80, samSeqDict, null); // End of chrom
    gc.zoomIn();
    assertEquals(16561, (int) gc.getFrom());
    assertEquals(16571, (int) gc.getTo());
  }

  @Test
  public void canPrepareRuler() throws InvalidGenomicCoordsException, IOException {
    GenomicCoords gc = new GenomicCoords("chr1:101-110", 80, samSeqDict, null);

    assertEquals(10, gc.getMapping().size());
    assertEquals(101.0, gc.getMapping().get(0), 0.01);
    assertEquals(102.0, gc.getMapping().get(1), 0.01);
    assertTrue(gc.isSingleBaseResolution());
  }

  @Test
  public void canPrintRuler()
      throws InvalidGenomicCoordsException,
          IOException,
          InvalidColourException,
          InvalidConfigException {

    GenomicCoords gc = new GenomicCoords("chr1:1-81", 80, samSeqDict, null);
    assertEquals(80, gc.printableGenomicRuler(10, true).length());
    assertTrue(gc.printableGenomicRuler(10, true).startsWith("1 "));

    // Can round labels
    gc = new GenomicCoords("chr1:2001234-3006789", 80, samSeqDict, null);
    String[] labels = gc.printableGenomicRuler(13, true).split(" +");
    for (int i = 0; i < labels.length - 1; i++) { // Do not test last label as it might be truncated
      assertTrue(labels[i].endsWith("000"));
    }

    assertTrue((gc.printableGenomicRuler(10, false).contains("[")));

    // Starts with 1 with large span
    gc = new GenomicCoords("chr1:1-10000", 80, null, null);
    assertTrue(gc.printableGenomicRuler(10, true).startsWith("1 "));
  }

  @Test
  public void canCenterAndExtendGenomicCoords() throws InvalidGenomicCoordsException, IOException {
    GenomicCoords gc = new GenomicCoords("chr1:10000-20000", 80, null, null);
    int size = 100;
    double slop = 5.0;
    gc.centerAndExtendGenomicCoords(gc, size, slop);
    assertEquals(9550, (int) gc.getFrom());
    assertEquals(10550, (int) gc.getTo());

    gc = new GenomicCoords("chr1:10000-20000", 80, null, null);
    size = 3;
    gc.centerAndExtendGenomicCoords(gc, size, 3.3); // Extended size smaller then windowSize
    assertTrue(gc.getUserWindowSize() <= gc.getGenomicWindowSize());

    gc = new GenomicCoords("chr1:1-300", 80, null, null);
    size = 100;
    gc.centerAndExtendGenomicCoords(gc, size, 5.0);
    assertEquals(1, (int) gc.getFrom());
  }

  @Test
  public void canPutFeatureInMidOfWindow() throws InvalidGenomicCoordsException, IOException {
    GenomicCoords gc = new GenomicCoords("chr1:100-1000", 80, null, null);
    System.err.println(Utils.getTerminalWidth());
    int size = 1; // not relevant
    double slop = 0;
    gc.centerAndExtendGenomicCoords(gc, size, slop);
    System.err.println(gc);
    assertEquals(60, (int) gc.getFrom());
    assertEquals(139, (int) gc.getTo());
  }
}
