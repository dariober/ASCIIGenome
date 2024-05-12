package tracks;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import coloring.Config;
import coloring.Xterm256;
import com.google.common.base.Splitter;
import exceptions.InvalidColourException;
import exceptions.InvalidCommandLineException;
import exceptions.InvalidConfigException;
import exceptions.InvalidGenomicCoordsException;
import exceptions.InvalidRecordException;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Pattern;
import org.junit.Before;
import org.junit.Test;
import samTextViewer.GenomicCoords;

public class TrackIntervalFeatureTest {

  @Before
  public void prepareConfig() throws IOException, InvalidConfigException {
    new Config(null);
    new Xterm256();
  }

  @Test
  public void canReadOddFilename()
      throws InvalidGenomicCoordsException,
          IOException,
          ClassNotFoundException,
          InvalidRecordException,
          SQLException {
    GenomicCoords gc = new GenomicCoords("1:1-100000", 80, null, null);
    TrackIntervalFeature tif = new TrackIntervalFeature("test_data/odd[filename].vcf.gz", gc);
    assertEquals("1", (tif.getGc().getChrom()));
    tif.close();
  }

  @Test
  public void canCloseFiles()
      throws InvalidGenomicCoordsException,
          IOException,
          ClassNotFoundException,
          InvalidRecordException,
          SQLException {
    GenomicCoords gc = new GenomicCoords("chr1:1-100000", 80, null, null);
    TrackIntervalFeature tif = new TrackIntervalFeature("test_data/hg19_genes_head.gtf", gc);
    tif.close();

    tif = new TrackIntervalFeature("test_data/wgEncodeDukeDnase8988T.fdr01peaks.hg19.bb", gc);
    tif.close();
  }

  @Test
  public void canColorGTFFeaturesByRegex()
      throws InvalidGenomicCoordsException,
          IOException,
          ClassNotFoundException,
          InvalidRecordException,
          SQLException,
          InvalidColourException {

    GenomicCoords gc = new GenomicCoords("chr1:1-100000", 80, null, null);
    TrackIntervalFeature tif = new TrackIntervalFeature("test_data/hg19_genes_head.gtf", gc);
    tif.printToScreen(); // This is to populate the ideograms.

    List<Argument> colorForRegex = new ArrayList<Argument>();
    colorForRegex.add(new Argument("DDX11L1", "216", false));
    tif.setColorForRegex(colorForRegex);
    assertTrue(tif.printToScreen().contains("216"));

    colorForRegex.clear();
    colorForRegex.add(new Argument("WASH7P", "233", false)); // 233:grey7 (almost black)
    tif.setColorForRegex(colorForRegex);
    assertTrue(tif.printToScreen().contains("233"));
    assertTrue(tif.printToScreen().contains("216"));
    assertTrue(tif.printToScreen().contains("253")); // Foreground color

    // Reset default
    tif.setColorForRegex(null);
    assertTrue(!tif.printToScreen().contains("216"));
  }

  @Test
  public void canColorFeaturesByAwk()
      throws InvalidGenomicCoordsException,
          IOException,
          ClassNotFoundException,
          InvalidRecordException,
          SQLException,
          InvalidColourException {

    GenomicCoords gc = new GenomicCoords("chr1:1-100000", 80, null, null);
    TrackIntervalFeature tif = new TrackIntervalFeature("test_data/hg19_genes_head.gtf", gc);
    tif.printToScreen(); // This is to populate the ideograms.

    List<Argument> colorForRegex = new ArrayList<Argument>();
    colorForRegex.add(new Argument("'$5 > 13000'", "216", false));
    tif.setColorForRegex(colorForRegex);
    assertTrue(tif.printToScreen().contains("216"));

    colorForRegex = new ArrayList<Argument>();
    colorForRegex.add(new Argument("$5 > 0", "100", false));
    tif.setColorForRegex(colorForRegex);
    assertTrue(!tif.printToScreen().contains("216"));
  }

  @Test
  public void canReturnFeaturesAsRawStrings()
      throws InvalidGenomicCoordsException,
          IOException,
          ClassNotFoundException,
          InvalidRecordException,
          SQLException {

    GenomicCoords gc = new GenomicCoords("chr1:1-100000", 80, null, null);
    TrackIntervalFeature tif = new TrackIntervalFeature("test_data/hg19_genes_head.gtf", gc);

    assertEquals(24, tif.getRecordsAsStrings().size());

    // Interval with no features
    gc = new GenomicCoords("FOO:1-100000", 80, null, null);
    tif = new TrackIntervalFeature("test_data/hg19_genes_head.gtf", gc);
    assertEquals(0, tif.getRecordsAsStrings().size());

    //        tif.setPrintRawLineCount(10);
    //        assertEquals(10, tif.getRecordsAsStrings().size());
    //
    //        tif.setPrintRawLineCount(0);
    //        assertEquals(0, tif.getRecordsAsStrings().size());
    //
    //        tif.setPrintRawLineCount(-1); // Return all
    //        assertEquals(24, tif.getRecordsAsStrings().size());

  }

  @Test
  public void canHandleGFFWithoutSupeFeatures()
      throws InvalidGenomicCoordsException,
          IOException,
          ClassNotFoundException,
          InvalidRecordException,
          SQLException,
          InvalidColourException,
          InvalidConfigException {
    // We have a GFF with only exons. Since there are no "transcripts", there is nothing to group
    // by.
    // See also issue #74.

    GenomicCoords gc = new GenomicCoords("chr1:11800-20000", 80, null, null);
    TrackIntervalFeature tif = new TrackIntervalFeature("test_data/issue74.gff3.gz", gc);
    tif.setNoFormat(true);
    assertEquals(10, tif.intervalFeatureList.size());
  }

  @Test
  public void canPrintChromsomeNames()
      throws ClassNotFoundException,
          IOException,
          InvalidGenomicCoordsException,
          InvalidRecordException,
          SQLException {

    String intervalFileName = "test_data/hg19_genes.gtf.gz";
    GenomicCoords gc = new GenomicCoords("7:5527151-5530709", 80, null, null);
    TrackIntervalFeature tif = new TrackIntervalFeature(intervalFileName, gc);

    assertTrue(tif.getChromosomeNames().size() > 10);

    tif = new TrackIntervalFeature("test_data/wgEncodeDukeDnase8988T.fdr01peaks.hg19.bb", gc);
    assertTrue(tif.getChromosomeNames().size() > 10);
  }

  @Test
  public void transcriptGTFToOneLine()
      throws InvalidGenomicCoordsException,
          IOException,
          ClassNotFoundException,
          InvalidRecordException,
          SQLException,
          InvalidColourException {

    // Expect:
    // ccccccccccc-----cccccae------------------------------------eeeeeeeee

    String intervalFileName = "test_data/hg19.gencode_genes_v19.gtf.gz";
    GenomicCoords gc = new GenomicCoords("chr7:5568562-5572120", 80, null, null);
    TrackIntervalFeature tif = new TrackIntervalFeature(intervalFileName, gc);
    tif.setNoFormat(true);
    System.err.println(tif.printToScreen());
    assertTrue(tif.printToScreen().trim().startsWith("ccc"));
    assertTrue(tif.printToScreen().trim().endsWith("eee"));
    assertEquals(6, tif.getIntervalFeatureList().size());
  }

  @Test
  public void canHideTitle()
      throws ClassNotFoundException,
          InvalidColourException,
          IOException,
          InvalidGenomicCoordsException,
          InvalidRecordException,
          SQLException {
    // See issue #42
    String intervalFileName = "test_data/hg19.gencode_genes_v19.gtf.gz";
    GenomicCoords gc = new GenomicCoords("chr7", 80, null, null);
    TrackIntervalFeature tif = new TrackIntervalFeature(intervalFileName, gc);
    tif.setHideTitle(true);
    assertEquals("", tif.getTitle());
  }

  @Test
  public void canPrintGFFRegionWithAndWithoutTranscripts()
      throws InvalidGenomicCoordsException,
          IOException,
          ClassNotFoundException,
          InvalidRecordException,
          SQLException,
          InvalidColourException {

    String intervalFileName = "test_data/Homo_sapiens.GRCh38.86.chromosome.7.gff3.gz";
    GenomicCoords gc = new GenomicCoords("7:1-100000", 80, null, null);
    TrackIntervalFeature tif = new TrackIntervalFeature(intervalFileName, gc);
    tif.setNoFormat(true);
    System.err.println(tif.printToScreen());
    assertTrue(tif.printToScreen().contains("||||||"));
    assertTrue(tif.printToScreen().contains("eee"));
  }

  @Test
  public void transcriptGFFToOneLine()
      throws InvalidGenomicCoordsException,
          IOException,
          ClassNotFoundException,
          InvalidRecordException,
          SQLException,
          InvalidColourException {

    String intervalFileName = "test_data/Homo_sapiens.GRCh38.86.ENST00000331789.gff3";
    GenomicCoords gc = new GenomicCoords("7:5527151-5530709", 80, null, null);
    TrackIntervalFeature tif = new TrackIntervalFeature(intervalFileName, gc);
    tif.setNoFormat(true);
    System.out.println("PRINTING:" + tif.printToScreen());
    assertTrue(tif.printToScreen().startsWith("uuuuu"));
    assertTrue(tif.printToScreen().endsWith("www"));

    tif.setNoFormat(false);
    assertTrue(tif.printToScreen().trim().startsWith("["));
  }

  @Test
  public void canAddNameToGFFTranscript()
      throws InvalidGenomicCoordsException,
          IOException,
          ClassNotFoundException,
          InvalidRecordException,
          SQLException,
          InvalidConfigException,
          InvalidColourException {

    String intervalFileName = "test_data/Homo_sapiens.GRCh38.86.ENST00000331789.gff3";

    GenomicCoords gc = new GenomicCoords("7:5527151-5530709", 80, null, null);
    TrackIntervalFeature tif = new TrackIntervalFeature(intervalFileName, gc);
    tif.setNoFormat(true);

    assertTrue(tif.printToScreen().contains("ACTB-001"));

    tif.setNoFormat(false);
    assertTrue(tif.printToScreen().trim().startsWith("["));
  }

  @Test
  public void canChangeFeatureName() throws Exception {

    // Complete GFF transscript
    String intervalFileName = "test_data/Homo_sapiens.GRCh38.86.ENST00000331789.gff3";

    GenomicCoords gc = new GenomicCoords("7:5527151-5530709", 80, null, null);
    TrackIntervalFeature tif = new TrackIntervalFeature(intervalFileName, gc);

    tif.setNoFormat(true);
    tif.setFeatureName("ID");
    assertTrue(tif.printToScreen().contains("transcript:ENS"));

    tif.setFeatureName("Name");
    assertTrue(tif.printToScreen().contains("ACTB-001"));
  }

  @Test
  public void canChangeFeatureName2() throws Exception {
    // GFF with features that are not transcript (e.g. "chromosome")
    String intervalFileName = "test_data/Homo_sapiens.GRCh38.86.chromosome.7.gff3.gz";
    GenomicCoords gc = new GenomicCoords("7:1-1000", 80, null, null);
    TrackIntervalFeature tif = new TrackIntervalFeature(intervalFileName, gc);
    tif.setNoFormat(true);
    assertTrue(tif.printToScreen().contains("chromosome:7"));

    tif.setFeatureName("Alias");
    System.err.println(tif.printToScreen());
    assertTrue(tif.printToScreen().contains("CM000669"));
  }

  @Test
  public void canAddNameToGTFTranscript()
      throws InvalidGenomicCoordsException,
          IOException,
          ClassNotFoundException,
          InvalidRecordException,
          SQLException,
          InvalidConfigException,
          InvalidColourException {

    String intervalFileName = "test_data/hg19_genes_head.gtf.gz";

    GenomicCoords gc = new GenomicCoords("chr1:11874-20000", 80, null, null);
    TrackIntervalFeature tif = new TrackIntervalFeature(intervalFileName, gc);
    tif.setNoFormat(true);
    assertTrue(tif.printToScreen().contains("NR_046018"));
  }

  @Test
  public void canSetColumnForBedFeature()
      throws InvalidGenomicCoordsException,
          IOException,
          ClassNotFoundException,
          InvalidRecordException,
          SQLException,
          InvalidColourException {
    String intervalFileName = "test_data/refSeq.hg19.short.sort.bed";

    GenomicCoords gc = new GenomicCoords("chr1:8404074-8404223", 80, null, null);
    TrackIntervalFeature tif = new TrackIntervalFeature(intervalFileName, gc);
    tif.setNoFormat(true);
    assertTrue(tif.printToScreen().contains("NM_001080397_utr3_8_0_chr1_8404074_f"));
    tif.setFeatureName("5"); // 1-based!
    assertTrue(tif.printToScreen().contains("_0_"));

    // Bigbed
    intervalFileName = "test_data/wgEncodeDukeDnase8988T.fdr01peaks.hg19.bb";

    gc = new GenomicCoords("chr1:1200328-1200477", 80, null, null);
    tif = new TrackIntervalFeature(intervalFileName, gc);
    tif.setNoFormat(true);
    assertTrue(tif.printToScreen().replaceAll("\\|", "").trim().isEmpty());

    tif.setFeatureName("5");
    assertTrue(tif.printToScreen().contains("_1000_"));

    // Name not to be shown
    tif.setFeatureName("-na");
    assertTrue(tif.printToScreen().replaceAll("\\|", "").trim().isEmpty());
  }

  @Test
  public void canReadBigBed()
      throws InvalidGenomicCoordsException,
          IOException,
          ClassNotFoundException,
          InvalidRecordException,
          SQLException {

    String filename = "test_data/wgEncodeDukeDnase8988T.fdr01peaks.hg19.bb";

    GenomicCoords gc = new GenomicCoords("chr1:1-800170", 80, null, null);
    TrackIntervalFeature tif = new TrackIntervalFeature(filename, gc);
    tif.setNoFormat(true);
    assertEquals(12, tif.getIntervalFeatureList().size());
    assertEquals(564665 + 1, tif.getIntervalFeatureList().get(0).getFrom());
  }

  @Test
  public void canReadTabix()
      throws IOException,
          InvalidGenomicCoordsException,
          ClassNotFoundException,
          InvalidRecordException,
          SQLException {

    String bgzFn =
        "test_data/refSeq.hg19.short.sort.bed.gz"; // "test_data/refSeq.hg19.short.sort.bed.gz";
    GenomicCoords gc = new GenomicCoords("chr1:16000000-20000000", 80, null, null);
    TrackIntervalFeature tif = new TrackIntervalFeature(bgzFn, gc);
    List<IntervalFeature> xset = tif.getFeaturesInInterval(gc.getChrom(), gc.getFrom(), gc.getTo());
    assertEquals(3, xset.size());
  }

  @Test
  public void canReadFromHTTP()
      throws IOException,
          InvalidGenomicCoordsException,
          ClassNotFoundException,
          InvalidRecordException,
          SQLException {
    GenomicCoords gc = new GenomicCoords("chr1:1-1000", 80, null, null);
    TrackIntervalFeature tif =
        new TrackIntervalFeature(
            "https://raw.githubusercontent.com/dariober/ASCIIGenome/master/test_data/refSeq.bed",
            gc);
    assertEquals("http", tif.getFilename().substring(0, 4));
    assertEquals(2, tif.getIntervalFeatureList().size());
  }

  @Test
  public void canReadTabixGTFFromHTTP()
      throws IOException,
          InvalidGenomicCoordsException,
          ClassNotFoundException,
          InvalidRecordException,
          SQLException {
    // If this file does not exist, put any valid tabix file and its index on Dropbox/Public and use
    // the dropbox link here.
    GenomicCoords gc = new GenomicCoords("chr1:64000-74000", 80, null, null);

    TrackIntervalFeature tif =
        new TrackIntervalFeature(
            "https://raw.githubusercontent.com/dariober/ASCIIGenome/master/test_data/hg19_genes_head.gtf.gz",
            gc);

    // We check the working file is on the remote server.
    assertEquals("http", tif.getFilename().substring(0, 4));
    assertEquals(
        "http",
        tif.getWorkFilename()
            .substring(
                0,
                4)); // Check we are using the remote file as working file. I.e. no need to download
                     // and index.
    assertEquals(4, tif.getIntervalFeatureList().size());
  }

  @Test
  public void canReadTabixVCFFromHTTP()
      throws IOException,
          InvalidGenomicCoordsException,
          ClassNotFoundException,
          InvalidRecordException,
          SQLException {
    String bgzFn =
        "https://raw.githubusercontent.com/dariober/ASCIIGenome/master/test_data/CHD.exon.2010_03.sites.vcf.gz";
    GenomicCoords gc = new GenomicCoords("1:1-2000000", 80, null, null);
    TrackIntervalFeature tif = new TrackIntervalFeature(bgzFn, gc);
    assertEquals(3, tif.getIntervalFeatureList().size());
    assertEquals("http", tif.getWorkFilename().substring(0, 4));
  }

  @Test
  public void canReadUnsortedVCFFromHTTP()
      throws IOException,
          InvalidGenomicCoordsException,
          ClassNotFoundException,
          InvalidRecordException,
          SQLException {
    GenomicCoords gc = new GenomicCoords("1:1-1142000", 80, null, null);
    TrackIntervalFeature tif =
        new TrackIntervalFeature(
            "https://raw.githubusercontent.com/dariober/ASCIIGenome/master/test_data/CHD.exon.2010_03.sites.unsorted.vcf",
            gc);
    assertEquals("http", tif.getFilename().substring(0, 4));
    assertEquals(3, tif.getIntervalFeatureList().size());
  }

  @Test
  public void canReadTabixVCFFromLocal()
      throws IOException,
          InvalidGenomicCoordsException,
          ClassNotFoundException,
          InvalidRecordException,
          SQLException {
    String bgzFn = "test_data/CHD.exon.2010_03.sites.vcf.gz";
    GenomicCoords gc = new GenomicCoords("1:1-2000000", 80, null, null);
    TrackIntervalFeature tif = new TrackIntervalFeature(bgzFn, gc);
    assertEquals(3, tif.getIntervalFeatureList().size());
  }

  @Test
  public void canReadBgzFileExtension()
      throws ClassNotFoundException,
          IOException,
          InvalidGenomicCoordsException,
          InvalidRecordException,
          SQLException {

    GenomicCoords gc = new GenomicCoords("1:1-200000000", 80, null, null);

    // .bgz, without index
    String intervalFileName = "test_data/bgz_noindex.vcf.bgz";
    TrackIntervalFeature tif = new TrackIntervalFeature(intervalFileName, gc);
    assertTrue(tif.getIntervalFeatureList().size() > 0);

    // .bgz, with index
    intervalFileName = "test_data/bgz_index.vcf.bgz";
    tif = new TrackIntervalFeature(intervalFileName, gc);
    assertTrue(tif.getFeaturesInInterval("1", 1, 200000000).size() > 0);
  }

  @Test
  public void canPrintGenotypeMatrix()
      throws InvalidGenomicCoordsException,
          IOException,
          ClassNotFoundException,
          InvalidRecordException,
          SQLException,
          InvalidColourException {

    GenomicCoords gc = new GenomicCoords("1:577583-759855", 80, null, null);
    String intervalFileName = "test_data/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz";
    TrackIntervalFeature tif = new TrackIntervalFeature(intervalFileName, gc);
    tif.setNoFormat(true);
    assertTrue(tif.printToScreen().contains("HG00096"));
  }

  @Test
  public void canPrintGtfFeatures()
      throws IOException,
          InvalidGenomicCoordsException,
          ClassNotFoundException,
          InvalidRecordException,
          SQLException,
          InvalidConfigException,
          InvalidColourException {

    String intervalFileName = "test_data/refSeq.bed";
    GenomicCoords gc = new GenomicCoords("chr1:1-70", 80, null, null);
    TrackIntervalFeature tif = new TrackIntervalFeature(intervalFileName, gc);
    tif.setNoFormat(true);
    assertEquals(2, tif.getIntervalFeatureList().size());

    assertEquals("||||", tif.printToScreen().substring(0, 4));
    tif.setNoFormat(false);
    assertTrue(
        tif.printToScreen().trim().startsWith("[")); // trim is necessary to remove escape \033
    assertTrue(tif.printToScreen().length() > 100);
  }

  @Test
  public void canConstructTrack()
      throws InvalidGenomicCoordsException,
          IOException,
          ClassNotFoundException,
          InvalidRecordException,
          SQLException,
          InvalidColourException {
    String intervalFileName = "test_data/refSeq.bed";
    GenomicCoords gc = new GenomicCoords("chr1:1-70", 80, null, null);
    TrackIntervalFeature tif = new TrackIntervalFeature(intervalFileName, gc);
    tif.setNoFormat(true);
    String exp = "||||_1_|||          ||||_2_|||                                        ";
    assertEquals(exp, tif.printToScreen());

    gc = new GenomicCoords("chr1:1-70", 80, null, null);
    tif = new TrackIntervalFeature(intervalFileName, gc);
    tif.setNoFormat(true);
    assertTrue(tif.printToScreen().startsWith("||||"));
    assertTrue(tif.printToScreen().endsWith("    "));
  }

  @Test
  public void canAssignFeatureText()
      throws InvalidGenomicCoordsException,
          IOException,
          ClassNotFoundException,
          InvalidRecordException,
          SQLException,
          InvalidColourException {
    String intervalFileName = "test_data/hg19_genes_head.gtf";
    GenomicCoords gc = new GenomicCoords("chr1:11874-12227", 80, null, null);
    TrackIntervalFeature tif = new TrackIntervalFeature(intervalFileName, gc);
    tif.setNoFormat(true);
    tif.setFeatureName("-na");
    assertTrue(tif.printToScreen().startsWith("EEEEE"));
  }

  @Test
  public void canStackFeatures()
      throws InvalidGenomicCoordsException,
          IOException,
          ClassNotFoundException,
          InvalidRecordException,
          SQLException,
          InvalidColourException {
    String intervalFileName = "test_data/overlapped.bed";
    GenomicCoords gc = new GenomicCoords("chr1:1-70", 80, null, null);
    TrackIntervalFeature tif = new TrackIntervalFeature(intervalFileName, gc);
    tif.setNoFormat(true);
    assertEquals(70 * 2 + 1, tif.printToScreen().length());

    String exp =
        ""
            + ">>>>>>>>>>     >>>>>>>>>>>>>>>          <<<<<<<<<<                    \n"
            + "     >>>>>>>>>>>>>>>                                                  ";
    assertEquals(exp, tif.printToScreen());
  }

  @Test
  public void canStackFeaturesInOneLine()
      throws InvalidGenomicCoordsException,
          IOException,
          ClassNotFoundException,
          InvalidRecordException,
          SQLException,
          InvalidColourException {
    String intervalFileName = "test_data/overlapped.bed";
    GenomicCoords gc = new GenomicCoords("chr1:1-70", 80, null, null);
    TrackIntervalFeature tif = new TrackIntervalFeature(intervalFileName, gc);
    tif.setNoFormat(true);
    tif.setFeatureDisplayMode(FeatureDisplayMode.ONELINE);

    String exp = ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>          <<<<<<<<<<";

    assertEquals(exp, tif.printToScreen().trim());
  }

  @Test
  public void canConstructTrackfromGtf()
      throws IOException,
          InvalidGenomicCoordsException,
          ClassNotFoundException,
          InvalidRecordException,
          SQLException {

    String intervalFileName = "test_data/hg19_genes.gtf.gz";
    GenomicCoords gc = new GenomicCoords("chr1:1-13000", 80, null, null);
    TrackIntervalFeature tif = new TrackIntervalFeature(intervalFileName, gc);

    gc = new GenomicCoords("chr7:5566000-5571000", 80, null, null);
    tif = new TrackIntervalFeature(intervalFileName, gc);
    System.out.println(tif);
  }

  @Test
  public void canShowHideFeatureByRegexFilters()
      throws InvalidGenomicCoordsException,
          IOException,
          ClassNotFoundException,
          InvalidRecordException,
          SQLException {

    String intervalFileName = "test_data/hg19_genes_head.gtf.gz";
    GenomicCoords gc = new GenomicCoords("chr1:10000-100000", 80, null, null);
    TrackIntervalFeature tif = new TrackIntervalFeature(intervalFileName, gc);
    tif.setNoFormat(true);

    tif.setShowHideRegex(
        Pattern.compile(Filter.DEFAULT_SHOW_REGEX.getValue()), Pattern.compile("\texon\t"));
    assertEquals(3, tif.getIntervalFeatureList().size());

    tif.setShowHideRegex(Pattern.compile("WASH7P"), Pattern.compile("^$"));
    assertTrue(tif.getIntervalFeatureList().size() == 11);
  }

  @Test
  public void canShowAndHide_getFeaturesInInterval()
      throws IOException,
          InvalidGenomicCoordsException,
          ClassNotFoundException,
          InvalidRecordException,
          SQLException {

    String intervalFileName = "test_data/hg19_genes_head.gtf.gz";
    GenomicCoords gc = new GenomicCoords("chr1:10000-100000", 80, null, null);
    TrackIntervalFeature tif = new TrackIntervalFeature(intervalFileName, gc);

    tif.setShowHideRegex(Pattern.compile("start_codon"), Pattern.compile("OR4F"));
    List<IntervalFeature> subset = tif.getFeaturesInInterval("chr1", 1, 500000000);
    assertEquals(40, subset.size());
  }

  @Test
  public void canStackFeaturesInBlocksToOneLine()
      throws InvalidGenomicCoordsException,
          IOException,
          ClassNotFoundException,
          InvalidRecordException,
          SQLException,
          InvalidColourException {
    // Test to reproduce issue#80
    String intervalFileName = "test_data/ovl.gff";
    GenomicCoords gc = new GenomicCoords("1:1-200", 80, null, null);
    TrackIntervalFeature tif = new TrackIntervalFeature(intervalFileName, gc);
    tif.setNoFormat(true);
    tif.setFeatureDisplayMode(FeatureDisplayMode.ONELINE);
    assertTrue(tif.printToScreen().trim().contains("|||   |||"));
  }

  @Test
  public void canApplyAwk_getFeaturesInInterval()
      throws IOException,
          InvalidGenomicCoordsException,
          ClassNotFoundException,
          InvalidRecordException,
          SQLException {

    String intervalFileName = "test_data/hg19_genes_head.gtf.gz";
    GenomicCoords gc = new GenomicCoords("chr1:10000-100000", 80, null, null);
    TrackIntervalFeature tif = new TrackIntervalFeature(intervalFileName, gc);

    String awk = "'$3 == \"start_codon\" && $9 !~ \"OR4F\"'";
    tif.setAwk(awk); // Note use single quotes

    assertEquals("-F '\\t' " + awk, tif.getAwk()); // Check -F arg has been prepended.

    List<IntervalFeature> subset = tif.getFeaturesInInterval("chr1", 1, 500000000);
    assertEquals(40, subset.size());

    tif.setAwk("-F \\t '($5 - $4) > 1000'"); // Filter for feature size > x
    subset = tif.getFeaturesInInterval("chr1", 1, 500000000);
    assertEquals(23, subset.size());

    tif.setAwk("  "); // Remove filter w/o args.
    subset = tif.getFeaturesInInterval("chr1", 1, 500000000);
    assertEquals(1000, subset.size());

    // Invalid script: Ugly stackTrace printed. All records returned
    boolean pass = false;
    try {
      tif.setAwk("foo_bar()");
    } catch (IOException e) {
      pass = true;
    }
    assertTrue(pass);
    assertEquals("", tif.getAwk()); // Faulty script has been removed.
    subset = tif.getFeaturesInInterval("chr1", 1, 500000000);
    assertEquals(1000, subset.size());
  }

  @Test
  public void canApplyAwkAndGrep_getFeaturesInInterval()
      throws IOException,
          InvalidGenomicCoordsException,
          ClassNotFoundException,
          InvalidRecordException,
          SQLException {

    String intervalFileName = "test_data/hg19_genes_head.gtf.gz";
    GenomicCoords gc = new GenomicCoords("chr1:10000-100000", 80, null, null);
    TrackIntervalFeature tif = new TrackIntervalFeature(intervalFileName, gc);

    tif.setAwk("'$3 == \"start_codon\"");
    tif.setShowHideRegex(
        Pattern.compile(Filter.DEFAULT_SHOW_REGEX.getValue()), Pattern.compile("OR4F"));
    List<IntervalFeature> subset = tif.getFeaturesInInterval("chr1", 1, 500000000);
    assertEquals(40, subset.size());
  }

  @Test
  public void canShowAndHide_coordsOfNextFeature()
      throws IOException,
          InvalidGenomicCoordsException,
          ClassNotFoundException,
          InvalidRecordException,
          SQLException {

    GenomicCoords gc = new GenomicCoords("chr1:1", 80, null, null);
    String intervalFileName = "test_data/hg19_genes_head.gtf.gz";
    TrackIntervalFeature tif = new TrackIntervalFeature(intervalFileName, gc);
    tif.setShowHideRegex(Pattern.compile(".*exon.*"), Pattern.compile(".*DDX11L1.*"));
    GenomicCoords curr = tif.coordsOfNextFeature(gc, false);
    assertEquals(14362, (int) curr.getFrom());
  }

  @Test
  public void canShowAndHide_findNextRegex()
      throws IOException,
          InvalidGenomicCoordsException,
          ClassNotFoundException,
          InvalidRecordException,
          SQLException {

    GenomicCoords gc = new GenomicCoords("chr1:1", 80, null, null);
    String intervalFileName = "test_data/hg19_genes_head.gtf.gz";
    TrackIntervalFeature tif = new TrackIntervalFeature(intervalFileName, gc);
    tif.setShowHideRegex(Pattern.compile(".*exon.*"), Pattern.compile(".*DDX11L1.*"));
    GenomicCoords curr = tif.findNextMatch(gc, Pattern.compile(".*gene_id.*"));
    assertEquals(14362, (int) curr.getFrom());
  }

  @Test
  public void canReadFileWithComments()
      throws IOException,
          InvalidGenomicCoordsException,
          ClassNotFoundException,
          InvalidRecordException,
          SQLException {

    GenomicCoords gc = new GenomicCoords("chr1:1-1000000", 80, null, null);
    String intervalFileName = "test_data/refSeq.bed";
    TrackIntervalFeature tif = new TrackIntervalFeature(intervalFileName, gc);
    assertEquals(2, tif.getFeaturesInInterval("chr1", 0, 100000).size());
  }

  @Test
  public void canFindNextFeatureOnChromGivenRegex()
      throws IOException,
          InvalidGenomicCoordsException,
          ClassNotFoundException,
          InvalidRecordException,
          SQLException {

    GenomicCoords gc = new GenomicCoords("chr1:1-1000000", 80, null, null);
    String intervalFileName = "test_data/refSeq.hg19.short.sort-2.bed";
    TrackIntervalFeature tif = new TrackIntervalFeature(intervalFileName, gc);

    IntervalFeature x = tif.findNextRegexInGenome(Pattern.compile(".*NM_.*"), "chr1", 20000000);
    assertTrue(x.getRaw().contains("NM_013943_utr3_5_0_chr1_25167429_f"));
    x = tif.findNextRegexInGenome(Pattern.compile(".*NM_.*"), "chr1", 80000000);
    assertTrue(x.getRaw().contains("NM_001080397_utr3_8_0_chr1_8404074_f"));

    x = tif.findNextRegexInGenome(Pattern.compile("NotPresent"), "chr1", 1);
    assertEquals(null, x);
  }

  @Test
  public void canFindIndel()
      throws IOException,
          InvalidGenomicCoordsException,
          ClassNotFoundException,
          InvalidRecordException,
          SQLException {

    GenomicCoords gc = new GenomicCoords("1:113050000", 80, null, null);
    String intervalFileName = "test_data/CEU.exon.2010_06.genotypes.vcf.gz";
    TrackIntervalFeature tif = new TrackIntervalFeature(intervalFileName, gc);

    IntervalFeature x = tif.findNextRegexInGenome(Pattern.compile(".*113054374.*"), "1", 113050000);
    assertTrue(x.getRaw().contains("\t113054374\t"));
    assertEquals(113054374, x.getFrom());
  }

  @Test
  /**
   * This should address issues #50 where a feature starting at the begining of the chrom is ignored
   */
  public void nextCanMoveToStartOfChrom()
      throws InvalidGenomicCoordsException,
          IOException,
          ClassNotFoundException,
          InvalidRecordException,
          SQLException {

    GenomicCoords gc = new GenomicCoords("chr1:2000-3000", 80, null, null);
    String intervalFileName = "test_data/refSeqZero.bed";
    TrackIntervalFeature tif = new TrackIntervalFeature(intervalFileName, gc);

    GenomicCoords newGc = tif.coordsOfNextFeature(gc, false);
    assertEquals(
        1,
        (int)
            newGc
                .getFrom()); // MEMO: Start of chrom is 0 in bed format but 1 in ASCIIGenome format

    // Backwards
    gc = new GenomicCoords("chr1:500-1000", 80, null, null);
    tif = new TrackIntervalFeature(intervalFileName, gc);
    newGc = tif.coordsOfNextFeature(gc, true);
    assertEquals(1, (int) newGc.getFrom());
  }

  @Test
  public void canGetCoordsOfNextFeature()
      throws IOException,
          InvalidGenomicCoordsException,
          ClassNotFoundException,
          InvalidRecordException,
          SQLException {

    GenomicCoords gc = new GenomicCoords("chr1:8000000-20000000", 80, null, null);
    TrackIntervalFeature tif =
        new TrackIntervalFeature("test_data/refSeq.hg19.short.sort-2.bed", gc);

    GenomicCoords newGc = tif.coordsOfNextFeature(gc, false);
    assertEquals(25167428 + 1, (int) newGc.getFrom());
    assertEquals(25167428 + gc.getGenomicWindowSize(), (int) newGc.getTo());

    // Next feature is on next chrom, current chrom not in file at all.
    gc = new GenomicCoords("foo:1-10000", 80, null, null);
    newGc = tif.coordsOfNextFeature(gc, false);
    assertEquals("chr1", newGc.getChrom());

    gc = new GenomicCoords("chr1:100000000-101000000", 80, null, null);
    newGc = tif.coordsOfNextFeature(gc, false);
    assertEquals("chr3", newGc.getChrom());

    gc = new GenomicCoords("chr1:10000000-10001000", 80, null, null);
    newGc = tif.coordsOfNextFeature(gc, true);
    assertEquals(8404074, (int) newGc.getFrom());

    newGc = tif.coordsOfNextFeature(newGc, true);
    assertEquals(67208779, (int) newGc.getFrom());
  }

  @Test
  public void canGetCoordsOfPreviousFeature()
      throws IOException,
          InvalidGenomicCoordsException,
          ClassNotFoundException,
          InvalidRecordException,
          SQLException {

    GenomicCoords gc = new GenomicCoords("chr1:8000000-20000000", 80, null, null);
    TrackIntervalFeature tif =
        new TrackIntervalFeature("test_data/refSeq.hg19.short.sort-2.bed", gc);

    gc = new GenomicCoords("chr1:10000000-10001000", 80, null, null);
    GenomicCoords newGc = tif.coordsOfNextFeature(gc, true);
    assertEquals(8404074, (int) newGc.getFrom());

    newGc = tif.coordsOfNextFeature(newGc, true);
    assertEquals(67208779, (int) newGc.getFrom());

    // Exactly at the start of a feature and move to previous one. This is the feature:
    // chrM hg19_wgEncodeGencodeBasicV19 exon 15957 16024 ...
    gc = new GenomicCoords("chrM:15957-17259", 80, null, null);
    tif = new TrackIntervalFeature("test_data/hg19.gencode_genes_v19.gtf.gz", gc);
    newGc = tif.coordsOfNextFeature(gc, true);
    assertEquals(
        15889, (int) newGc.getFrom()); // chrM hg19_wgEncodeGencodeBasicV19 exon 15889 15954

    gc = new GenomicCoords("chrM:14672-14898", 80, null, null);
    tif = new TrackIntervalFeature("test_data/hg19.gencode_genes_v19.gtf.gz", gc);
    newGc = tif.coordsOfNextFeature(gc, true);
    assertEquals(14150, (int) newGc.getFrom());
  }

  @Test
  public void canFindAllRegex()
      throws IOException,
          InvalidGenomicCoordsException,
          ClassNotFoundException,
          InvalidRecordException,
          SQLException {

    GenomicCoords gc = new GenomicCoords("chr18:1-10000", 80, null, null);
    TrackIntervalFeature tif = new TrackIntervalFeature("test_data/hg19_genes_head.gtf", gc);

    GenomicCoords matched =
        tif.genomicCoordsAllChromMatchInGenome(Pattern.compile(".*\"WASH7P\".*"), gc);
    assertEquals("chr1", matched.getChrom());
    assertEquals(14362, (int) matched.getFrom());
    assertEquals(29370, (int) matched.getTo());

    // No match
    gc = new GenomicCoords("chr18:1-10000", 80, null, null);
    tif = new TrackIntervalFeature("test_data/hg19_genes_head.gtf", gc);
    matched = tif.genomicCoordsAllChromMatchInGenome(Pattern.compile(".*\"FOOBAR\".*"), gc);
    assertEquals("chr18", matched.getChrom());
    assertEquals(1, (int) matched.getFrom());
    assertEquals(10000, (int) matched.getTo());
  }

  @Test
  public void canFetchInterval()
      throws IOException,
          InvalidGenomicCoordsException,
          ClassNotFoundException,
          InvalidRecordException,
          SQLException {

    GenomicCoords gc = new GenomicCoords("chr18:1-10000", 80, null, null);
    TrackIntervalFeature tif = new TrackIntervalFeature("test_data/refSeq.hg19.short.bed", gc);

    List<IntervalFeature> interval = tif.getFeaturesInInterval("chr1", 20000000, 40000000);
    assertEquals(25167428 + 1, interval.get(0).getFrom()); // Note adding 1 because bed is 0-based
    assertEquals(33586132, interval.get(interval.size() - 1).getTo());

    // Nothing to fetch: Range not in bed
    tif = new TrackIntervalFeature("test_data/refSeq.hg19.short.bed", gc);
    interval = tif.getFeaturesInInterval("chr1", 500000000, 600000000);
    assertEquals(0, interval.size());
    // Nothing to fetch: chrom not in bed:
    interval = tif.getFeaturesInInterval("chrNonSense", 1, 10);
    assertEquals(0, interval.size());
  }

  @Test
  public void canPrintRawLines()
      throws ClassNotFoundException,
          IOException,
          InvalidGenomicCoordsException,
          InvalidRecordException,
          SQLException,
          InvalidColourException,
          InvalidConfigException,
          InvalidCommandLineException {

    GenomicCoords gc = new GenomicCoords("chr1:1-40000", 80, null, null);
    TrackIntervalFeature tif = new TrackIntervalFeature("test_data/hg19_genes_head.gtf", gc);
    tif.setPrintMode(PrintRawLine.CLIP);

    tif.setPrintRawLineCount(-1);
    assertEquals(20, tif.printLines().split("\n").length);

    tif.setPrintRawLineCount(5);
    assertEquals(5 + 1, tif.printLines().split("\n").length); // +1 for the string of omitted count.
  }

  @Test
  public void canPrintNormalizedVcfLines()
      throws ClassNotFoundException,
          IOException,
          InvalidGenomicCoordsException,
          InvalidRecordException,
          SQLException,
          InvalidColourException,
          InvalidConfigException,
          InvalidCommandLineException {

    GenomicCoords gc = new GenomicCoords("1:645709-645975", 80, null, null);
    TrackIntervalFeature tif =
        new TrackIntervalFeature("test_data/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz", gc);
    tif.setPrintMode(PrintRawLine.FULL);
    tif.setNoFormat(true);
    tif.setPrintNormalizedVcf(true);

    String out = tif.printLines();
    assertEquals(3, out.split("\n").length);
    assertTrue(out.contains(" HG00096 | GT"));

    // VCF with without samples
    gc = new GenomicCoords("1:1105467-1105647", 80, null, null);
    tif = new TrackIntervalFeature("test_data/CHD.exon.2010_03.sites.vcf.gz", gc);
    tif.setPrintMode(PrintRawLine.FULL);
    tif.setNoFormat(true);
    tif.setPrintNormalizedVcf(true);

    assertEquals(1, tif.printLines().split("\n").length);
  }

  @Test
  public void canPrintFormattedVepAnnotation()
      throws InvalidGenomicCoordsException,
          IOException,
          ClassNotFoundException,
          InvalidRecordException,
          SQLException,
          InvalidColourException,
          InvalidCommandLineException {
    GenomicCoords gc = new GenomicCoords("chr1:14327-14836", 80, null, null);
    TrackIntervalFeature tif = new TrackIntervalFeature("test_data/vep.vcf", gc);
    tif.setPrintMode(PrintRawLine.FULL);
    tif.setNoFormat(true);
    String woVep = tif.printLines();
    tif.setPrintFormattedVep("");
    String printed = tif.printLines();
    assertTrue(Splitter.on("\n").splitToList(tif.printLines()).size() > 50);
    assertTrue(printed.contains("Consequence "));
    assertTrue(!printed.contains("SWISSPROT"));

    // INFO tag not found: Do nothing
    tif.setPrintFormattedVep("csq_na");
    printed = tif.printLines();
    assertEquals(woVep, printed);

    // Only ask for some headers, case insensitive
    tif.setPrintFormattedVep("CSQ,ConseQUENCE,Allele");
    printed = tif.printLines();
    assertTrue(printed.contains("Consequence"));
    assertTrue(printed.contains("Allele"));
    assertTrue(!printed.contains("IMPACT"));

    // Omit CSQ
    tif.setPrintFormattedVep("CSQ,null");
    printed = tif.printLines();
    assertTrue(printed.contains("CSQ=... "));

    // No effect without VEP tag
    gc = new GenomicCoords("1:1105468-34435998", 80, null, null);
    tif = new TrackIntervalFeature("test_data/CHD.exon.2010_03.sites.vcf", gc);
    tif.setPrintMode(PrintRawLine.FULL);
    tif.setNoFormat(true);
    tif.setPrintFormattedVep(null);
    woVep = tif.printLines();
    tif.setPrintFormattedVep("");
    String withVep = tif.printLines();
    assertEquals(woVep, withVep);

    // No effect on non-VCF track
    gc = new GenomicCoords("chr18:1-10000", 80, null, null);
    tif = new TrackIntervalFeature("test_data/refSeq.hg19.short.bed", gc);
    tif.setPrintMode(PrintRawLine.FULL);
    tif.setNoFormat(true);
    tif.setPrintFormattedVep(null);
    woVep = tif.printLines();
    tif.setPrintFormattedVep("");
    withVep = tif.printLines();
    assertEquals(woVep, withVep);
  }

  @Test
  public void canPrintMappingOfFeaturesToScreen()
      throws IOException,
          InvalidGenomicCoordsException,
          ClassNotFoundException,
          InvalidRecordException,
          SQLException {
    List<Double> rulerMap = new ArrayList<Double>();
    for (int i = 14000; i < 14400; i += 10) {
      rulerMap.add((double) i);
    }
    System.out.println(rulerMap);

    GenomicCoords gc = new GenomicCoords("chr18:1-10000", 80, null, null);
    TrackIntervalFeature tif = new TrackIntervalFeature("test_data/refSeq.hg19.short.bed", gc);
    System.out.println(tif.getFeaturesInInterval("chr1", 0, 1000000000).get(0));
  }

  @Test
  public void canProcessIndelAtWindowBoundary()
      throws IOException,
          InvalidGenomicCoordsException,
          ClassNotFoundException,
          InvalidRecordException,
          SQLException,
          InvalidColourException {

    // VCF record:
    // 1 113054374 . CTTG C
    // The genomic start is at 113054374, which is inside the GenomicCoords interval.
    // However, the printable coordinates start at 113054374 + 1 because the first
    // base is equal to reference. This means that at the boundary we have a feature
    // that is in this interval but not visible.

    GenomicCoords gc = new GenomicCoords("1:113054305-113054375", 70, null, null);
    TrackIntervalFeature tif =
        new TrackIntervalFeature("test_data/CEU.exon.2010_06.genotypes.vcf.gz", gc);
    tif.setNoFormat(true);

    // VCF record is in interval and visible:
    assertEquals(1, tif.getIntervalFeatureList().size());
    System.err.println(tif.printToScreen());
    assertTrue(tif.printToScreen().contains("D")); // Deletion

    // Now the end coordinate equals the start of the deletion. Since the first base of the
    // deletion is equal to ref, there is nothing to print on screen:
    gc = new GenomicCoords("1:113054305-113054374", 70, null, null);
    tif = new TrackIntervalFeature("test_data/CEU.exon.2010_06.genotypes.vcf.gz", gc);
    tif.setNoFormat(true);
    // tif.getGenotypeMatrix().setnMaxSamples(0);

    assertEquals(1, tif.getIntervalFeatureList().size()); // Feature is in interval
    assertTrue(tif.printToScreen().length() > 70); // But nothing to print
    assertTrue(!tif.printToScreen().contains("D")); // No deletion visible
  }

  @Test
  public void canReadVCFTabix()
      throws IOException,
          InvalidGenomicCoordsException,
          ClassNotFoundException,
          InvalidRecordException,
          SQLException,
          InvalidColourException {

    GenomicCoords gc = new GenomicCoords("chr18:1-10000", 80, null, null);
    TrackIntervalFeature tif =
        new TrackIntervalFeature("test_data/CHD.exon.2010_03.sites.vcf.gz", gc);

    List<IntervalFeature> xset = tif.getFeaturesInInterval("1", 1, 10000000);
    assertEquals(9, xset.size());
    IntervalFeature x = xset.get(1);
    assertEquals("1", x.getChrom());
    assertEquals(1108138, x.getFrom());
    System.err.println(tif.printToScreen());
  }

  @Test
  public void canReadUnsortedVCF()
      throws IOException,
          InvalidGenomicCoordsException,
          ClassNotFoundException,
          InvalidRecordException,
          SQLException {

    GenomicCoords gc = new GenomicCoords("chr18:1-10000", 80, null, null);
    TrackIntervalFeature tif =
        new TrackIntervalFeature("test_data/CHD.exon.2010_03.sites.unsorted.vcf", gc);
    List<IntervalFeature> xset = tif.getFeaturesInInterval("1", 1, 10000000);
    assertEquals(9, xset.size());
    IntervalFeature x = xset.get(1);
    assertEquals("1", x.getChrom());
    assertEquals(1108138, x.getFrom());
  }

  @Test
  public void canReadFeaturesOfLengthOne()
      throws IOException,
          InvalidGenomicCoordsException,
          ClassNotFoundException,
          InvalidRecordException,
          SQLException {

    GenomicCoords gc = new GenomicCoords("chr18:1-10000", 80, null, null);
    TrackIntervalFeature tif = new TrackIntervalFeature("test_data/refSeqZero.bed", gc);

    List<IntervalFeature> xset = tif.getFeaturesInInterval("chr1", 1, 100);
    assertEquals(1, xset.size());
  }

  // @Test
  // public void canReadFromURL() throws IOException, InvalidGenomicCoordsException,
  // ClassNotFoundException, InvalidRecordException, SQLException{

  //     System.err.println("canReadFromURL: This can take  a while...");
  //     String urlStr=
  // "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeHaibTfbs/wgEncodeHaibTfbsA549Atf3V0422111Etoh02PkRep1.broadPeak.gz";

  //     GenomicCoords gc= new GenomicCoords("chr18:1-10000", 80, null, null);
  //     TrackIntervalFeature tif= new TrackIntervalFeature(urlStr, gc);
  //
  //     List<IntervalFeature> xset = tif.getFeaturesInInterval("chr1", 1, 1000000);
  //     assertEquals(2, xset.size());
  //     assertEquals(878407+1, xset.get(0).getFrom());
  //
  // }

}
