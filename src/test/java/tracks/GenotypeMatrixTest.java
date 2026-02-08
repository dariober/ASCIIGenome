package tracks;

import static org.junit.Assert.*;

import colouring.Config;
import colouring.Xterm256;
import exceptions.InvalidColourException;
import exceptions.InvalidCommandLineException;
import exceptions.InvalidConfigException;
import exceptions.InvalidGenomicCoordsException;
import exceptions.InvalidRecordException;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.List;
import org.junit.Before;
import org.junit.Test;
import samTextViewer.GenomicCoords;

public class GenotypeMatrixTest {

  @Before
  public void config() throws IOException, InvalidConfigException {
    new Config(null);
    new Xterm256();
  }

  @Test
  public void bigData()
      throws InvalidGenomicCoordsException,
          IOException,
          ClassNotFoundException,
          InvalidRecordException,
          SQLException,
          InvalidColourException,
          InvalidCommandLineException {

    VCFFileReader reader =
        new VCFFileReader(new File("test_data/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz"));
    VCFHeader vcfHeader = reader.getFileHeader();
    reader.close();
    GenomicCoords gc = new GenomicCoords("1:5934301-6000000", 80, null, null);
    TrackVCF vcf = new TrackVCF("test_data/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz", gc);

    List<VCFFeature> linf = vcf.getIntervalFeatureList();

    GenotypeMatrix gm = new GenotypeMatrix();
    gm.setPyScriptFilter("{GT} == 10 or 1>2");
    gm.printToScreen(true, linf, 80, vcfHeader);
  }

  @Test
  public void canDetectAmbigousTags()
      throws InvalidGenomicCoordsException,
          IOException,
          ClassNotFoundException,
          InvalidRecordException,
          SQLException,
          InvalidColourException,
          IOException {

    VCFFileReader reader = new VCFFileReader(new File("test_data/info_formats.vcf.gz"));
    VCFHeader vcfHeader = reader.getFileHeader();
    reader.close();
    GenomicCoords gc = new GenomicCoords("1:17822074-17822184", 80, null, null);
    TrackVCF vcf = new TrackVCF("test_data/info_formats.vcf.gz", gc);
    List<VCFFeature> linf = vcf.getIntervalFeatureList();
    GenotypeMatrix gm = new GenotypeMatrix();
    // No filter
    assertEquals(2, gm.printToScreen(true, linf, 80, vcfHeader).split("\n").length);

    gm.setPyScriptFilter("{XA}[1] > 0");
    boolean pass = false;
    try {
      gm.printToScreen(true, linf, 80, vcfHeader);
    } catch (Exception e) {
      if (e.getMessage().contains("XA")) {
        pass = true;
      }
      ;
    }
    assertTrue(pass);
  }

  @Test
  public void canFilterWithPython()
      throws InvalidGenomicCoordsException,
          ClassNotFoundException,
          InvalidRecordException,
          SQLException,
          InvalidColourException,
          IOException {

    VCFFileReader reader = new VCFFileReader(new File("test_data/info_formats.vcf.gz"));
    VCFHeader vcfHeader = reader.getFileHeader();
    reader.close();

    GenomicCoords gc = new GenomicCoords("1:17822074-17822184", 80, null, null);
    TrackVCF vcf = new TrackVCF("test_data/info_formats.vcf.gz", gc);

    List<VCFFeature> linf = vcf.getIntervalFeatureList();

    GenotypeMatrix gm = new GenotypeMatrix();

    // No filter
    assertEquals(2, gm.printToScreen(true, linf, 80, vcfHeader).split("\n").length);

    // Filter by array in FORMAT
    gm.setPyScriptFilter("{FMT/XA}[1] > 0");
    String y = gm.printToScreen(true, linf, 80, vcfHeader);
    assertTrue(y.contains("sample1"));
    assertTrue(!y.contains("sample2"));

    // One sample satisfies DP at at least one record.
    gm.setPyScriptFilter("{DP} >= 100");
    String x = gm.printToScreen(true, linf, 80, vcfHeader);
    assertEquals(1, x.split("\n").length);

    // Filter by array in INFO
    gm.setPyScriptFilter("{INFO/XA}[0] > 0");
    x = gm.printToScreen(true, linf, 80, vcfHeader);
    assertTrue(x.contains("sample1"));

    // Exclude all samples.
    gm.setPyScriptFilter("{DP} >= 10000");
    x = gm.printToScreen(true, linf, 80, vcfHeader);
    assertTrue(x.isEmpty());

    // Exclude all sample II: Non sense comparison
    gm.setPyScriptFilter("{POS} == 'G' and {DP} > 5");
    x = gm.printToScreen(true, linf, 80, vcfHeader);
    assertTrue(x.isEmpty());

    // Remove filter
    gm.setPyScriptFilter("1 > 2"); // First remove all samples
    x = gm.printToScreen(true, linf, 80, vcfHeader);
    assertTrue(x.isEmpty());

    gm.setPyScriptFilter(null); // Then remove filter and check samples are back
    x = gm.printToScreen(true, linf, 80, vcfHeader);
    assertTrue(!x.isEmpty());

    gm.setPyScriptFilter(""); // Empty also works to remove filter
    x = gm.printToScreen(true, linf, 80, vcfHeader);
    assertTrue(!x.isEmpty());

    // Handle missing values
    gm.setPyScriptFilter("{ID} == '.'"); // Use literal '.', not the null value.
    x = gm.printToScreen(true, linf, 80, vcfHeader);
    assertTrue(x.trim().length() > 50);

    // Filter by ALT: Note array slicing
    // NB: Although all samples satisfy ALT, only sample2 also satisfies DP
    gm.setPyScriptFilter("{ALT}[0] == 'G' and {DP} > 5");
    x = gm.printToScreen(true, linf, 80, vcfHeader);
    assertEquals(1, x.split("\n").length);
    assertTrue(x.contains("sample2"));

    // Same using POS
    gm.setPyScriptFilter("{POS} == 17822094 and {DP} > 5");
    x = gm.printToScreen(true, linf, 80, vcfHeader);
    assertEquals(1, x.split("\n").length);
    assertTrue(x.contains("sample2"));
  }

  @Test
  public void flagFilterInPy()
      throws InvalidGenomicCoordsException,
          IOException,
          SQLException,
          ClassNotFoundException,
          InvalidRecordException,
          InvalidColourException {
    VCFFileReader reader = new VCFFileReader(new File("test_data/info_formats.vcf.gz"));
    VCFHeader vcfHeader = reader.getFileHeader();
    reader.close();

    GenomicCoords gc = new GenomicCoords("1:17822074-17822184", 80, null, null);
    TrackVCF vcf = new TrackVCF("test_data/info_formats.vcf.gz", gc);

    List<VCFFeature> linf = vcf.getIntervalFeatureList();

    GenotypeMatrix gm = new GenotypeMatrix();

    // Filter by FLAG field type
    gm.setPyScriptFilter("{XB} and {DP} == 100");
    String x = gm.printToScreen(true, linf, 80, vcfHeader);
    assertTrue(x.contains("sample2"));
    assertTrue(!x.contains("sample1"));
  }

  @Test
  public void genotypeAsNumericInPy()
      throws ClassNotFoundException,
          IOException,
          InvalidGenomicCoordsException,
          InvalidRecordException,
          SQLException,
          InvalidColourException,
          IOException {
    VCFFileReader reader = new VCFFileReader(new File("test_data/info_formats.vcf.gz"));
    VCFHeader vcfHeader = reader.getFileHeader();
    reader.close();
    GenomicCoords gc = new GenomicCoords("1:17822074-17822184", 80, null, null);
    TrackVCF vcf = new TrackVCF("test_data/info_formats.vcf.gz", gc);
    List<VCFFeature> linf = vcf.getIntervalFeatureList();

    GenotypeMatrix gm = new GenotypeMatrix();

    gm.setPyScriptFilter("{GT} == '1/0'");
    String x = gm.printToScreen(true, linf, 80, vcfHeader);
    assertTrue(x.contains("sample1"));
    assertFalse(x.contains("sample2"));

    gm.setPyScriptFilter("{GT} == './.'");
    x = gm.printToScreen(true, linf, 80, vcfHeader);
    assertTrue(x.contains("sample1"));
    assertFalse(x.contains("sample2"));
  }

  @Test
  public void genotypeKeywordInPy()
      throws ClassNotFoundException,
          IOException,
          InvalidGenomicCoordsException,
          InvalidRecordException,
          SQLException,
          InvalidColourException,
          IOException {

    VCFFileReader reader = new VCFFileReader(new File("test_data/info_formats.vcf.gz"));
    VCFHeader vcfHeader = reader.getFileHeader();
    reader.close();
    GenomicCoords gc = new GenomicCoords("1:17822074-17822184", 80, null, null);
    TrackVCF vcf = new TrackVCF("test_data/info_formats.vcf.gz", gc);
    List<VCFFeature> linf = vcf.getIntervalFeatureList();

    GenotypeMatrix gm = new GenotypeMatrix();

    gm.setPyScriptFilter("{HOM}");
    String x = gm.printToScreen(true, linf, 80, vcfHeader);
    assertTrue(x.contains("sample2"));
    assertFalse(x.contains("sample1"));

    gm.setPyScriptFilter("{NO_CALL}");
    x = gm.printToScreen(true, linf, 80, vcfHeader);
    assertTrue(x.contains("sample1"));
    assertFalse(x.contains("sample2"));

    gm.setPyScriptFilter("{HET_NON_REF}");
    x = gm.printToScreen(true, linf, 80, vcfHeader);
    assertTrue(x.contains("sample2"));
    assertFalse(x.contains("sample1"));

    gm.setPyScriptFilter("{CALLED} and {POS} == 17822094");
    x = gm.printToScreen(true, linf, 80, vcfHeader);
    assertTrue(x.contains("sample2"));
    assertFalse(x.contains("sample1"));
  }

  @Test
  public void invalidPyScripts()
      throws ClassNotFoundException,
          IOException,
          InvalidGenomicCoordsException,
          InvalidRecordException,
          SQLException,
          InvalidColourException,
          IOException,
          InvalidCommandLineException {

    VCFFileReader reader = new VCFFileReader(new File("test_data/info_formats.vcf.gz"));
    VCFHeader vcfHeader = reader.getFileHeader();
    reader.close();
    GenomicCoords gc = new GenomicCoords("1:17822074-17822184", 80, null, null);
    TrackVCF vcf = new TrackVCF("test_data/info_formats.vcf.gz", gc);
    List<VCFFeature> linf = vcf.getIntervalFeatureList();

    GenotypeMatrix gm = new GenotypeMatrix();

    // Script must return boolean
    gm.setPyScriptFilter("{DP} > 5 and 10 + 3");
    gm.printToScreen(true, linf, 80, vcfHeader);
    assertNull(gm.getPyScriptFilter()); // After faulty script reset to null.

    // Invalid syntax. E.g. when {TAG} does not exist.
    gm.setPyScriptFilter("{FOOBAR} > 10");
    gm.printToScreen(true, linf, 80, vcfHeader);
    assertNull(gm.getPyScriptFilter());

    // After exception the invalid filter has been removed
    gm.printToScreen(true, linf, 80, vcfHeader);
  }

  @Test
  public void canFilterGenotypeWithPy()
      throws IOException,
          ClassNotFoundException,
          IOException,
          InvalidGenomicCoordsException,
          InvalidRecordException,
          SQLException,
          InvalidColourException,
          InvalidCommandLineException {

    VCFFileReader reader = new VCFFileReader(new File("test_data/info_formats.vcf.gz"));
    VCFHeader vcfHeader = reader.getFileHeader();
    reader.close();
    GenomicCoords gc = new GenomicCoords("1:17822074-17822184", 80, null, null);
    TrackVCF vcf = new TrackVCF("test_data/info_formats.vcf.gz", gc);
    List<VCFFeature> linf = vcf.getIntervalFeatureList();
    GenotypeMatrix gm = new GenotypeMatrix();
    // Missing alleles
    gm.setPyScriptFilter("{GT} == \"./.\"");
    String x = gm.printToScreen(true, linf, 80, vcfHeader);
    assertTrue(x.contains("sample1"));
    assertFalse(x.contains("sample2"));

    // Using FMT
    // Missing alleles
    gm.setPyScriptFilter("{FMT/GT} == \"./.\"");
    x = gm.printToScreen(true, linf, 80, vcfHeader);
    assertTrue(x.contains("sample1"));
    assertFalse(x.contains("sample2"));

    gm = new GenotypeMatrix();
    gm.setPyScriptFilter("{FMT/XA}[0] > 0");
    x = gm.printToScreen(true, linf, 80, vcfHeader);
    assertTrue(x.contains("sample1"));
    assertFalse(x.contains("sample2"));

    gm = new GenotypeMatrix();
    gm.setPyScriptFilter("{POS} == 17822092 and {FMT/XA}[0] == 0");
    x = gm.printToScreen(true, linf, 80, vcfHeader);
    assertTrue(x.contains("sample2"));
    assertFalse(x.contains("sample1"));

    gm = new GenotypeMatrix();
    gm.setPyScriptFilter("{POS} == 17822092 and {FMT/XA}[1] == 2");
    x = gm.printToScreen(true, linf, 80, vcfHeader);
    assertTrue(x.contains("sample1"));
    assertFalse(x.contains("sample2"));

    gm = new GenotypeMatrix();
    gm.setPyScriptFilter("{POS} == 17822092 and {FMT/XA}[1] == 3");
    x = gm.printToScreen(true, linf, 80, vcfHeader);
    assertFalse(x.contains("sample1"));
    assertFalse(x.contains("sample2"));
  }

  @Test
  public void canInitMatrix()
      throws IOException,
          InvalidGenomicCoordsException,
          IOException,
          ClassNotFoundException,
          InvalidRecordException,
          SQLException,
          InvalidColourException {

    VCFFileReader reader =
        new VCFFileReader(new File("test_data/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz"));
    VCFHeader vcfHeader = reader.getFileHeader();
    reader.close();

    GenomicCoords gc = new GenomicCoords("1:572807-755079", 80, null, null);
    TrackVCF vcf = new TrackVCF("test_data/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz", gc);

    List<VCFFeature> linf = vcf.getIntervalFeatureList();

    GenotypeMatrix gm = new GenotypeMatrix();

    // No feature at all
    assertTrue(gm.printToScreen(true, new ArrayList<VCFFeature>(), 80, vcfHeader).isEmpty());

    String x = gm.printToScreen(true, linf, 80, vcfHeader);

    // Check sample name
    assertTrue(x.startsWith("HG00096"));

    String[] rows = x.split("\n");
    assertEquals(3, rows.length);
    assertEquals(80, rows[0].length());

    // Check genotype coding
    assertTrue(rows[0].contains("?")); // Missing genotype
    assertTrue(rows[0].contains(".")); // HOM ref
    assertTrue(rows[1].contains("E")); // HET
    assertTrue(rows[2].contains("0")); // HOM alt
  }

  @Test
  public void overalappingSymbols()
      throws IOException,
          InvalidGenomicCoordsException,
          IOException,
          ClassNotFoundException,
          InvalidRecordException,
          SQLException,
          InvalidColourException {

    VCFFileReader reader =
        new VCFFileReader(new File("test_data/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz"));
    VCFHeader vcfHeader = reader.getFileHeader();
    reader.close();

    // Same genotype: no ambiguity.
    GenomicCoords gc = new GenomicCoords("1:199882-200100", 80, null, null);
    TrackVCF vcf = new TrackVCF("test_data/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz", gc);
    List<VCFFeature> linf = vcf.getIntervalFeatureList();

    GenotypeMatrix gm = new GenotypeMatrix();
    String x = gm.printToScreen(false, linf, 80, vcfHeader);
    assertTrue(x.contains("O")); // Two genotype at the same position: Both the same
    assertTrue(x.contains("*")); // Different
  }

  @Test
  public void canSelectSamplesByRegex() throws Exception {

    GenomicCoords gc = new GenomicCoords("1:572807-755079", 80, null, null);
    TrackVCF vcf = new TrackVCF("test_data/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz", gc);
    List<VCFFeature> linf = vcf.getIntervalFeatureList();

    GenotypeMatrix gm = new GenotypeMatrix();
    gm.setSelectSampleRegex("96|99");
    // gm.makeMatrix(linf, 80, null);

    String x = gm.printToScreen(true, linf, 80, null);
    String[] rows = x.split("\n");
    assertEquals(2, rows.length);

    // Exclude all samples
    gm.setSelectSampleRegex("^$");
    // gm.makeMatrix(linf, 80, null);
    assertTrue(gm.printToScreen(true, linf, 80, null).isEmpty());
  }

  @Test
  public void behaviourWithHaploid() {
    // TODO
  }

  @Test
  public void canHandleMissingGenotypes() {
    // TODO
  }

  @Test
  public void canHandleIntervalWithNoFeatures()
      throws IOException,
          InvalidGenomicCoordsException,
          IOException,
          ClassNotFoundException,
          InvalidRecordException,
          SQLException,
          InvalidColourException {
    GenomicCoords gc = new GenomicCoords("2:1-1000", 80, null, null);
    TrackVCF vcf = new TrackVCF("test_data/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz", gc);
    List<VCFFeature> linf = vcf.getIntervalFeatureList();

    GenotypeMatrix gm = new GenotypeMatrix();
    // gm.makeMatrix(linf, 80, null);
    String x = gm.printToScreen(true, linf, 80, null);
    assertTrue(x.isEmpty());
  }

  @Test
  public void canHandleVCFWithNoSamples() throws Exception {
    GenomicCoords gc = new GenomicCoords("1:1-10000000", 80, null, null);
    TrackVCF vcf = new TrackVCF("test_data/CHD.exon.2010_03.sites.vcf.gz", gc);
    List<VCFFeature> linf = vcf.getIntervalFeatureList();
    assertTrue(
        linf.size() > 0); // Make sure we do have some features otherwise the test is meaningless.

    GenotypeMatrix gm = new GenotypeMatrix();
    String x = gm.printToScreen(true, linf, 80, null);
    assertTrue(x.isEmpty());
  }

  @Test
  public void sampleOrderIsTheSameAsInVcf() throws Exception {
    GenomicCoords gc = new GenomicCoords("1:199982-200052", 70, null, null);
    TrackVCF vcf = new TrackVCF("test_data/sample_order.vcf", gc);
    List<VCFFeature> linf = vcf.getIntervalFeatureList();

    GenotypeMatrix gm = new GenotypeMatrix();
    String x = gm.printToScreen(true, linf, 80, vcf.getVcfHeader());

    assertTrue(x.split("\n")[0].startsWith("HG03"));
    assertTrue(x.split("\n")[1].startsWith("HG01"));
    assertTrue(x.split("\n")[2].startsWith("HG02"));
  }
}
