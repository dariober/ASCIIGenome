package tracks;

import static org.junit.Assert.*;

import colouring.Config;
import colouring.Xterm256;
import com.google.common.base.Splitter;
import exceptions.InvalidColourException;
import exceptions.InvalidCommandLineException;
import exceptions.InvalidConfigException;
import exceptions.InvalidGenomicCoordsException;
import exceptions.InvalidRecordException;

import java.io.IOException;
import java.sql.SQLException;
import java.util.List;
import java.util.regex.Pattern;

import org.junit.Before;
import org.junit.Test;
import samTextViewer.GenomicCoords;

public class TrackVCFTest {

  @Before
  public void prepareConfig() throws IOException, InvalidConfigException {
    new Config(null);
    new Xterm256();
  }

  @Test
  public void canReadBcf() throws IOException, InvalidGenomicCoordsException, SQLException, ClassNotFoundException, InvalidRecordException {
    GenomicCoords gc = new GenomicCoords("chr1:1-12160", 80, null, null);
    TrackVCF tif = new TrackVCF("test_data/gnomad.exomes.v4.1.sites.chr1.bcf", gc);
    assertEquals(8, tif.getFeatureList().size());
    assertTrue(tif.getFeatureList().get(0).getRaw().startsWith("chr1\t11994\t.\tT\tC\t.\tAC0;AS_VQSR\tAC=0;AC_XX=0;AC_XY=0;"));
    assertEquals(24, tif.getReader().getChromosomes().size());
  }

  //  @Test
//  public void canReadBcfWithoutIndex() throws SQLException, InvalidGenomicCoordsException, IOException, ClassNotFoundException, InvalidRecordException {
//    GenomicCoords gc = new GenomicCoords("chr1:1-12160", 80, null, null);
//    TrackVCF tif = new TrackVCF("test_data/gnomad.exomes.v4.1.sites.chr1.no_index.bcf", gc);
//  }

  @Test
  public void canReadOddFilename()
      throws InvalidGenomicCoordsException,
          IOException,
          ClassNotFoundException,
          InvalidRecordException,
          SQLException {
    GenomicCoords gc = new GenomicCoords("1:1-100000", 80, null, null);
    TrackVCF tif = new TrackVCF("test_data/odd[filename].vcf.gz", gc);
    assertEquals("1", (tif.getGc().getChrom()));
    assertEquals(22, tif.getReader().getChromosomes().size());
    tif.close();
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
    TrackVCF tif = new TrackVCF(bgzFn, gc);
    assertEquals(3, tif.getFeatureList().size());
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
    TrackVCF tif =
        new TrackVCF(
            "https://raw.githubusercontent.com/dariober/ASCIIGenome/master/test_data/CHD.exon.2010_03.sites.unsorted.vcf",
            gc);
    assertEquals("http", tif.getFilename().substring(0, 4));
    assertEquals(3, tif.getFeatureList().size());
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
    TrackVCF tif = new TrackVCF(bgzFn, gc);
    assertEquals(3, tif.getFeatureList().size());
    assertEquals("1	1105468	.	G	A	.	PASS	AA=g;AC=2;AN=174;DP=1238", tif.getFeatureList().get(0).getRaw());
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
    TrackVCF tif = new TrackVCF(intervalFileName, gc);
    assertFalse(tif.getFeatureList().isEmpty());

    // .bgz, with index
    intervalFileName = "test_data/bgz_index.vcf.bgz";
    tif = new TrackVCF(intervalFileName, gc);
    assertFalse(tif.getFeaturesInInterval("1", 1, 200000000).isEmpty());
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
    String fn = "test_data/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz";
    TrackVCF tif = new TrackVCF(fn, gc);
    tif.setNoFormat(true);
    assertTrue(tif.printToScreen().contains("HG00096"));
  }

  @Test
  public void canFindIndel()
      throws IOException,
          InvalidGenomicCoordsException,
          ClassNotFoundException,
          InvalidRecordException,
          SQLException {

    GenomicCoords gc = new GenomicCoords("1:113050000", 80, null, null);
    String fn = "test_data/CEU.exon.2010_06.genotypes.vcf.gz";
    TrackVCF tif = new TrackVCF(fn, gc);
    VCFFeature x = tif.findNextRegexInGenome(Pattern.compile(".*113054374.*"), "1", 113050000);
    assertTrue(x.getRaw().contains("\t113054374\t"));
    assertEquals(113054374, x.getFrom());
  }

  @Test
  public void canChangeFeatureDisplayMode()
      throws IOException,
      InvalidGenomicCoordsException,
      ClassNotFoundException,
      InvalidRecordException,
      SQLException, InvalidColourException {

    GenomicCoords gc = new GenomicCoords("1:632191-1685839", 200, null, null);
    String fn = "test_data/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf";
    TrackVCF tif = new TrackVCF(fn, gc);
    tif.setNoFormat(true);
    assertEquals(9, Splitter.on("\n").splitToList(tif.printToScreen()).size());

    tif.setFeatureDisplayMode(FeatureDisplayMode.COLLAPSED);
    assertEquals(6, Splitter.on("\n").splitToList(tif.printToScreen()).size());

    tif.setFeatureDisplayMode(FeatureDisplayMode.ONELINE);
    assertEquals(4, Splitter.on("\n").splitToList(tif.printToScreen()).size());
  }

  @Test
  public void canPrintNormalizedVcfLines()
      throws ClassNotFoundException,
          IOException,
          InvalidGenomicCoordsException,
          InvalidRecordException,
          SQLException,
          InvalidColourException,
          InvalidCommandLineException {

    GenomicCoords gc = new GenomicCoords("1:645709-645975", 80, null, null);
    TrackVCF tif = new TrackVCF("test_data/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz", gc);
    tif.setPrintMode(PrintRawLine.FULL);
    tif.setNoFormat(true);
    tif.setPrintNormalizedVcf(true);

    String out = tif.printLines();
    assertEquals(3, out.split("\n").length);
    assertTrue(out.contains(" HG00096 | GT"));

    // VCF with without samples
    gc = new GenomicCoords("1:1105467-1105647", 80, null, null);
    tif = new TrackVCF("test_data/CHD.exon.2010_03.sites.vcf.gz", gc);
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
    TrackVCF tvcf = new TrackVCF("test_data/vep.vcf", gc);
    tvcf.setPrintMode(PrintRawLine.FULL);
    tvcf.setNoFormat(true);
    String woVep = tvcf.printLines();
    tvcf.setPrintFormattedVep("");
    String printed = tvcf.printLines();
    assertTrue(Splitter.on("\n").splitToList(tvcf.printLines()).size() > 50);
    assertTrue(printed.contains("Consequence "));
    assertFalse(printed.contains("SWISSPROT"));

    // INFO tag not found: Do nothing
    tvcf.setPrintFormattedVep("csq_na");
    printed = tvcf.printLines();
    assertEquals(woVep, printed);

    // Only ask for some headers, case insensitive
    tvcf.setPrintFormattedVep("CSQ,ConseQUENCE,Allele");
    printed = tvcf.printLines();
    assertTrue(printed.contains("Consequence"));
    assertTrue(printed.contains("Allele"));
    assertFalse(printed.contains("IMPACT"));

    // Omit CSQ
    tvcf.setPrintFormattedVep("CSQ,null");
    printed = tvcf.printLines();
    assertTrue(printed.contains("CSQ=..."));

    // No effect without VEP tag
    gc = new GenomicCoords("1:1105468-34435998", 80, null, null);
    tvcf = new TrackVCF("test_data/CHD.exon.2010_03.sites.vcf", gc);
    tvcf.setPrintMode(PrintRawLine.FULL);
    tvcf.setNoFormat(true);
    tvcf.setPrintFormattedVep(null);
    woVep = tvcf.printLines();
    tvcf.setPrintFormattedVep("");
    String withVep = tvcf.printLines();
    assertEquals(woVep, withVep);

    // No effect on non-VCF track
    gc = new GenomicCoords("chr18:1-10000", 80, null, null);
    TrackIntervalFeature tif = new TrackIntervalFeature("test_data/refSeq.hg19.short.bed", gc);
    tif.setPrintMode(PrintRawLine.FULL);
    tif.setNoFormat(true);
    tif.setPrintFormattedVep(null);
    woVep = tif.printLines();
    tif.setPrintFormattedVep("");
    withVep = tif.printLines();
    assertEquals(woVep, withVep);
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
    TrackVCF tif = new TrackVCF("test_data/CHD.exon.2010_03.sites.vcf.gz", gc);

    List<VCFFeature> xset = tif.getFeaturesInInterval("1", 1, 10000000);
    assertEquals(9, xset.size());
    VCFFeature x = xset.get(1);
    assertEquals("1", x.getChrom());
    assertEquals(1108138, x.getFrom());
  }

    @Test
    public void canShowTitle()
        throws ClassNotFoundException,
        InvalidColourException,
        IOException,
        InvalidGenomicCoordsException,
        InvalidRecordException,
        SQLException {
        // See issue #42
        String intervalFileName = "test_data/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz";
        GenomicCoords gc = new GenomicCoords("1:200000-200261", 80, null, null);
        TrackVCF tif = new TrackVCF(intervalFileName, gc);
        tif.setNoFormat(true);
        assertTrue(tif.getTitle().contains("N: 2"));
    }

  @Test
  public void canReadUnsortedVCF()
      throws IOException,
          InvalidGenomicCoordsException,
          ClassNotFoundException,
          InvalidRecordException,
          SQLException {

    GenomicCoords gc = new GenomicCoords("chr18:1-10000", 80, null, null);
    TrackVCF tif = new TrackVCF("test_data/CHD.exon.2010_03.sites.unsorted.vcf", gc);
    List<VCFFeature> xset = tif.getFeaturesInInterval("1", 1, 10000000);
    assertEquals(9, xset.size());
    VCFFeature x = xset.get(1);
    assertEquals("1", x.getChrom());
    assertEquals(1108138, x.getFrom());
  }
}
