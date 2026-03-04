package sortBgzipIndex;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import com.google.common.io.Files;
import exceptions.InvalidGenomicCoordsException;
import exceptions.InvalidRecordException;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.tribble.index.tabix.TabixFormat;
import htsjdk.tribble.readers.TabixReader.Iterator;
import htsjdk.tribble.util.TabixUtils;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import org.junit.Test;
import samTextViewer.Utils;
import utils.FlexibleTabixReader;

public class MakeTabixFileTest {

  public void vcfTester(String inVcf) {
    // Credit: https://www.biostars.org/p/262943/
    VCFFileReader r = new VCFFileReader(new File(inVcf));
    CloseableIterator<VariantContext> t = r.iterator();
    while (t.hasNext()) {
      t.next();
    }
    t.close();
    r.close();
  }

  @Test
  public void testSpeed() throws SQLException, IOException, ClassNotFoundException {
    String infile = "/home/dario/Downloads/schistosoma_mansoni.PRJEA36577.WBPS16.annotations.gff3";
    File outfile = new File("tmp.gff.gz");
    // outfile.deleteOnExit();

    File expectedTbi = new File(outfile.getAbsolutePath() + TabixUtils.STANDARD_INDEX_EXTENSION);
    // expectedTbi.deleteOnExit();
    long t0 = System.currentTimeMillis();
    new MakeTabixIndex(infile, outfile, TabixFormat.GFF, '\t');
    long t1 = System.currentTimeMillis();
    System.out.println("DONE " + (t1 -t0)/1000.0);
  }

  @Test
  public void canIndexGtf()
      throws ClassNotFoundException, IOException, InvalidRecordException, SQLException {
    String infile = "test_data/hg19_genes.gtf.gz";
    File outfile = new File("test_data/hg19_genes.gtf.tmp.gz");
    outfile.deleteOnExit();

    File expectedTbi = new File(outfile.getAbsolutePath() + TabixUtils.STANDARD_INDEX_EXTENSION);
    expectedTbi.deleteOnExit();
    new MakeTabixIndex(infile, outfile, TabixFormat.GFF, '\t');
  }

  @Test
  public void canIndexCustomSep()
      throws ClassNotFoundException, IOException, InvalidRecordException, SQLException {
    String infile = "test_data/generic.csv";
    File outfile = new File("test_data/generic.csv.gz");
    outfile.deleteOnExit();

    File expectedTbi = new File(outfile.getAbsolutePath() + TabixUtils.STANDARD_INDEX_EXTENSION);
    expectedTbi.deleteOnExit();
    TabixFormat tabixFormat = new TabixFormat(TabixFormat.GENERIC_FLAGS, 5, 1, 2, '#', 1);
    new MakeTabixIndex(infile, outfile, tabixFormat, ',');
    try (BufferedReader br = Utils.reader(outfile.getAbsolutePath())) {
      assertEquals("# header", br.readLine());
      assertEquals("start,end,V1,V2,chrom,V3", br.readLine());
      assertEquals("5566777,5566778,1,0.1,chr7,-99", br.readLine());
    }

    infile = "test_data/generic.tsv";
    outfile = new File("test_data/generic.tsv.gz");
    outfile.deleteOnExit();
    expectedTbi = new File(outfile.getAbsolutePath() + TabixUtils.STANDARD_INDEX_EXTENSION);
    expectedTbi.deleteOnExit();
    tabixFormat = new TabixFormat(TabixFormat.GENERIC_FLAGS, 5, 1, 2, '#', 1);
    new MakeTabixIndex(infile, outfile, tabixFormat, '\t');
    try (BufferedReader br = Utils.reader(outfile.getAbsolutePath())) {
      assertEquals("# header", br.readLine());
      assertEquals("start\tend\tV1\tV2\tchrom\tV3", br.readLine());
      assertEquals("5566777\t5566778\t1\t0.1\tchr7\t-99", br.readLine());
    }

    infile = "test_data/generic.txt";
    outfile = new File("test_data/generic.txt.gz");
    outfile.deleteOnExit();
    expectedTbi = new File(outfile.getAbsolutePath() + TabixUtils.STANDARD_INDEX_EXTENSION);
    expectedTbi.deleteOnExit();
    tabixFormat = new TabixFormat(TabixFormat.GENERIC_FLAGS, 5, 1, 2, '#', 1);
    new MakeTabixIndex(infile, outfile, tabixFormat, ' ');
    try (BufferedReader br = Utils.reader(outfile.getAbsolutePath())) {
      assertEquals("# header", br.readLine());
      assertEquals("start end V1 V2 chrom V3", br.readLine());
      assertEquals("5566777 5566778 1 0.1 chr7 -99", br.readLine());
    }

    infile = "test_data/test.bedGraph";
    outfile = new File("test_data/test.bedGraph.tmp.gz");
    outfile.deleteOnExit();
    expectedTbi = new File(outfile.getAbsolutePath() + TabixUtils.STANDARD_INDEX_EXTENSION);
    expectedTbi.deleteOnExit();
    tabixFormat = new TabixFormat(TabixFormat.ZERO_BASED, 1, 2, 3, '#', 0);
    new MakeTabixIndex(infile, outfile, tabixFormat, '\t');
    FlexibleTabixReader tr = new FlexibleTabixReader(outfile.getAbsolutePath());
    Iterator qry = tr.query("chr1:1-20");
    assertEquals("chr1\t0\t1\t1\t0", qry.next());
    assertEquals("chr1\t5\t10\t-1\t1", qry.next());
    assertEquals("chr1\t15\t20\t5\t2", qry.next());
    assertNull(qry.next());
  }

  @Test
  public void canQueryCustomSepZeroBased()
      throws ClassNotFoundException, IOException, SQLException {
    String infile = "test_data/generic2.csv";
    File outfile = new File("test_data/generic2.csv.tmp.gz");
    outfile.deleteOnExit();
    File expectedTbi = new File(outfile.getAbsolutePath() + TabixUtils.STANDARD_INDEX_EXTENSION);
    expectedTbi.deleteOnExit();
    TabixFormat tabixFormat = new TabixFormat(TabixFormat.ZERO_BASED, 3, 1, 2, '#', 1);
    new MakeTabixIndex(infile, outfile, tabixFormat, ',');
    FlexibleTabixReader tr = new FlexibleTabixReader(outfile.getAbsolutePath());
    tr.setColumnSeparator(',');
    Iterator qry = tr.query("chr1:1-20");
    assertEquals("0,1,chr1,1,0", qry.next());
    assertEquals("5,10,chr1,-1,1", qry.next());
    assertEquals("15,20,chr1,5,2", qry.next());
    assertNull(qry.next());
  }

  @Test
  public void canReadSpaceSeparatedBed()
      throws ClassNotFoundException, IOException, InvalidRecordException, SQLException {
    String infile = "test_data/space_sep.bedgraph";
    File outfile = new File("test_data/tmp.bedgraph.gz");
    outfile.deleteOnExit();

    File expectedTbi = new File(outfile.getAbsolutePath() + TabixUtils.STANDARD_INDEX_EXTENSION);
    expectedTbi.deleteOnExit();

    TabixFormat tabixFormat = new TabixFormat(TabixFormat.GENERIC_FLAGS, 1, 2, 3, '#', 0);
    new MakeTabixIndex(infile, outfile, tabixFormat, ' ');

    try (BufferedReader br = Utils.reader(outfile.getAbsolutePath())) {
      assertEquals("# comment line", br.readLine());
      assertEquals("chr1 0 10 1", br.readLine());
    }
  }

  @Test
  public void canReadSpaceSeparatedSorted()
      throws ClassNotFoundException, IOException, InvalidRecordException, SQLException {
    String infile = "test_data/space_sorted.bedgraph";
    File outfile = new File("test_data/tmp.bedgraph.gz");
    outfile.deleteOnExit();

    File expectedTbi = new File(outfile.getAbsolutePath() + TabixUtils.STANDARD_INDEX_EXTENSION);
    expectedTbi.deleteOnExit();

    TabixFormat tabixFormat = new TabixFormat(TabixFormat.ZERO_BASED, 1, 2, 3, '#', 0);
    new MakeTabixIndex(infile, outfile, tabixFormat, ' ');

    try (FlexibleTabixReader reader = new FlexibleTabixReader(outfile.getAbsolutePath())) {
      reader.setColumnSeparator(' ');
      Iterator iter = reader.query("chr1:1-21");
      assertEquals("chr1 0 10 1", iter.next());
      assertEquals("chr1 20 30 2", iter.next());
      assertNull(iter.next());
    }
  }

  @Test
  public void canReadSpaceSeparatedOnelineInput()
      throws ClassNotFoundException, IOException, SQLException {
    String infile = "test_data/space_oneline.bedgraph";
    File outfile = new File("test_data/tmp.bedgraph.gz");
    outfile.deleteOnExit();

    File expectedTbi = new File(outfile.getAbsolutePath() + TabixUtils.STANDARD_INDEX_EXTENSION);
    expectedTbi.deleteOnExit();
    TabixFormat tabixFormat = new TabixFormat(TabixFormat.ZERO_BASED, 1, 2, 3, '#', 0);
    new MakeTabixIndex(infile, outfile, tabixFormat, ' ');
    try (FlexibleTabixReader reader = new FlexibleTabixReader(outfile.getAbsolutePath())) {
      reader.setColumnSeparator(' ');
      Iterator iter = reader.query("chr1");
      assertEquals("chr1 0 10 1", iter.next());
      assertNull(iter.next());
    }
  }

  @Test
  public void canReadGffWithFasta()
      throws ClassNotFoundException, IOException, InvalidRecordException, SQLException {
    String infile = "test_data/features_fasta.gff";
    File outfile = new File("test_data/tmp.gff.gz");
    outfile.deleteOnExit();

    File expectedTbi = new File(outfile.getAbsolutePath() + TabixUtils.STANDARD_INDEX_EXTENSION);
    expectedTbi.deleteOnExit();

    new MakeTabixIndex(infile, outfile, TabixFormat.GFF, '\t');
    BufferedReader br = Utils.reader(outfile.getAbsolutePath());
    assertEquals("##FASTA", br.readLine());
    assertEquals("# bla bla", br.readLine());
    assertTrue(br.readLine().startsWith("1\tASCIIGenome\tbookmark\t1\t"));
  }

  @Test
  public void canHandleBgzExtension() throws Exception {
    // This effectively tests IOUtils
    String infile = "test_data/bgz_noindex.vcf.bgz";
    File outfile = new File("test_data/tmp.vcf.gz");
    outfile.deleteOnExit();

    File expectedTbi = new File(outfile.getAbsolutePath() + TabixUtils.STANDARD_INDEX_EXTENSION);
    expectedTbi.deleteOnExit();

    new MakeTabixIndex(infile, outfile, TabixFormat.VCF, '\t');
    vcfTester(outfile.getCanonicalPath());
  }

  @Test
  public void canCompressAndIndexHeaderlessVCF()
      throws ClassNotFoundException, IOException, InvalidRecordException, SQLException {

    String infile = "test_data/noheader.vcf";
    File outfile = new File("test_data/noheader.vcf.gz");
    outfile.deleteOnExit();

    File expectedTbi = new File(outfile.getAbsolutePath() + TabixUtils.STANDARD_INDEX_EXTENSION);
    expectedTbi.deleteOnExit();

    new MakeTabixIndex(infile, outfile, TabixFormat.VCF, '\t');

    assertTrue(outfile.exists());
    assertTrue(outfile.length() > 200);
    assertTrue(expectedTbi.exists());
    assertTrue(expectedTbi.length() > 100);

    FlexibleTabixReader tbx = new FlexibleTabixReader(outfile.getAbsolutePath());
    Iterator x = tbx.query("1", 1, 10000000);
    assertTrue(x.next().startsWith("1"));
  }

  @Test
  public void testRealFileSizeVCF()
      throws ClassNotFoundException, IOException, InvalidRecordException, SQLException {

    // See test_data/README.md for this file. This is fairly large and we want to check it is
    // processed
    // in a reasonable amount of time.
    String infile = "test_data/ALL.wex.union_illumina_wcmc_bcm_bc_bi.20110521.snps.exome.sites.vcf";

    File outfile = new File("deleteme.vcf.gz");
    outfile.deleteOnExit();
    File expectedTbi = new File(outfile.getAbsolutePath() + TabixUtils.STANDARD_INDEX_EXTENSION);
    expectedTbi.deleteOnExit();

    long t0 = System.currentTimeMillis();
    new MakeTabixIndex(infile, outfile, TabixFormat.VCF, '\t');
    long t1 = System.currentTimeMillis();

    assertTrue(outfile.exists());
    assertTrue(outfile.length() > 1000);
    assertTrue((t1 - t0) < 20000); // Should be << than 20 sec, ~2 sec

    // Check you can read ok
    this.vcfTester(outfile.getAbsolutePath());

    try (FlexibleTabixReader tr = new FlexibleTabixReader(outfile.getAbsolutePath())) {
      Iterator qry = tr.query("1:69270-69270");
      assertEquals("1\t69270\t.\tA\tG\t.\tFAIL\tCalledBy=WCMC,BCM,BC", qry.next());
      assertNull(qry.next());
    }
    try (FlexibleTabixReader tr = new FlexibleTabixReader(outfile.getAbsolutePath())) {
      Iterator qry = tr.query("1:69270-69428");
      assertEquals("1\t69270\t.\tA\tG\t.\tFAIL\tCalledBy=WCMC,BCM,BC", qry.next());
      assertEquals("1\t69428\t.\tT\tG\t.\tFAIL\tCalledBy=BCM", qry.next());
      assertNull(qry.next());
    }
    try (FlexibleTabixReader tr = new FlexibleTabixReader(outfile.getAbsolutePath())) {
      Iterator qry = tr.query("1:69270-69428");
      assertEquals("1\t69270\t.\tA\tG\t.\tFAIL\tCalledBy=WCMC,BCM,BC", qry.next());
      assertEquals("1\t69428\t.\tT\tG\t.\tFAIL\tCalledBy=BCM", qry.next());
      assertNull(qry.next());
    }
    try (FlexibleTabixReader tr = new FlexibleTabixReader(outfile.getAbsolutePath())) {
      Iterator qry = tr.query("Y:22942929");
      assertEquals("Y\t22942929\t.\tT\tG\t.\tFAIL\tCalledBy=BI", qry.next());
      assertNull(qry.next());
    }
  }

  @Test
  public void canCompressAndIndexVCF_CEU()
      throws ClassNotFoundException, IOException, InvalidRecordException, SQLException {

    String infile = "test_data/CEU.exon.2010_06.genotypes.vcf";

    File outfile = new File("deleteme.vcf.gz");
    outfile.deleteOnExit();
    File expectedTbi = new File(outfile.getAbsolutePath() + TabixUtils.STANDARD_INDEX_EXTENSION);
    expectedTbi.deleteOnExit();

    new MakeTabixIndex(infile, outfile, TabixFormat.VCF);

    assertTrue(outfile.exists());
    assertTrue(outfile.length() > 1000);

    // Check you can read ok
    this.vcfTester(outfile.getAbsolutePath());

    try (FlexibleTabixReader tr = new FlexibleTabixReader(outfile.getAbsolutePath())) {
      Iterator qry = tr.query("1:113054374-113054374");
      assertTrue(qry.next().startsWith("1\t113054374\t.\tCTTG\tC\t23\tPASS"));
      assertNull(qry.next());
    }
  }

  // @Test
  public void handlingInvalidFile()
      throws ClassNotFoundException, IOException, InvalidRecordException, SQLException {

    String infile = "test_data/chr7.fa";
    File outfile = new File("deleteme.gtf.gz");
    outfile.deleteOnExit();
    File expectedTbi = new File(outfile.getAbsolutePath() + TabixUtils.STANDARD_INDEX_EXTENSION);
    expectedTbi.deleteOnExit();

    new MakeTabixIndex(infile, outfile, TabixFormat.BED, '\t');
  }

  @Test
  public void canHandleEmptyFile()
      throws ClassNotFoundException,
          IOException,
          InvalidRecordException,
          SQLException,
          InvalidGenomicCoordsException {
    String infile = "test_data/empty.bed";
    File outfile = new File("test_data/empty.tmp.bed.gz");
    outfile.deleteOnExit();
    File expectedTbi = new File(outfile.getAbsolutePath() + TabixUtils.STANDARD_INDEX_EXTENSION);
    expectedTbi.deleteOnExit();

    new MakeTabixIndex(infile, outfile, TabixFormat.BED, '\t');

    assertTrue(outfile.exists());
    assertTrue(expectedTbi.exists());

    assertEquals(28, outfile.length()); // Checked with `bgzip empty.bed && ls -l empty.bed.gz`
    assertTrue(expectedTbi.length() > 0);
  }

  @Test
  public void canCompressAndIndexSortedFile()
      throws IOException, InvalidRecordException, ClassNotFoundException, SQLException {
    String infile = "test_data/overlapped.bed";
    File outfile = new File("test_data/tmp.bed.gz");
    outfile.deleteOnExit();

    File expectedTbi = new File(outfile.getAbsolutePath() + TabixUtils.STANDARD_INDEX_EXTENSION);
    expectedTbi.deleteOnExit();

    new MakeTabixIndex(infile, outfile, TabixFormat.BED, '\t');

    assertTrue(outfile.exists());
    assertTrue(outfile.length() > 80);
    assertTrue(expectedTbi.exists());
    assertTrue(expectedTbi.length() > 80);

    FlexibleTabixReader tbx = new FlexibleTabixReader(outfile.getAbsolutePath());
    Iterator x = tbx.query("chr1", 1, 1000000);
    assertTrue(x.next().startsWith("chr1"));
  }

  @Test
  public void canCompressAndIndexSortedGzipFile()
      throws IOException, InvalidRecordException, ClassNotFoundException, SQLException {

    String infile = "test_data/hg19_genes.gtf.gz";
    File outfile = new File("test_data/tmp2.bed.gz");
    outfile.deleteOnExit();

    File expectedTbi = new File(outfile.getAbsolutePath() + TabixUtils.STANDARD_INDEX_EXTENSION);
    expectedTbi.deleteOnExit();

    new MakeTabixIndex(infile, outfile, TabixFormat.GFF, '\t');

    assertTrue(outfile.exists());
    assertTrue(outfile.length() > 7000000);
    assertTrue(expectedTbi.exists());

    FlexibleTabixReader tbx = new FlexibleTabixReader(outfile.getAbsolutePath());
    Iterator x = tbx.query("chr1", 1, 1000000);
    assertTrue(x.next().startsWith("chr1"));
  }

  @Test
  public void canCompressAndIndexUnsortedFile()
      throws IOException, InvalidRecordException, ClassNotFoundException, SQLException {

    String infile = "test_data/refSeq.hg19.bed.gz";
    File outfile = new File("test_data/tmp3.bed.gz");
    outfile.deleteOnExit();

    File expectedTbi = new File(outfile.getAbsolutePath() + TabixUtils.STANDARD_INDEX_EXTENSION);
    expectedTbi.deleteOnExit();

    new MakeTabixIndex(infile, outfile, TabixFormat.BED, '\t');

    assertTrue(outfile.exists());
    assertTrue(outfile.length() > 200000);
    assertTrue(expectedTbi.exists());
    assertTrue(expectedTbi.length() > 100000);
  }

  @Test
  public void canCompressAndIndexSortedURL()
      throws IOException, InvalidRecordException, ClassNotFoundException, SQLException {

    String infile =
        "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeDacMapabilityConsensusExcludable.bed.gz";
    File outfile = new File("test_data/tmp4.bed.gz");
    outfile.deleteOnExit();

    File expectedTbi = new File(outfile.getAbsolutePath() + TabixUtils.STANDARD_INDEX_EXTENSION);
    expectedTbi.deleteOnExit();

    new MakeTabixIndex(infile, outfile, TabixFormat.BED, '\t');

    assertTrue(outfile.exists());
    assertTrue(outfile.length() > 1000);
    assertTrue(expectedTbi.exists());
    assertTrue(expectedTbi.length() > 1000);
  }

  @Test
  public void canCompressAndIndexUnsortedURL()
      throws IOException, InvalidRecordException, ClassNotFoundException, SQLException {

    String infile =
        "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeSydhTfbs/wgEncodeSydhTfbsGm12878P300bStdPk.narrowPeak.gz";
    File outfile = new File("test_data/tmp5.bed.gz");
    outfile.deleteOnExit();

    File expectedTbi = new File(outfile.getAbsolutePath() + TabixUtils.STANDARD_INDEX_EXTENSION);
    expectedTbi.deleteOnExit();

    new MakeTabixIndex(infile, outfile, TabixFormat.BED, '\t');

    assertTrue(outfile.exists());
    assertTrue(outfile.length() > 1000);
    assertTrue(expectedTbi.exists());
    assertTrue(expectedTbi.length() > 1000);
  }

  @Test
  public void canCompressAndIndexVCF()
      throws ClassNotFoundException, IOException, InvalidRecordException, SQLException {

    String infile = "test_data/CHD.exon.2010_03.sites.unsorted.vcf";
    File outfile = new File("test_data/tmp6.bed.gz");
    outfile.deleteOnExit();

    File expectedTbi = new File(outfile.getAbsolutePath() + TabixUtils.STANDARD_INDEX_EXTENSION);
    expectedTbi.deleteOnExit();

    new MakeTabixIndex(infile, outfile, TabixFormat.VCF, '\t');

    assertTrue(outfile.exists());
    assertTrue(outfile.length() > 1000);
    assertTrue(expectedTbi.exists());
    assertTrue(expectedTbi.length() > 1000);

    FlexibleTabixReader tbx = new FlexibleTabixReader(outfile.getAbsolutePath());
    Iterator x = tbx.query("1", 20000000, 30000000);
    assertTrue(x.next().startsWith("1"));

    // Check you can read ok
    this.vcfTester(outfile.getAbsolutePath());
  }

  @Test
  public void blockCompressFileInPlace()
      throws IOException, ClassNotFoundException, InvalidRecordException, SQLException {
    // Test we can compress and index a file and overwrite the original file.

    File testFile = new File("deleteme.bed.gz");
    testFile.deleteOnExit();
    Files.copy(new File("test_data/refSeq.hg19.bed.gz"), testFile);

    File expectedTbi = new File(testFile.getAbsolutePath() + TabixUtils.STANDARD_INDEX_EXTENSION);
    expectedTbi.deleteOnExit();

    new MakeTabixIndex(testFile.getAbsolutePath(), testFile, TabixFormat.BED, '\t');

    assertTrue(testFile.exists());
    assertTrue(testFile.length() > 200000);
    assertTrue(expectedTbi.exists());
    assertTrue(expectedTbi.length() > 100000);
  }

  @Test
  public void canProcessBedgraphWithTrackLine()
      throws ClassNotFoundException, IOException, InvalidRecordException, SQLException {

    String testFile = "test_data/test2.bedGraph";

    String outfile = testFile + ".tmp.bedGraph.gz";
    File expectedTbi = new File(outfile + TabixUtils.STANDARD_INDEX_EXTENSION);
    new File(outfile).deleteOnExit();
    expectedTbi.deleteOnExit();

    new MakeTabixIndex(testFile, new File(outfile), TabixFormat.BED, '\t');
    assertTrue(expectedTbi.exists());
    assertTrue(expectedTbi.length() > 50);
  }
}
