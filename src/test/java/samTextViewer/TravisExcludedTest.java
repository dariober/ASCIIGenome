package samTextViewer;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import java.io.File;
import java.io.IOException;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;
import org.junit.Test;

/** Tests that for some strange reason fail on Travis */
public class TravisExcludedTest {

  static SamReaderFactory srf = SamReaderFactory.make();
  static SamReader samReader = srf.open(new File("test_data/ds051.short.bam"));
  public static SAMSequenceDictionary samSeqDict =
      samReader.getFileHeader().getSequenceDictionary();

  public static String fastaFile = "test_data/chr7.fa";

  @Test
  public void canTestForExistingURLFile() {
    String urlStr =
        "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeHaibTfbs/wgEncodeHaibTfbsA549Atf3V0422111Etoh02PkRep1.broadPeak.gz";
    assertTrue(Utils.urlFileExists(urlStr));
    assertFalse(Utils.urlFileExists(urlStr + "foobar"));

    assertTrue(
        Utils.urlFileExists(
            "ftp://ftp.ensembl.org/pub/release-84/gtf/homo_sapiens/Homo_sapiens.GRCh38.84.abinitio.gtf.gz"));
    assertFalse(
        Utils.urlFileExists("ftp://ftp.ensembl.org/pub/release-84/gtf/homo_sapiens/foobar"));
  }

  @Test
  public void canGlobFiles() throws IOException {
    ArrayList<String> cmdInput = Utils.tokenize("test_data/ear*{bam,tdf} README.*", " ");
    List<String> globbed = Utils.globFiles(cmdInput);
    cmdInput =
        Utils.tokenize(
            "ftp://ftp.ensembl.org/pub/release-84/gtf/homo_sapiens/Homo_sapiens.GRCh38.84.abinitio.gtf.gz",
            " ");
    globbed = Utils.globFiles(cmdInput);
    assertEquals(1, globbed.size());
  }

  @Test
  public void canReformatFileName() {
    String x = "License.md";
    assertEquals("License.md", Utils.reformatFileName(x, false));
    assertTrue(Utils.reformatFileName(x, true).length() > x.length());

    x = Paths.get(System.getProperty("user.dir"), "/test_data/ds051.actb.bam").toString();
    assertEquals("test_data/ds051.actb.bam", Utils.reformatFileName(x, false));
    assertEquals(x, Utils.reformatFileName(x, true));

    x = "test_data/../test_data/ds051.actb.bam";
    assertEquals("test_data/ds051.actb.bam", Utils.reformatFileName(x, false));

    x = "test_data/../test_data/foobar"; // Non existent file
    assertEquals("test_data/foobar", Utils.reformatFileName(x, false));

    // URLs
    x =
        "http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeAwgTfbsUniform/wgEncodeAwgTfbsBroadDnd41CtcfUniPk.narrowPeak.gz";
    assertEquals(x, Utils.reformatFileName(x, false));

    x = "ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr18.fa.gz";
    assertEquals(x, Utils.reformatFileName(x, false));
    assertEquals(x, Utils.reformatFileName(x, true));

    x = x + "foobar";
    assertEquals(x, Utils.reformatFileName(x, true));
  }
}
