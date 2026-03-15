package aligner;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import colouring.Config;
import colouring.ConfigKey;
import com.google.common.base.Joiner;
import com.google.common.base.Splitter;
import exceptions.InvalidColourException;
import exceptions.InvalidConfigException;
import exceptions.InvalidGenomicCoordsException;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Stream;
import org.broad.igv.util.FileUtils;
import org.junit.Test;
import samTextViewer.GenomicCoords;
import samTextViewer.Utils;
import tracks.TrackReads;

public class SassyTest {

  private final Path sassy = Path.of("test_data/sassy/sassy-x86_64-unknown-linux-gnu");

  @Test
  public void testPattern() throws Exception {
    Path tempWorkDir = Utils.createTempDir("tmp.sassy.", true);

    GenomicCoords gc = new GenomicCoords("chr7:20001-25000", 80, null, "test_data/chr7.fa");

    Sassy sassy = new Sassy(gc, tempWorkDir);
    sassy.setExecPath(this.sassy);
    sassy.search(Splitter.on(" ").splitToList("-p AGRTGA -k 1"));
    sassy.writeSAMFile(Paths.get(sassy.getWorkDir().toString(), "tmp.sam"));

    // Test fasta reference is correct
    String fastaRef = sassy.getSamFileHeader().getProgramRecord("sassy").getCommandLine().replaceAll(".* ", "");
    String fastaSeq = new String(FileUtils.readFully(fastaRef));
    assertTrue(fastaSeq.startsWith(">chr7\nTGCCAGGAAAGCAGACACATCATCAAAATCCATTACAGAGGCTATAGTTCAGCCAAAGCT"));
    assertTrue(fastaSeq.endsWith("agacagccatcaaggaggtg\n"));
  }

  @Test
  public void testPatternFastaDupNames() throws Exception {
    Path tempWorkDir = Utils.createTempDir("tmp.sassy.", true);

    GenomicCoords gc = new GenomicCoords("chr7:10000-20000", 80, null, "test_data/chr7.fa");

    Sassy sassy = new Sassy(gc, tempWorkDir);
    sassy.setExecPath(this.sassy);
    boolean pass = false;
    try {
      sassy.search(Splitter.on(" ").splitToList("-f test_data/pattern_dups.fa.gz -k 1"));
    } catch (Exception e) {
      assertEquals("Duplicate FASTA sequence name: q1", e.getMessage().trim());
      pass = true;
    }
    assertTrue(pass);
  }

  @Test
  public void testPatternFasta() throws Exception {
    Path tempWorkDir = Utils.createTempDir("tmp.sassy.", true);

    GenomicCoords gc = new GenomicCoords("chr7:10000-20000", 80, null, "test_data/chr7.fa");

    Sassy sassy = new Sassy(gc, tempWorkDir);
    sassy.setExecPath(this.sassy);
    sassy.search(Splitter.on(" ").splitToList("-f test_data/pattern.fa.gz -k 1"));

    assertTrue(sassy.getSamFileHeader().getProgramRecord("sassy").getCommandLine().contains("sassy search -f test_data/pattern.fa.gz -k 1"));

    Path outsam = Paths.get(sassy.getWorkDir().toString(), "tmp.sam");
    sassy.writeSAMFile(outsam);
    String sam = new String(FileUtils.readFully(outsam.toString()));
    assertTrue(sam.contains("q1\t0\tchr7\t11251\t255\t50=\t*\t0\t0\tAAGAGGGCTACATTATTTATGAAACAGATACTGTTAACTCAGTCACCAGA\t*\tNM:i:0"));
  }

  @Test
  public void testPatternFile() throws Exception {
    Path tempWorkDir = Utils.createTempDir("tmp.sassy.", true);

    GenomicCoords gc = new GenomicCoords("chr7:10000-20000", 80, null, "test_data/chr7.fa");

    Sassy sassy = new Sassy(gc, tempWorkDir);
    sassy.setExecPath(this.sassy);
    sassy.search(Splitter.on(" ").splitToList("-l test_data/pattern.txt -k 1"));

    Path outsam = Paths.get(sassy.getWorkDir().toString(), "tmp.sam");
    sassy.writeSAMFile(outsam);
    String sam = new String(FileUtils.readFully(outsam.toString()));
    assertTrue(sam.contains("1\t0\tchr7\t11251\t255\t50=\t*\t0\t0\tAAGAGGGCTACATTATTTATGAAACAGATACTGTTAACTCAGTCACCAGA\t*\tNM:i:0"));
  }

  @Test
  public void testConfigPath()
      throws IOException, InvalidConfigException, InvalidColourException, InvalidGenomicCoordsException, InterruptedException {
    new Config(null);
    Config.set(ConfigKey.sassy, this.sassy.toString());

    Path tempWorkDir = Utils.createTempDir("tmp.sassy.", true);
    GenomicCoords gc = new GenomicCoords("chr7:10000-20000", 80, null, "test_data/chr7.fa");

    Sassy sassy = new Sassy(gc, tempWorkDir);
    sassy.setExecPath(Paths.get(Config.get(ConfigKey.sassy)));
    sassy.search(Splitter.on(" ").splitToList("-p AAGRTGAG -k 1"));

  }
}
