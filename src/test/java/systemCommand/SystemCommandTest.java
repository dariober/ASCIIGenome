package systemCommand;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.IOException;
import java.util.List;
import java.util.stream.Stream;
import org.junit.Test;

public class SystemCommandTest {
  @Test
  public void canExecSystemCommand() throws IOException, InterruptedException {
    String[] input = {"chr1\t1\t10\tfoo", "chr1\t1\t100\tbar", "chr1\t1\t200\tspam"};

    SystemCommand sysCmd = new SystemCommand();
    try (Stream<String> s =
        sysCmd.streamLinesThroughSystemCommand(
            Stream.of(input), null, "awk -v FOO=5 '$3 > FOO' | grep bar | cat")) {
      List<String> res = s.toList();
      assertEquals(1, res.size());
      assertEquals("chr1\t1\t100\tbar", res.get(0));
    }

    boolean pass = false;
    try (Stream<String> s =
        sysCmd.streamLinesThroughSystemCommand(Stream.of(input), null, "foobar")) {
    } catch (RuntimeException e) {
      pass = true;
    }
    assertTrue(pass);
  }

  @Test
  public void canStreamLinesThroughSystemCommand() throws IOException, InterruptedException {
    String[] input = {"chr1\t1\t10\tfoo", "chr1\t1\t100\tbar", "chr1\t1\t200\tspam"};

    SystemCommand sysCmd = new SystemCommand();
    sysCmd.deleteTempFile(); // Check this is harmless if tmp file is null
    try (Stream<String> stream =
        sysCmd.streamLinesThroughSystemCommand(
            Stream.of(input), "##foo\n##bar\n##chr1\t2\t300", "grep bar | sed 's/bar/BAR/'")) {
      List<String> res = stream.toList();
      assertEquals(2, res.size());
      assertEquals("##BAR", res.get(0));
      assertEquals("chr1\t1\t100\tBAR", res.get(1));
    }

    try (Stream<String> stream =
        sysCmd.streamLinesThroughSystemCommand(
            Stream.of(input), null, "grep bar | sed 's/bar/BAR/'")) {
      List<String> res = stream.toList();
      assertEquals(1, res.size());
      assertEquals("chr1\t1\t100\tBAR", res.get(0));
    }

    try (Stream<String> stream =
        sysCmd.streamLinesThroughSystemCommand(Stream.of(input), null, "")) {
      List<String> res = stream.toList();
      assertTrue(res.isEmpty());
    }

    try (Stream<String> stream =
        sysCmd.streamLinesThroughSystemCommand(Stream.of(new String[] {}), null, "")) {
      List<String> res = stream.toList();
      assertTrue(res.isEmpty());
    }

    try (Stream<String> stream =
        sysCmd.streamLinesThroughSystemCommand(Stream.of(new String[] {}), null, "cat")) {
      List<String> res = stream.toList();
      assertTrue(res.isEmpty());
    }

    try (Stream<String> stream =
        sysCmd.streamLinesThroughSystemCommand(Stream.of(new String[] {}), "##HEADER", "cat")) {
      List<String> res = stream.toList();
      assertEquals("##HEADER", res.get(0));
    }

    boolean pass = false;
    try (Stream<String> ignored =
        sysCmd.streamLinesThroughSystemCommand(
            Stream.of(input), "##foo\n##bar\n##chr1\t2\t300", "grep bar | foo")) {
    } catch (RuntimeException e) {
      pass = true;
      assertTrue(e.getMessage().contains("foo: command not found"));
    }
    assertTrue(pass);
  }

  @Test
  public void canFilterUsingAwk() throws IOException {

    String[] in3 = {"chr1\t1\t100", "chr1\t10\t100", "chr1\t2\t100"};
    SystemCommand sysCmd = new SystemCommand();
    try (Stream<String> s =
        sysCmd.streamLinesThroughAwk(Stream.of(in3), "{VAR=5} $2 > VAR && $1")) {
      List<String> results = s.toList();
      assertEquals("chr1\t10\t100", results.get(0));
      assertEquals(1, results.size());
    }

    String[] in = {"chr1\t1\t100", "chr1\t1\t100", "chr1\t2\t100"};
    try (Stream<String> s = sysCmd.streamLinesThroughAwk(Stream.of(in), "{VAR=5} $2 > VAR && $1")) {
      List<String> results = s.toList();
      assertEquals(0, results.size());
    }

    String[] in4 = {"chr1\t10\t100", "chr1\t10\t100", "chr1\t2\t100"};
    try (Stream<String> s =
        sysCmd.streamLinesThroughAwk(Stream.of(in4), "{VAR=5} $2 > VAR && $1")) {
      List<String> results = s.toList();
      assertEquals(2, results.size());
      assertEquals(in4[0], results.get(0));
      assertEquals(in4[1], results.get(1));
    }

    String[] in2 = {};
    try (Stream<String> s =
        sysCmd.streamLinesThroughAwk(Stream.of(in2), "{VAR=5} $2 > VAR && $1")) {
      List<String> results = s.toList();
      assertEquals(0, results.size());
    }
    sysCmd.deleteTempFile();
  }

  @Test
  public void canHandleQuotesInAwkScript() throws IOException {
    // NB: A single \ needs to be \\ in Java strings
    String[] in = {"chr'1"};
    SystemCommand sysCmd = new SystemCommand();
    try (Stream<String> s = sysCmd.streamLinesThroughAwk(Stream.of(in), "$1 == \"chr'1\"")) {
      List<String> results = s.toList();
      assertEquals(1, results.size());
    }

    // Double quote. Note three '\': Two to represent a single \ and one to escape
    // the double quote in Java string
    String[] in2 = {"chr\"1"};
    try (Stream<String> s = sysCmd.streamLinesThroughAwk(Stream.of(in2), "$1 == \"chr\\\"1\"")) {
      List<String> results = s.toList();
      assertEquals(1, results.size());
    }

    String[] in3 = {"chr\\1"};
    try (Stream<String> s = sysCmd.streamLinesThroughAwk(Stream.of(in3), "$1 == \"chr\\\\1\"")) {
      List<String> results = s.toList();
      assertEquals(1, results.size());
    }
  }

  @Test
  public void canDetectBrokenAwkScript() throws InterruptedException, IOException {
    SystemCommand sysCmd = new SystemCommand();
    String[] in = {"'chr1\t10\t100'"};
    for (int i = 0; i < 3; i++) {
      // Broken awk script:
      boolean pass = false;
      try (Stream<String> s = sysCmd.streamLinesThroughAwk(Stream.of(in), "'print {'")) {

      } catch (RuntimeException e) {
        pass = true;
      }
      assertTrue(pass);
      Thread.sleep(1000);
    }
  }

  @Test
  public void canFilterGtfAwkFunc() throws IOException {
    // NB: You need to set the field separator to \\t.
    SystemCommand sysCmd = new SystemCommand();
    String[] gtf = {
      "GL873520\tchr1\tstop_codon\t8064\t8066\t0.000000\t-\t.\tgene_id 100; trax_id \"ACTB\";"
    };
    try (Stream<String> s =
        sysCmd.streamLinesThroughAwk(Stream.of(gtf), "getGtfTag(\"gene_id\") == 100")) {
      List<String> results = s.toList();
      assertEquals(gtf[0], results.get(0));
    }
    try (Stream<String> s =
        sysCmd.streamLinesThroughAwk(Stream.of(gtf), "getGtfTag(\"trax_id\") == \"ACTB\"")) {
      List<String> results = s.toList();
      assertEquals(gtf[0], results.get(0));
    }
    // Empty string if key not found
    try (Stream<String> s =
        sysCmd.streamLinesThroughAwk(Stream.of(gtf), "getGtfTag(\"SPAM\") == \"\"")) {
      List<String> results = s.toList();
      assertEquals(gtf[0], results.get(0));
    }

    // No attributes at all: Empty string returned
    gtf = new String[] {"GL873520\tchr1\tstop_codon\t8064\t8066\t0.000000\t-\t.\t."};
    try (Stream<String> s =
        sysCmd.streamLinesThroughAwk(Stream.of(gtf), "getGtfTag(\"gene_id\") == \"\"")) {
      List<String> results = s.toList();
      assertEquals(gtf[0], results.get(0));
    }

    // No attribute column at all (this would be an invalid GTF, by the way)
    gtf = new String[] {"GL873520\tchr1\tstop_codon\t8064\t8066\t0.000000\t-\t."};
    try (Stream<String> s =
        sysCmd.streamLinesThroughAwk(Stream.of(gtf), "getGtfTag(\"gene_id\") == \"\"")) {
      List<String> results = s.toList();
      assertEquals(gtf[0], results.get(0));
    }
  }

  @Test
  public void canFilterGffAwkFunc() throws IOException {
    SystemCommand sysCmd = new SystemCommand();
    // NB: You need to set the field separator to \\t.
    String[] x = {
      ".|.|.|.|.|.|.|.|Tag=100; ID = foo : bar ; Alias=spam,bar;".replaceAll("\\|", "\t")
    };
    try (Stream<String> s =
        sysCmd.streamLinesThroughAwk(Stream.of(x), "getGffTag(\"Tag\") == 100")) {
      List<String> results = s.toList();
      assertEquals(x[0], results.get(0));
    }
    try (Stream<String> s =
        sysCmd.streamLinesThroughAwk(Stream.of(x), "getGffTag(\"ID\") == \"foo : bar\"")) {
      List<String> results = s.toList();
      assertEquals(x[0], results.get(0));
    }
    try (Stream<String> s =
        sysCmd.streamLinesThroughAwk(Stream.of(x), "getGffTag(\"Alias\") == \"spam,bar\"")) {
      List<String> results = s.toList();
      assertEquals(x[0], results.get(0));
    }
    try (Stream<String> s =
        sysCmd.streamLinesThroughAwk(Stream.of(x), "getGffTag(\"Alias\") == \"spam,bar\"")) {
      List<String> results = s.toList();
      assertEquals(x[0], results.get(0));
    }
    try (Stream<String> s =
        sysCmd.streamLinesThroughAwk(Stream.of(x), "getGffTag(\"Alias\", 1) == \"spam\"")) {
      List<String> results = s.toList();
      assertEquals(x[0], results.get(0));
    }
    try (Stream<String> s =
        sysCmd.streamLinesThroughAwk(Stream.of(x), "getGffTag(\"Alias\", 2) == \"bar\"")) {
      List<String> results = s.toList();
      assertEquals(x[0], results.get(0));
    }
    try (Stream<String> s =
        sysCmd.streamLinesThroughAwk(Stream.of(x), "getGffTag(\"Alias\", 99) == \"\"")) {
      List<String> results = s.toList();
      assertEquals(x[0], results.get(0));
    }
    try (Stream<String> s =
        sysCmd.streamLinesThroughAwk(Stream.of(x), "getGffTag(\"SPAM\") == \"\"")) {
      List<String> results = s.toList();
      assertEquals(x[0], results.get(0));
    }

    // NB: Missing tag i.e., empty string, evaluates to 0!!
    try (Stream<String> s =
        sysCmd.streamLinesThroughAwk(Stream.of(x), "getGffTag(\"SPAM\") == 0")) {
      List<String> results = s.toList();
      assertEquals(x[0], results.get(0));
    }

    x = new String[] {".|.|.|.|.|.|.|.|.".replaceAll("\\|", "\t")};
    try (Stream<String> s =
        sysCmd.streamLinesThroughAwk(Stream.of(x), "getGffTag(\"Alias\") == \"\"")) {
      List<String> results = s.toList();
      assertEquals(x[0], results.get(0));
    }

    x = new String[] {".|.|.|.|.|.|.|.".replaceAll("\\|", "\t")};
    try (Stream<String> s =
        sysCmd.streamLinesThroughAwk(Stream.of(x), "getGffTag(\"Alias\") == \"\"")) {
      List<String> results = s.toList();
      assertEquals(x[0], results.get(0));
    }

    // The value of the gff tag includes double-quotes. So it is "X" not X
    x =
        new String[] {
          ".|.|.|.|.|.|.|.|Tag=\"X\"".replaceAll("\\|", "\t")
        }; // Double quotes are not stripped
    try (Stream<String> s =
        sysCmd.streamLinesThroughAwk(Stream.of(x), "getGffTag(\"Tag\") == \"\\\"X\\\"\"")) {
      List<String> results = s.toList();
      assertEquals(x[0], results.get(0));
    }
  }

  @Test
  public void canFilterInfoVcfWithAwkFunc() throws IOException {
    SystemCommand sysCmd = new SystemCommand();
    String[] vcf = {
      "chr1 75888 . A T . . IMPRECISE;SVTYPE=DEL;DP=20,30;SVLEN=-32945;FOLD_CHANGE=0.723657;FOLD_CHANGE_LOG=-0.466623;PROBES=21"
          .replaceAll(" ", "\t")
    };
    try (Stream<String> s =
        sysCmd.streamLinesThroughAwk(Stream.of(vcf), "getInfoTag(\"IMPRECISE\") == 1")) {
      List<String> results = s.toList();
      assertEquals(vcf[0], results.get(0));
    }
    try (Stream<String> s =
        sysCmd.streamLinesThroughAwk(Stream.of(vcf), "getInfoTag(\"ABSENT_TAG\") == 0")) {
      List<String> results = s.toList();
      assertEquals(vcf[0], results.get(0));
    }
    try (Stream<String> s =
        sysCmd.streamLinesThroughAwk(Stream.of(vcf), "getInfoTag(\"FOLD_CHANGE\") <= 0.723657")) {
      List<String> results = s.toList();
      assertEquals(vcf[0], results.get(0));
    }
    try (Stream<String> s =
        sysCmd.streamLinesThroughAwk(Stream.of(vcf), "getInfoTag(\"FOLD_CHANGE\") > 0.723657")) {
      List<String> results = s.toList();
      assertEquals(0, results.size());
    }

    // Split list of values
    try (Stream<String> s =
        sysCmd.streamLinesThroughAwk(Stream.of(vcf), "getInfoTag(\"DP\") == \"20,30\"")) {
      List<String> results = s.toList();
      assertEquals(vcf[0], results.get(0));
    }
    try (Stream<String> s =
        sysCmd.streamLinesThroughAwk(Stream.of(vcf), "getInfoTag(\"DP\", 1) == 20")) {
      List<String> results = s.toList();
      assertEquals(vcf[0], results.get(0));
    }
    try (Stream<String> s =
        sysCmd.streamLinesThroughAwk(Stream.of(vcf), "getInfoTag(\"DP\", 2) == 30")) {
      List<String> results = s.toList();
      assertEquals(vcf[0], results.get(0));
    }

    // Out of range
    try (Stream<String> s =
        sysCmd.streamLinesThroughAwk(Stream.of(vcf), "getInfoTag(\"DP\", 3) == \"\"")) {
      List<String> results = s.toList();
      assertEquals(vcf[0], results.get(0));
    }

    // Missing INFO
    vcf = new String[] {"chr1 75888 . A T . . .".replaceAll(" ", "\t")};
    try (Stream<String> s =
        sysCmd.streamLinesThroughAwk(Stream.of(vcf), "getInfoTag(\"FOO\") == 0")) {
      List<String> results = s.toList();
      assertEquals(vcf[0], results.get(0));
    }

    // INFO column not there at all
    vcf = new String[] {"chr1 75888 . A T . .".replaceAll(" ", "\t")};
    try (Stream<String> s =
        sysCmd.streamLinesThroughAwk(Stream.of(vcf), "getInfoTag(\"FOO\") == 0")) {
      List<String> results = s.toList();
      assertEquals(vcf[0], results.get(0));
    }
  }

  @Test
  public void canFilterFormatVcfWithAwkFunc() throws IOException {
    SystemCommand sysCmd = new SystemCommand();
    String[] vcf = {"chr1 75888 . A T . . . GT:GQ 11:21,10 22:99,100".replaceAll(" ", "\t")};
    try (Stream<String> s =
        sysCmd.streamLinesThroughAwk(Stream.of(vcf), "getFmtTag(\"GT\") == 11")) {
      List<String> results = s.toList();
      assertEquals(vcf[0], results.get(0)); // Default sample_idx= 1
    }
    try (Stream<String> s =
        sysCmd.streamLinesThroughAwk(Stream.of(vcf), "getFmtTag(\"GT\", 1) == 11")) {
      List<String> results = s.toList();
      assertEquals(vcf[0], results.get(0));
    }
    try (Stream<String> s =
        sysCmd.streamLinesThroughAwk(Stream.of(vcf), "getFmtTag(\"GT\", 2) == 22")) {
      List<String> results = s.toList();
      assertEquals(vcf[0], results.get(0));
    }
    try (Stream<String> s =
        sysCmd.streamLinesThroughAwk(Stream.of(vcf), "getFmtTag(\"GT\", 2) < 22")) {
      List<String> results = s.toList();
      assertEquals(0, results.size());
    }

    // Get value from list
    try (Stream<String> s =
        sysCmd.streamLinesThroughAwk(Stream.of(vcf), "getFmtTag(\"GQ\", 2) == \"99,100\"")) {
      List<String> results = s.toList();
      assertEquals(vcf[0], results.get(0));
    }
    try (Stream<String> s =
        sysCmd.streamLinesThroughAwk(Stream.of(vcf), "getFmtTag(\"GQ\", 2, 1) == \"99\"")) {
      List<String> results = s.toList();
      assertEquals(vcf[0], results.get(0));
    }
    try (Stream<String> s =
        sysCmd.streamLinesThroughAwk(Stream.of(vcf), "getFmtTag(\"GQ\", 2, 2) == \"100\"")) {
      List<String> results = s.toList();
      assertEquals(vcf[0], results.get(0));
    }
    // Out range
    try (Stream<String> s =
        sysCmd.streamLinesThroughAwk(Stream.of(vcf), "getFmtTag(\"GQ\", 2, 3) == \"\"")) {
      List<String> results = s.toList();
      assertEquals(vcf[0], results.get(0));
    }

    // Tag not found
    try (Stream<String> s =
        sysCmd.streamLinesThroughAwk(Stream.of(vcf), "getFmtTag(\"ABSENT\", 1) == \"\"")) {
      List<String> results = s.toList();
      assertEquals(vcf[0], results.get(0));
    }

    // Invalid indexes
    try (Stream<String> s =
        sysCmd.streamLinesThroughAwk(Stream.of(vcf), "getFmtTag(\"GT\", 99) == \"\"")) {
      List<String> results = s.toList();
      assertEquals(vcf[0], results.get(0));
    }
    try (Stream<String> s =
        sysCmd.streamLinesThroughAwk(Stream.of(vcf), "getFmtTag(\"GT\", \"foobar\") == \"\"")) {
      List<String> results = s.toList();
      assertEquals(vcf[0], results.get(0));
    }
  }

  @Test
  public void canFilterSamTagWithAwk() throws IOException {
    SystemCommand sysCmd = new SystemCommand();
    String[] rec = {
      "read\t0\tchr7\t5566778\t50\t5M\t*\t0\t0\tCTCAT\tIIIII\tMD:Z:75\tRG:Z:1\tXG:i:0\tNH:i:1"
          + "\tNM:i:0\tXM:i:0\tXN:i:0\tXO:i:0\tAS:i:0\tYT:Z:UU"
    };

    // Filter for NH tag value
    try (Stream<String> s = sysCmd.streamLinesThroughAwk(Stream.of(rec), "getSamTag(\"NH\") > 0")) {
      List<String> results = s.toList();
      assertEquals(rec[0], results.get(0));
    }
    try (Stream<String> s =
        sysCmd.streamLinesThroughAwk(Stream.of(rec), "getSamTag(\"NH\") > 10")) {
      List<String> results = s.toList();
      assertEquals(0, results.size());
    }

    // Missing tag
    try (Stream<String> s = sysCmd.streamLinesThroughAwk(Stream.of(rec), "getSamTag(\"ZZ\") > 0")) {
      List<String> results = s.toList();
      assertEquals(0, results.size());
    }

    // Missing tag searched but not used
    try (Stream<String> s =
        sysCmd.streamLinesThroughAwk(Stream.of(rec), "{getSamTag(\"ZZ\"); print $0}")) {
      List<String> results = s.toList();
      assertEquals(rec[0], results.get(0));
    }
    try (Stream<String> s = sysCmd.streamLinesThroughAwk(Stream.of(rec), "getSamTag(\"NM\") > 0")) {
      List<String> results = s.toList();
      assertEquals(0, results.size());
    }

    // Tags missing altogether returns empty string
    rec = new String[] {"read\t0\tchr7\t5566778\t50\t5M\t*\t0\t0\tCTCAT\tIIIII"};
    try (Stream<String> s =
        sysCmd.streamLinesThroughAwk(Stream.of(rec), "getSamTag(\"NM\") == \"\"")) {
      List<String> results = s.toList();
      assertEquals(rec[0], results.get(0));
    }

    long t0 = System.currentTimeMillis();
    int i = 0;
    while (i < 1000) {
      try (Stream<String> s =
          sysCmd.streamLinesThroughAwk(Stream.of(rec), "getSamTag(\"NH\") > 0")) {
        List<String> results = s.toList();
        assertEquals(0, results.size());
      }
      i++;
    }
    long t1 = System.currentTimeMillis();
    assertTrue((t1 - t0) < 10000); // It can filter reasonably fast (?)
  }

  @Test
  public void canFilterSamBitflagWithAwk() throws IOException {
    SystemCommand sysCmd = new SystemCommand();

    String[] rec = new String[] {"read\t16"};
    try (Stream<String> s = sysCmd.streamLinesThroughAwk(Stream.of(rec), "bitset(16) == 1")) {
      List<String> results = s.toList();
      assertEquals(rec[0], results.get(0));
    }

    rec = new String[] {"read\t17"};
    try (Stream<String> s = sysCmd.streamLinesThroughAwk(Stream.of(rec), "bitset(16) == 1")) {
      List<String> results = s.toList();
      assertEquals(rec[0], results.get(0));
    }

    rec = new String[] {"read\t1"};
    try (Stream<String> s = sysCmd.streamLinesThroughAwk(Stream.of(rec), "bitset(16) == 0")) {
      List<String> results = s.toList();
      assertEquals(rec[0], results.get(0));
    }

    rec = new String[] {"read\t3585"};
    try (Stream<String> s = sysCmd.streamLinesThroughAwk(Stream.of(rec), "bitset(1025) == 1")) {
      List<String> results = s.toList();
      assertEquals(rec[0], results.get(0));
    }

    rec = new String[] {"read\t3584"};
    try (Stream<String> s = sysCmd.streamLinesThroughAwk(Stream.of(rec), "bitset(1025) == 0")) {
      List<String> results = s.toList();
      assertEquals(rec[0], results.get(0));
    }

    rec = new String[] {"read\t3840"};
    try (Stream<String> s = sysCmd.streamLinesThroughAwk(Stream.of(rec), "bitset(3840) == 1")) {
      List<String> results = s.toList();
      assertEquals(rec[0], results.get(0));
    }

    rec = new String[] {"read\t3840"};
    try (Stream<String> s = sysCmd.streamLinesThroughAwk(Stream.of(rec), "bitset(3841) == 0")) {
      List<String> results = s.toList();
      assertEquals(rec[0], results.get(0));
    }

    rec = new String[] {"read\t3840"};
    try (Stream<String> s = sysCmd.streamLinesThroughAwk(Stream.of(rec), "bitset(0) == 1")) {
      List<String> results = s.toList();
      assertEquals(rec[0], results.get(0));
    }

    rec = new String[] {"read\t0"};
    try (Stream<String> s = sysCmd.streamLinesThroughAwk(Stream.of(rec), "bitset(0) == 1")) {
      List<String> results = s.toList();
      assertEquals(rec[0], results.get(0));
    }

    // Invalid input: Always return 0 (false)
    rec = new String[] {"read\tfoo"};
    try (Stream<String> s = sysCmd.streamLinesThroughAwk(Stream.of(rec), "bitset(0) == 0")) {
      List<String> results = s.toList();
      assertEquals(rec[0], results.get(0));
    }

    rec = new String[] {"read\t16"};
    try (Stream<String> s = sysCmd.streamLinesThroughAwk(Stream.of(rec), "bitset(\"foo\") == 0")) {
      List<String> results = s.toList();
      assertEquals(rec[0], results.get(0));
    }

    rec = new String[] {"read\t16"};
    try (Stream<String> s = sysCmd.streamLinesThroughAwk(Stream.of(rec), "bitset(-16) == 0")) {
      List<String> results = s.toList();
      assertEquals(rec[0], results.get(0));
    }

    rec = new String[] {"read\t-16"};
    try (Stream<String> s = sysCmd.streamLinesThroughAwk(Stream.of(rec), "bitset(0) == 0")) {
      List<String> results = s.toList();
      assertEquals(rec[0], results.get(0));
    }

    rec = new String[] {"read\t16"};
    try (Stream<String> s = sysCmd.streamLinesThroughAwk(Stream.of(rec), "bitset(17) == 0")) {
      List<String> results = s.toList();
      assertEquals(rec[0], results.get(0));
    }
  }
}
