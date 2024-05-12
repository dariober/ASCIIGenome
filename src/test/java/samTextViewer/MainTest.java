package samTextViewer;

import static org.junit.Assert.*;

import com.google.common.base.Joiner;
import com.itextpdf.text.DocumentException;
import exceptions.BamIndexNotFoundException;
import exceptions.InvalidColourException;
import exceptions.InvalidCommandLineException;
import exceptions.InvalidConfigException;
import exceptions.InvalidGenomicCoordsException;
import exceptions.InvalidRecordException;
import faidx.UnindexableFastaFileException;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.List;
import org.junit.Test;

public class MainTest {

  @Test
  public void canStartFromCram()
      throws ClassNotFoundException,
          IOException,
          InvalidGenomicCoordsException,
          InvalidCommandLineException,
          InvalidRecordException,
          BamIndexNotFoundException,
          SQLException,
          DocumentException,
          UnindexableFastaFileException,
          InvalidColourException,
          InvalidConfigException {
    String[] args =
        new String[] {
          "-ni",
          "-nf",
          "--debug",
          "2",
          "-fa",
          "test_data/chr7.fa",
          "--exec",
          "goto chr7:5567419-5567599",
          "test_data/ds051.actb.cram"
        };
    String out = Joiner.on("\n").join(this.runMain(args));
    assertTrue(out.contains("chr7:5567419-5567599") && out.contains("<<<<<"));
  }

  @Test
  /*You should really test this in InteractiveInputTest.java but setting it up is a bit of a mess */
  public void canGoToNextChromosome()
      throws ClassNotFoundException,
          IOException,
          InvalidGenomicCoordsException,
          InvalidCommandLineException,
          InvalidRecordException,
          BamIndexNotFoundException,
          SQLException,
          DocumentException,
          UnindexableFastaFileException,
          InvalidColourException,
          InvalidConfigException {
    String[] args =
        new String[] {"-ni", "-nf", "--exec", "nextChrom -s u", "test_data/ds051.actb.bam"};
    String out = Joiner.on("\n").join(this.runMain(args));
    System.err.println(out);
    assertTrue(out.contains("chr8:1-"));

    args =
        new String[] {
          "-ni", "-nf", "--exec", "goto chrM && nextChrom -s u", "test_data/ds051.actb.bam"
        };
    out = Joiner.on("\n").join(this.runMain(args));
    assertTrue(out.contains("chr1:1-"));

    // Wrap around
    args =
        new String[] {
          "-ni", "-nf", "--exec", "goto chrY && nextChrom -s u", "test_data/ds051.actb.bam"
        };
    out = Joiner.on("\n").join(this.runMain(args));
    assertTrue(out.contains("chrM:1-"));

    // One chrom - stay there.
    args =
        new String[] {"-ni", "-nf", "--exec", "nextChrom -s u", "test_data/refSeq.hg19.short.bed"};
    out = Joiner.on("\n").join(this.runMain(args));
    assertTrue(out.contains("chr1:67208779"));

    // No tracks, no genome
    args = new String[] {"-ni", "-nf", "--exec", "nextChrom"};
    out = Joiner.on("\n").join(this.runMain(args));
    assertTrue(out.contains("Undefined_contig:1-"));

    // Only genome
    args = new String[] {"-ni", "-nf", "-fa", "test_data/seq_cg.fa", "--exec", "nextChrom"};
    out = Joiner.on("\n").join(this.runMain(args));
    assertTrue(out.contains("seq:1-"));
  }

  @Test
  /*You should really test this in InteractiveInputTest.java but setting it up is a bit of a mess */
  public void canGoToNextChromosomeRegex()
      throws ClassNotFoundException,
          IOException,
          InvalidGenomicCoordsException,
          InvalidCommandLineException,
          InvalidRecordException,
          BamIndexNotFoundException,
          SQLException,
          DocumentException,
          UnindexableFastaFileException,
          InvalidColourException,
          InvalidConfigException {
    String[] args =
        new String[] {"-ni", "-nf", "--exec", "nextChrom M", "test_data/ds051.actb.bam"};
    String out = Joiner.on("\n").join(this.runMain(args));
    assertTrue(out.contains("chrM:1-"));

    args =
        new String[] {
          "-ni", "-nf", "--exec", "nextChrom -min 249000000 chr1", "test_data/ds051.actb.bam"
        };
    out = Joiner.on("\n").join(this.runMain(args));
    assertTrue(out.contains("chr1:1-"));

    args =
        new String[] {
          "-ni", "-nf", "--exec", "nextChrom -min 135000000 chr1", "test_data/ds051.actb.bam"
        };
    out = Joiner.on("\n").join(this.runMain(args));
    System.err.println(out);
    assertTrue(out.contains("chr11:1-"));
  }

  @Test
  public void canSuggestCommand()
      throws ClassNotFoundException,
          IOException,
          InvalidGenomicCoordsException,
          InvalidCommandLineException,
          InvalidRecordException,
          BamIndexNotFoundException,
          SQLException,
          DocumentException,
          UnindexableFastaFileException,
          InvalidColourException,
          InvalidConfigException {
    String[] args = new String[] {"-ni", "-nf", "--exec", "prnt"};
    List<String> out = this.runMain(args);
    assertTrue(out.get(1).contains("Maybe you mean print?"));
  }

  @Test
  public void canSetConfig()
      throws ClassNotFoundException,
          IOException,
          InvalidGenomicCoordsException,
          InvalidCommandLineException,
          InvalidRecordException,
          BamIndexNotFoundException,
          SQLException,
          DocumentException,
          UnindexableFastaFileException,
          InvalidColourException,
          InvalidConfigException {
    String[] args =
        new String[] {"-ni", "-nf", "--exec", "setConfig nucs f", "test_data/ds051.short.bam"};
    List<String> out = this.runMain(args);
    assertTrue(out.get(0).contains(">>>>>>>>>>>>>>>>>>>"));
  }

  @Test
  public void canFlipBooleanConfig()
      throws ClassNotFoundException,
          IOException,
          InvalidGenomicCoordsException,
          InvalidCommandLineException,
          InvalidRecordException,
          BamIndexNotFoundException,
          SQLException,
          DocumentException,
          UnindexableFastaFileException,
          InvalidColourException,
          InvalidConfigException {
    String[] args =
        new String[] {"-ni", "-nf", "--exec", "setConfig nucs", "test_data/ds051.short.bam"};
    List<String> out = this.runMain(args);
    assertTrue(out.get(0).contains(">>>>>>>>>>>>>>>>>>>"));
  }

  @Test
  public void doNotSetInvalidBooleanConfig()
      throws ClassNotFoundException,
          IOException,
          InvalidGenomicCoordsException,
          InvalidCommandLineException,
          InvalidRecordException,
          BamIndexNotFoundException,
          SQLException,
          DocumentException,
          UnindexableFastaFileException,
          InvalidColourException,
          InvalidConfigException {
    String[] args =
        new String[] {"-ni", "-nf", "--exec", "setConfig nucs 999", "test_data/ds051.short.bam"};
    List<String> out = this.runMain(args);
    System.err.println(out);
    assertTrue(out.get(1).contains("Unable to set"));
  }

  @Test
  public void doNotSetInvalidColourConfig()
      throws ClassNotFoundException,
          IOException,
          InvalidGenomicCoordsException,
          InvalidCommandLineException,
          InvalidRecordException,
          BamIndexNotFoundException,
          SQLException,
          DocumentException,
          UnindexableFastaFileException,
          InvalidColourException,
          InvalidConfigException {
    String[] args =
        new String[] {"-ni", "-nf", "--exec", "setConfig seq_a 999", "test_data/ds051.short.bam"};
    List<String> out = this.runMain(args);
    System.err.println(out);
    assertTrue(out.get(1).contains("Unable to set"));
  }

  @Test
  public void doNotSetInvalidIntegerConfig()
      throws ClassNotFoundException,
          IOException,
          InvalidGenomicCoordsException,
          InvalidCommandLineException,
          InvalidRecordException,
          BamIndexNotFoundException,
          SQLException,
          DocumentException,
          UnindexableFastaFileException,
          InvalidColourException,
          InvalidConfigException {
    String[] args =
        new String[] {
          "-ni", "-nf", "--exec", "setConfig shade_baseq foo", "test_data/ds051.short.bam"
        };
    List<String> out = this.runMain(args);
    System.err.println(out);
    assertTrue(out.get(1).contains("Unable to set"));
  }

  @Test
  public void canIgnoreComments()
      throws ClassNotFoundException,
          IOException,
          InvalidGenomicCoordsException,
          InvalidCommandLineException,
          InvalidRecordException,
          BamIndexNotFoundException,
          SQLException,
          DocumentException,
          UnindexableFastaFileException,
          InvalidColourException,
          InvalidConfigException {
    String[] args =
        new String[] {
          "-ni", "-nf", "--exec", "print && grep -i NCTNTCCN", "test_data/ds051.short.bam"
        };
    List<String> woComm = this.runMain(args);
    args =
        new String[] {
          "-ni",
          "-nf",
          "--exec",
          "print && grep -i NCTNTCCN // A comment",
          "test_data/ds051.short.bam"
        };
    List<String> withComm = this.runMain(args);
    assertEquals(woComm, withComm);

    args = new String[] {"-ni", "-nf", "--exec", "goto chr7", "test_data/ds051.short.bam"};
    woComm = this.runMain(args);

    args = new String[] {"-ni", "-nf", "--exec", "goto chr7//comment", "test_data/ds051.short.bam"};
    withComm = this.runMain(args);

    assertEquals(woComm, withComm);

    args =
        new String[] {"-ni", "-nf", "--exec", "goto chr7 // comment", "test_data/ds051.short.bam"};
    withComm = this.runMain(args);

    assertEquals(woComm, withComm);
  }

  /* H E L P E R S */

  /**
   * Execute main with the given array of arguments and return a list of length 2 containing 1)
   * stdout and 2) stderr.
   */
  private List<String> runMain(String[] args)
      throws ClassNotFoundException,
          IOException,
          InvalidGenomicCoordsException,
          InvalidCommandLineException,
          InvalidRecordException,
          BamIndexNotFoundException,
          SQLException,
          DocumentException,
          UnindexableFastaFileException,
          InvalidColourException,
          InvalidConfigException {

    PrintStream stdout = System.out;
    ByteArrayOutputStream baosOut = new ByteArrayOutputStream();
    System.setOut(new PrintStream(baosOut));

    PrintStream stderr = System.err;
    ByteArrayOutputStream baosErr = new ByteArrayOutputStream();
    System.setErr(new PrintStream(baosErr));

    Main.main(args);

    String out = baosOut.toString();
    System.setOut(stdout);

    String err = baosErr.toString();
    System.setErr(stderr);

    List<String> outErr = new ArrayList<String>();
    outErr.add(out);
    outErr.add(err);
    return outErr;
  }
}
