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
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import org.apache.commons.lang3.StringUtils;
import org.junit.Test;

public class MainTest {

  @Test
  public void canFilterGenotype()
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
          "--exec",
          "genotype -f '{GT} == \"1|2\"'",
          "test_data/info_formats.vcf.gz"
        };
    // Genotype filter selects for "sample2". "sample1" appears once right at the start before the
    // filtering is applied.
    String out = Joiner.on("\n").join(this.runMain(args));
    assertEquals(1, StringUtils.countMatches(out, "sample1"));
    assertEquals(2, StringUtils.countMatches(out, "sample2"));
  }

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
          "test_data/chr7.fa",
          "--exec",
          "goto chr7:5567419-5567599",
          "test_data/ds051.actb.cram"
        };
    String out = Joiner.on("\n").join(this.runMain(args));
    assertTrue(out.contains("chr7:5567419-5567599") && out.contains("<<<<<"));
  }

  @Test
  public void canStartFromCramVersion3_1()
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
          "test_data/chr7.fa",
          "--exec",
          "goto chr7:5567419-5567599",
          "test_data/ds051.actb.v3_1.cram"
        };
    String out = Joiner.on("\n").join(this.runMain(args));
    assertTrue(out.contains("chr7:5567419-5567599") && out.contains("<<<<<"));
  }

  @Test
  public void canHandlePairedFilterInSingleEndBam()
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
          "-r",
          "chr7:5566778-5566857",
          "-ni",
          "-nf",
          "--debug",
          "2",
          "--exec",
          "samtools -F 2",
          "test_data/ds051.actb.bam"
        };
    String out = Joiner.on("\n").join(this.runMain(args));
    assertTrue(out.contains("Reads: 292"));

    args =
        new String[] {
          "-r",
          "chr7:5566778-5566857",
          "-ni",
          "-nf",
          "--debug",
          "2",
          "--exec",
          "samtools -f 2",
          "test_data/ds051.actb.bam"
        };
    out = Joiner.on("\n").join(this.runMain(args));
    assertTrue(out.contains("Reads: 0"));
  }

  @Test
  /*Add these to canGoToNextChromosome in InteractiveInputTest.java */
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
    // One chrom - stay there.
    String[] args =
        new String[] {"-ni", "-nf", "--exec", "nextChrom -s u", "test_data/refSeq.hg19.short.bed"};
    String out = Joiner.on("\n").join(this.runMain(args));
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
    assertTrue(out.contains("chr11:1-"));
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
  public void canInitWithoutDictionary()
      throws UnindexableFastaFileException,
          SQLException,
          DocumentException,
          InvalidGenomicCoordsException,
          InvalidCommandLineException,
          InvalidColourException,
          IOException,
          InvalidConfigException,
          BamIndexNotFoundException,
          ClassNotFoundException,
          InvalidRecordException {
    String[] args =
        new String[] {"-ni", "-nf", "test_data/Homo_sapiens.GRCh38.86.ENST00000331789.gff3"};
    List<String> out = this.runMain(args);
    assertTrue(out.get(0).contains("ENST00000331789"));
  }

  @Test
  public void addHeaderIssue108()
      throws UnindexableFastaFileException,
          SQLException,
          DocumentException,
          InvalidGenomicCoordsException,
          InvalidCommandLineException,
          InvalidColourException,
          IOException,
          InvalidConfigException,
          BamIndexNotFoundException,
          ClassNotFoundException,
          InvalidRecordException {
    String[] args =
        new String[] {
          "-ni",
          "-nf",
          "--exec",
          "addHeader FOOBAR @2 && addHeader 'SPAMEGGS {-}'",
          "test_data/ds051.actb.bam"
        };
    String out = this.runMain(args).get(0);

    // Add SPAMEGGS only to the track(s) with header
    Matcher matcher = Pattern.compile(Pattern.quote("SPAMEGGS")).matcher(out);
    assertEquals(1, matcher.results().count());

    matcher = Pattern.compile(Pattern.quote("SPAMEGGS FOOBAR")).matcher(out);
    assertEquals(1, matcher.results().count());

    matcher = Pattern.compile(Pattern.quote("{-}")).matcher(out);
    assertEquals(0, matcher.results().count());
  }

  @Test
  public void showGenomeIssue107()
      throws UnindexableFastaFileException,
          SQLException,
          DocumentException,
          InvalidGenomicCoordsException,
          InvalidCommandLineException,
          InvalidColourException,
          IOException,
          InvalidConfigException,
          BamIndexNotFoundException,
          ClassNotFoundException,
          InvalidRecordException {
    String[] args = new String[] {"-ni", "-nf", "--exec", "show ge", "test_data/chr7.fa"};
    String out = this.runMain(args).toString();
    assertTrue(out.contains("chr7 159138663 |||||"));

    args =
        new String[] {"-ni", "-nf", "--exec", "seqRegex ACTCACTG && show ge", "test_data/chr7.fa"};
    out = this.runMain(args).toString();
    assertTrue(out.contains("chr7 159138663 |||||"));

    args = new String[] {"-ni", "-nf", "--exec", "show ge", "test_data/ds051.actb.bam"};
    out = this.runMain(args).toString();
    assertTrue(out.contains("chr1  249250621 |||||||||"));

    args =
        new String[] {
          "-ni", "-nf", "--exec", "show ge", "test_data/Homo_sapiens.GRCh38.86.ENST00000331789.gff3"
        };
    out = this.runMain(args).toString();
    assertTrue(out.contains("Number of contigs: 1"));
  }

  @Test
  public void canParseLinesWithBash()
      throws UnindexableFastaFileException,
          SQLException,
          DocumentException,
          InvalidGenomicCoordsException,
          InvalidCommandLineException,
          InvalidColourException,
          IOException,
          InvalidConfigException,
          BamIndexNotFoundException,
          ClassNotFoundException,
          InvalidRecordException {
    String[] args =
        new String[] {
          "-ni", "-nf", "--exec", "print -sys 'cut -f 1,3'", "test_data/hg19_genes_head.gtf.gz"
        };
    String out = this.runMain(args).toString();
    assertTrue(out.contains("chr1 | exon"));

    args =
        new String[] {
          "-ni",
          "-nf",
          "--exec",
          "print -sys \"cut -f 1,3 | sed 's/exon/EXON/'\"",
          "test_data/hg19_genes_head.gtf.gz"
        };
    out = this.runMain(args).toString();
    assertTrue(out.contains("chr1 | EXON"));

    args =
        new String[] {
          "-ni",
          "-nf",
          "--exec",
          "print -sys \"cut -f 1,3 | sed 's/exon//'\"",
          "test_data/hg19_genes_head.gtf.gz"
        };
    out = this.runMain(args).toString();
    assertFalse(out.contains("exon"));
  }

  @Test
  public void seqRegexIssue123()
      throws UnindexableFastaFileException,
          SQLException,
          DocumentException,
          InvalidGenomicCoordsException,
          InvalidCommandLineException,
          InvalidColourException,
          IOException,
          InvalidConfigException,
          BamIndexNotFoundException,
          ClassNotFoundException,
          InvalidRecordException {
    String[] args = new String[] {"-ni", "-nf", "--exec", "seqRegex ACTG", "test_data/seq_cg.fa"};
    String out = this.runMain(args).toString();
    assertTrue(out.contains("ACTGACTG"));

    args =
        new String[] {"-ni", "-nf", "--exec", "seqRegex ACTG && nextChrom", "test_data/seq_cg.fa"};
    out = this.runMain(args).toString();
    assertTrue(out.contains("ATCGATCGATCGATCG"));
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
