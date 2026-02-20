package samTextViewer;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import colouring.Config;
import exceptions.InvalidColourException;
import exceptions.InvalidConfigException;
import exceptions.InvalidGenomicCoordsException;
import java.io.IOException;
import org.biojava.nbio.core.sequence.transcription.Frame;
import org.junit.Before;
import org.junit.Test;

public class GenomicSequenceTest {

  @Before
  public void initConfig() throws IOException, InvalidConfigException {
    new Config(null);
  }

  @Test
  public void canGetGeneticCodeNames() throws InvalidGenomicCoordsException {
    String dna = "";
    GenomicSequence gs = new GenomicSequence(dna.getBytes(), 1, 10);
    assertEquals(17, gs.geneticCodeNames().size());
    assertEquals("UNIVERSAL", gs.geneticCodeNames().get(1));
  }

  @Test
  public void canTranslateSequenceIupac() throws InvalidColourException, InvalidGenomicCoordsException {
    String dna = "ATGCRGTAG";
    GenomicSequence gs = new GenomicSequence(dna.getBytes(), 4, 12);
    gs.setNoFormat(true);
    gs.setPrintCodon(PrintCodon.ALL);
    gs.setFrames(Frame.getAllFrames());
    assertTrue(gs.getPrintableSequence(dna.length()).contains("\n1M--X--* "));
  }

  @Test
  public void canTranslateSequence() throws InvalidColourException, InvalidGenomicCoordsException {
    String dna = "ATGCTGTAG";
    GenomicSequence gs = new GenomicSequence(dna.getBytes(), 4, 12);
    gs.setNoFormat(true);
    gs.setPrintCodon(PrintCodon.ALL);
    gs.setFrames(Frame.getAllFrames());

    assertTrue(gs.getPrintableSequence(dna.length()).contains("\n1M--L--* "));
    assertTrue(gs.getPrintableSequence(dna.length()).contains("\n2 C  C "));
    assertTrue(gs.getPrintableSequence(dna.length()).startsWith("3  A  V"));

    assertTrue(gs.getPrintableSequence(dna.length()).contains("\n H  Q  L1\n"));
    assertTrue(gs.getPrintableSequence(dna.length()).contains("\n   S  Y 2\n"));
    assertTrue(gs.getPrintableSequence(dna.length()).endsWith("\n  A  T  3\n"));

    gs = new GenomicSequence(dna.getBytes(), 2, 11);
    gs.setNoFormat(true);
    gs.setPrintCodon(PrintCodon.ALL);
    gs.setFrames(Frame.getAllFrames());

    assertTrue(gs.getPrintableSequence(dna.length()).contains("\n1  A  V"));
    assertTrue(gs.getPrintableSequence(dna.length()).contains("\n2M--L--*"));
    assertTrue(gs.getPrintableSequence(dna.length()).startsWith("3 C  C"));

    assertTrue(gs.getPrintableSequence(dna.length()).contains("\n  A  T  1\n"));
    assertTrue(gs.getPrintableSequence(dna.length()).contains("\n H  Q  L2\n"));
    assertTrue(gs.getPrintableSequence(dna.length()).endsWith("\n   S  Y 3\n"));
  }

  @Test
  public void canPrintSequenceSameAsWindowSize()
      throws InvalidColourException, InvalidGenomicCoordsException {
    String dna = "ATGCTGTAGATGCTGTAGATGCTGTAG";
    GenomicSequence gs = new GenomicSequence(dna.getBytes(), 1, dna.length());
    gs.setNoFormat(true);
    gs.setPrintCodon(PrintCodon.ALL);
    gs.setFrames(Frame.ONE);

    // Window size same as DNA length: Print DNA and aminoacids
    String seq = gs.getPrintableSequence(dna.length());
    assertTrue(seq.contains("1M--L--*"));
    assertTrue(seq.contains(dna));

    // Window size larger DNA length: Print DNA and aminoacids
    seq = gs.getPrintableSequence(dna.length() * 2);
    assertTrue(seq.contains("1M--L--*"));
    assertTrue(seq.contains(dna));
  }

  @Test
  public void testAdaptSequenceToWindowSizeForward()
      throws InvalidColourException, InvalidGenomicCoordsException {
    String dna =
        "TAAATG"
            + "n".repeat(300)
            + "atgTAA"
            + "n".repeat(150)
            + "ATGatgATG"
            + "n".repeat(150)
            + "ATGtaaATG"
            + "n".repeat(150)
            + "TAATAATAA"
            + "n".repeat(150)
            + "TAA";
    GenomicSequence gs = new GenomicSequence(dna.getBytes(), 1, dna.length());
    gs.setNoFormat(true);
    gs.setFrames(Frame.ONE);

    String out = gs.getPrintableSequence(dna.length() - 1);
    assertTrue(out.startsWith(" *  M>>>>>>"));
    assertTrue(out.contains(">>>>M>>*"));
    assertTrue(out.contains("*  *  * "));
    assertTrue(out.endsWith(" * \n"));

    out = gs.getPrintableSequence(30);
    assertEquals("$>>>>>>>>*    M>>>>$>>>>*    *\n", out);

    dna = "n" + dna;
    gs = new GenomicSequence(dna.getBytes(), 1, dna.length());
    gs.setNoFormat(true);
    gs.setFrames(Frame.TWO);

    out = gs.getPrintableSequence(dna.length() - 1);
    assertTrue(out.startsWith("  *  M>>>>"));
    assertTrue(out.endsWith(" * \n"));

    dna = "n" + dna;
    gs = new GenomicSequence(dna.getBytes(), 1, dna.length());
    gs.setNoFormat(true);
    gs.setFrames(Frame.THREE);

    out = gs.getPrintableSequence(dna.length() - 1);
    assertTrue(out.startsWith("   *  M>>>"));
    assertTrue(out.endsWith(" * \n"));
  }

  @Test
  public void testAdaptSequenceNotMultipleOfThreeForward()
      throws InvalidColourException, InvalidGenomicCoordsException {
    String dna =
        "TGATTGATCTGCCAAAAGGGGAAGAATGAGTCCAGCTAGAATCCAGGACTAACCAGCGGGTGAGCTTCAAGGAACAAAGGGCTTCCGCTGGGTCAGCCCACGAGAGGGAGCTGCCTGCAGGTACCTGGGAGGGCACAGCCACCGTGTCTGATGCT";
    GenomicSequence gs = new GenomicSequence(dna.getBytes(), 3, 1000);
    gs.setNoFormat(true);
    gs.setFrames(Frame.THREE);
    String out = gs.getPrintableSequence(dna.length() - 1);
    assertTrue(out.startsWith(" *  "));

    dna =
        "nTGATTGATCTGCCAAAAGGGGAAGAATGAGTCCAGCTAGAATCCAGGACTAACCAGCGGGTGAGCTTCAAGGAACAAAGGGCTTCCGCTGGGTCAGCCCACGAGAGGGAGCTGCCTGCAGGTACCTGGGAGGGCACAGCCACCGTGTCTGTTCCT";
    gs = new GenomicSequence(dna.getBytes(), 1, 1000);
    gs.setNoFormat(true);
    gs.setFrames(Frame.TWO);
    out = gs.getPrintableSequence(dna.length() - 1);
    assertTrue(out.startsWith("  *  "));
  }

  @Test
  public void testAdaptSequenceToWindowSizeReverse()
      throws InvalidColourException, InvalidGenomicCoordsException {
    String dna = "TTAnnnTTAnnnCATnnnnnnTTAnnnnnnnnnnnnnnnnnnCATnnnCATnnnCATnn";
    GenomicSequence gs = new GenomicSequence(dna.getBytes(), 1, dna.length());
    gs.setNoFormat(true);
    gs.setFrames(Frame.REVERSED_THREE);
    String out = gs.getPrintableSequence(dna.length() - 1);
    assertEquals(" *     *<<<<<M        *<<<<<<<<<<<<<<<<<<<M<<<<<M<<<<<M   \n", out);
  }

  @Test
  public void initMet() throws InvalidColourException, InvalidGenomicCoordsException {
    /*
    If the transcriptionEngine includes '.initMet(true)' then you tell biojava
    to start with Met if the first codon is a suitable alternative start codon.
    TTG is a possible initiation codon and biojava translates it as Met if
    .initMet(true). However, we set .initMet(false) and we expect Leu instead
    of Met.
     */
    String dna = "TTC CAA TA".replaceAll(" ", "");
    GenomicSequence gs = new GenomicSequence(dna.getBytes(), 1, dna.length());
    gs.setNoFormat(true);
    gs.setFrames(Frame.REVERSED_THREE);

    String out = gs.getPrintableSequence(dna.length());
    assertTrue(out.contains(" E  L "));
  }
}
