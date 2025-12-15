package samTextViewer;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
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
  public void canTranslateSequence() throws InvalidColourException, InvalidGenomicCoordsException {
    String dna = "ATGCTGTAG";
    GenomicSequence gs = new GenomicSequence(dna.getBytes(), 4, 12);
    gs.setNoFormat(true);
    gs.setPrintCodon(PrintCodon.ALL);
    gs.setFrames(Frame.getAllFrames());

    assertTrue(gs.getPrintableSequence(dna.length()).contains("\n1M  L  * "));
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
    assertTrue(gs.getPrintableSequence(dna.length()).contains("\n2M  L  *"));
    assertTrue(gs.getPrintableSequence(dna.length()).startsWith("3 C  C"));

    assertTrue(gs.getPrintableSequence(dna.length()).contains("\n  A  T  1\n"));
    assertTrue(gs.getPrintableSequence(dna.length()).contains("\n H  Q  L2\n"));
    assertTrue(gs.getPrintableSequence(dna.length()).endsWith("\n   S  Y 3\n"));
  }

  @Test
  public void canAdaptSequenceToWindowSize()
      throws InvalidColourException, InvalidGenomicCoordsException {
    String dna = "ATGCTGTAGATGCTGTAGATGCTGTAG";
    GenomicSequence gs = new GenomicSequence(dna.getBytes(), 1, dna.length());
    gs.setNoFormat(true);
    gs.setPrintCodon(PrintCodon.ALL);
    gs.setFrames(Frame.ONE);

    // Window size same as DNA length: Print DNA and aminoacids
    String seq = gs.getPrintableSequence(dna.length());
    assertTrue(seq.contains("1M  L  *"));
    assertTrue(seq.contains(dna));

    // Window size larger DNA length: Print DNA and aminoacids
    seq = gs.getPrintableSequence(dna.length() * 2);
    assertTrue(seq.contains("1M  L  *"));
    assertTrue(seq.contains(dna));

    // Window size smaller DNA length (1 char on screen spans more than 1 nt)
    seq = gs.getPrintableSequence(10);
    assertFalse(seq.contains("ATG"));
    assertFalse(seq.contains("M"));
  }
}
