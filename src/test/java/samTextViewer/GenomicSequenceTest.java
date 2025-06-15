package samTextViewer;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import colouring.Config;
import exceptions.InvalidColourException;
import exceptions.InvalidConfigException;
import exceptions.InvalidGenomicCoordsException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
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

    assertTrue(gs.getPrintableSequence().contains("\n1M  L  * "));
    assertTrue(gs.getPrintableSequence().contains("\n2 C  C "));
    assertTrue(gs.getPrintableSequence().startsWith("3  A  V"));

    assertTrue(gs.getPrintableSequence().contains("\n H  Q  L1\n"));
    assertTrue(gs.getPrintableSequence().contains("\n   S  Y 2\n"));
    assertTrue(gs.getPrintableSequence().endsWith("\n  A  T  3\n"));

    gs = new GenomicSequence(dna.getBytes(), 2, 11);
    gs.setNoFormat(true);
    gs.setPrintCodon(PrintCodon.ALL);
    gs.setFrames(Frame.getAllFrames());

    assertTrue(gs.getPrintableSequence().contains("\n1  A  V"));
    assertTrue(gs.getPrintableSequence().contains("\n2M  L  *"));
    assertTrue(gs.getPrintableSequence().startsWith("3 C  C"));

    assertTrue(gs.getPrintableSequence().contains("\n  A  T  1\n"));
    assertTrue(gs.getPrintableSequence().contains("\n H  Q  L2\n"));
    assertTrue(gs.getPrintableSequence().endsWith("\n   S  Y 3\n"));

    System.err.println(gs.getPrintableSequence());
  }
}