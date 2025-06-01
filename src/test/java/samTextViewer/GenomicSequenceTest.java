package samTextViewer;

import static org.junit.Assert.assertEquals;

import colouring.Config;
import exceptions.InvalidColourException;
import exceptions.InvalidConfigException;
import exceptions.InvalidGenomicCoordsException;
import java.io.IOException;
import java.util.ArrayList;
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
    GenomicSequence gs = new GenomicSequence(dna.getBytes());
    assertEquals(17, gs.geneticCodeNames().size());
    assertEquals("UNIVERSAL", gs.geneticCodeNames().get(1));
  }

  @Test
  public void canTranslateSequence() throws InvalidColourException, InvalidGenomicCoordsException {
    String dna = "ATGCTGTAG";
    GenomicSequence gs = new GenomicSequence(dna.getBytes());
    // gs.setGeneticCode("UNIVERSAL");
    gs.setNoFormat(false);
    gs.setPrintCodon(PrintCodon.ALL);

    gs.setFrames(Frame.getAllFrames());
    System.out.print(gs.getPrintableSequence() + "---\n");
    gs.setFrames(Frame.getForwardFrames());
    System.out.print(gs.getPrintableSequence() + "---\n");
    gs.setFrames(Frame.getReverseFrames());
    System.out.print(gs.getPrintableSequence() + "---\n");

    gs.setFrames(new ArrayList<>());
    System.out.print(gs.getPrintableSequence() + "---\n");
  }
}
