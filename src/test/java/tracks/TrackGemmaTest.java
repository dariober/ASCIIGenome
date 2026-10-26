package tracks;

import colouring.Config;
import exceptions.*;
import org.junit.Before;
import org.junit.Test;
import samTextViewer.GenomicCoords;

import java.io.IOException;
import java.sql.SQLException;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

public class TrackGemmaTest {
  @Before
  public void prepareConfig() throws IOException, InvalidConfigException {
    new Config(null);
  }

  // Format of Gemma prefix.assoc.txt (from manual):
  //  The 11 columns are: chromosome numbers, snp ids, base pair positions on the chromosome,
  //  number of missing individuals for a given snp, number of non-missing individuals for a given snp,
  //  minor allele, major allele, allele frequency, beta estimates (additive eﬀect), standard errors for beta,
  //  and p values from the Wald test.
  @Test
  public void canReadAndPlotProfile()
          throws InvalidGenomicCoordsException,
          IOException,
          ClassNotFoundException,
          InvalidRecordException,
          SQLException,
          InvalidColourException, InvalidCommandLineException {

    GenomicCoords gc = new GenomicCoords("4:36494-1678835", 80, null, null);
    TrackGemma tg = new TrackGemma("test_data/gemma_ebi.assoc.txt.gz", gc);
    tg.setNoFormat(true);
    assertTrue(tg.concatTitleAndTrack().contains("range[1.03 9.08]"));
    assertTrue(tg.concatTitleAndTrack().contains(":::::::::::::::::"));

    tg.setPrintMode(PrintRawLine.FULL);
    tg.setPrintRawLineCount(10);
    assertTrue(tg.printLines().contains("4:36494:G:A"));

    tg.setScoreColIdx(8); // Beta column
    tg.setDataTransformation(DataTransformation.IDENTITY);
    tg.setDataAggregationMethod(DataAggregationMethod.ABS_MAX);
    assertTrue(tg.concatTitleAndTrack().contains("range[-0.197 0.427]"));

    tg.getGc().setChrom("1");
    tg.getGc().setFrom(1005721);
    tg.getGc().setTo(1005772);
    tg.getGc().setTerminalWidth(80);
    tg.setDataTransformation(DataTransformation.IDENTITY);
    tg.update();
    assertTrue(tg.concatTitleAndTrack().contains("range[-0.00723 0.00244]"));
    assertTrue(tg.concatTitleAndTrack().contains(",                                                  :"));
  }
}
