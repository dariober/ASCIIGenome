package tracks;

import static org.junit.Assert.assertTrue;

import coloring.Config;
import exceptions.InvalidConfigException;
import exceptions.InvalidGenomicCoordsException;
import exceptions.InvalidRecordException;
import java.io.IOException;
import java.sql.SQLException;
import org.junit.Before;
import org.junit.Test;
import samTextViewer.GenomicCoords;

public class TrackWigglesTest {

  @Before
  public void prepareConfig() throws IOException, InvalidConfigException {
    new Config(null);
  }

  @Test
  public void canCloseReaders()
      throws ClassNotFoundException,
          IOException,
          InvalidRecordException,
          InvalidGenomicCoordsException,
          SQLException {
    GenomicCoords gc = new GenomicCoords("chr7:5540000-5570000", 80, null, null);
    TrackWiggles tw = new TrackWiggles("test_data/hg18_var_sample.wig.v2.1.30.tdf", gc);
    tw.close();
  }

  @Test
  public void canPrintChromosomeNames()
      throws InvalidGenomicCoordsException,
          IOException,
          ClassNotFoundException,
          InvalidRecordException,
          SQLException {

    GenomicCoords gc = new GenomicCoords("chr7:5540000-5570000", 80, null, null);
    TrackWiggles tw = new TrackWiggles("test_data/hg18_var_sample.wig.v2.1.30.tdf", gc);
    assertTrue(tw.getChromosomeNames().size() > 10);
  }

  @Test
  /** Snippet to extract totalCount from TDF, useful for normalizing signal. */
  public void canNomrmalizeTDFtoRPM()
      throws InvalidGenomicCoordsException,
          IOException,
          InvalidRecordException,
          ClassNotFoundException,
          SQLException {
    GenomicCoords gc = new GenomicCoords("chr7:5540000-5570000", 80, null, null);
    TrackWiggles tw = new TrackWiggles("test_data/ear045.oxBS.actb.tdf", gc);
    Float raw = tw.getScreenScores().get(0);
    tw.setRpm(true);
    Float rpm = tw.getScreenScores().get(0);
    assertTrue(rpm > raw);
  }
}
