package tracks;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;

import colouring.Config;
import exceptions.InvalidConfigException;
import exceptions.InvalidGenomicCoordsException;
import exceptions.InvalidRecordException;
import java.io.IOException;
import java.sql.SQLException;
import org.junit.Before;
import org.junit.Test;
import samTextViewer.GenomicCoords;
import utils.CsvFormat;

public class TrackCSVTest {

  @Before
  public void prepareConfig() throws IOException, InvalidConfigException {
    new Config(null);
  }

  @Test
  public void canConstructTrackCSV()
      throws ClassNotFoundException,
      IOException,
      InvalidRecordException,
      InvalidGenomicCoordsException,
      SQLException {

//    GenomicCoords gc = new GenomicCoords("chr1:1-20", 80, null, null);
//    TrackBedgraph tb = new TrackBedgraph("test_data/test.bedGraph", gc);
//     System.err.println(tb.getFeatureList());

    GenomicCoords gc2 = new GenomicCoords("chr1:1-20", 80, null, null);
    CsvFormat csv = new CsvFormat(0, 1, 2, 3, true, 0, '#', '\t');
    TrackCSV tb2 = new TrackCSV("test_data/test.bedGraph", gc2, csv);
    System.err.println(tb2.getFeatureList());
  }

  @Test
  public void canPrintChromosomeNames()
      throws InvalidGenomicCoordsException,
      IOException,
      ClassNotFoundException,
      InvalidRecordException,
      SQLException {

    GenomicCoords gc = new GenomicCoords("chr7:5540000-5570000", 80, null, null);

    TrackBedgraph tb = new TrackBedgraph("test_data/test.bedGraph", gc);
    System.out.println(tb.getTrackFormat());
    assertFalse(tb.getChromosomeNames().isEmpty());
    assertEquals("chr1", tb.getChromosomeNames().get(0));
  }

}
