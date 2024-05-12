package tracks;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import coloring.Config;
import exceptions.InvalidConfigException;
import exceptions.InvalidGenomicCoordsException;
import exceptions.InvalidRecordException;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.List;
import org.junit.Before;
import org.junit.Test;
import samTextViewer.GenomicCoords;

public class TrackBookmarkTest {

  @Before
  public void prepareConfig() throws IOException, InvalidConfigException {
    new Config(null);
  }

  @Test
  public void canColorByRegex() throws Exception {
    GenomicCoords gc = new GenomicCoords("chr1:1-100", 80, null, null);
    TrackBookmark bm = new TrackBookmark(gc, "book1");
    bm.printToScreen(); // This is to populate the ideograms.

    List<Argument> colorForRegex = new ArrayList<Argument>();
    colorForRegex.add(new Argument(".*", "216", false));
    bm.setColorForRegex(colorForRegex);
    assertTrue(bm.printToScreen().contains("216"));
  }

  @Test
  public void canAddIntervalsToTrackBookmark()
      throws InvalidGenomicCoordsException,
          IOException,
          ClassNotFoundException,
          InvalidRecordException,
          SQLException {

    GenomicCoords gc = new GenomicCoords("chr1:1-100", 80, null, null);
    TrackBookmark bm = new TrackBookmark(gc, "book1");
    bm.setNoFormat(true);

    gc = new GenomicCoords("chr1:400-500", 80, null, null);
    bm.setGc(gc); // Move track to this position w/o adding this pos to bookmarks.

    gc = new GenomicCoords("chr1:200-300", 80, null, null);
    bm.setGc(gc);
    bm.addBookmark(gc, "book2");

    gc = new GenomicCoords("chr1:1-1000", 80, null, null);
    bm.setGc(gc);

    assertTrue(bm.getIntervalFeatureList().size() == 2);

    // Can set name
    assertEquals("book2", bm.getIntervalFeatureList().get(1).getName());
  }

  @Test
  public void canRemoveBookmark()
      throws InvalidGenomicCoordsException,
          IOException,
          ClassNotFoundException,
          InvalidRecordException,
          SQLException {

    GenomicCoords gc = new GenomicCoords("chr1:1-100", 80, null, null);
    TrackBookmark bm = new TrackBookmark(gc, "book1");
    bm.setNoFormat(true);

    gc = new GenomicCoords("chr1:200-300", 80, null, null);
    bm.setGc(gc);
    bm.addBookmark(gc, "book2");

    gc = new GenomicCoords("chr1:1-1000", 80, null, null); // Zoom out to span both bookmarks.
    bm.setGc(gc);
    // Just check we created a correct test object with two intervals
    assertEquals(2, bm.getIntervalFeatureList().size());

    // Partial overlap not touched
    gc = new GenomicCoords("chr1:1-1000", 80, null, null);
    bm.setGc(gc);
    bm.removeBookmark(gc);
    assertEquals(2, bm.getIntervalFeatureList().size());

    // Go to first bookmark and remove it:
    gc = new GenomicCoords("chr1:1-100", 80, null, null);
    bm.setGc(gc);
    bm.removeBookmark(gc);
    bm.setGc(new GenomicCoords("chr1:1-1000", 80, null, null));
    assertEquals(1, bm.getIntervalFeatureList().size());

    // Go to second bookmark and remove it. No bookmarks left:
    gc = new GenomicCoords("chr1:200-300", 80, null, null);
    bm.setGc(gc);
    bm.removeBookmark(gc);
    assertTrue(bm.getIntervalFeatureList().size() == 0);
  }

  @Test
  public void canPrintBookmarksAsList()
      throws InvalidGenomicCoordsException,
          IOException,
          ClassNotFoundException,
          InvalidRecordException,
          SQLException {

    GenomicCoords gc = new GenomicCoords("chr1:1-100", 80, null, null);
    TrackBookmark bm = new TrackBookmark(gc, "book1");
    bm.setNoFormat(true);

    gc = new GenomicCoords("chr1:200-300", 80, null, null);
    bm.setGc(gc);
    bm.addBookmark(gc, "book2");

    gc = new GenomicCoords("chr2:2000-3000", 80, null, null);
    bm.setGc(gc);
    bm.addBookmark(gc, "foo bar");

    assertEquals(3, bm.asList().size());
    assertTrue(bm.asList().get(1).contains("chr1:200-300"));
    assertTrue(bm.asList().get(2).contains("foo bar"));
  }
}
