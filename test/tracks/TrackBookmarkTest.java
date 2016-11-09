package tracks;

import static org.junit.Assert.*;

import java.io.IOException;
import java.sql.SQLException;

import org.junit.Test;

import exceptions.InvalidGenomicCoordsException;
import exceptions.InvalidRecordException;
import samTextViewer.GenomicCoords;

public class TrackBookmarkTest {

	@Test
	public void canAddIntervalsToTrackBookmark() throws InvalidGenomicCoordsException, IOException, ClassNotFoundException, InvalidRecordException, SQLException {
		
		GenomicCoords gc= new GenomicCoords("chr1:1-100", null, null);
		TrackBookmark bm= new TrackBookmark(gc, "book1");
		bm.setNoFormat(true);
		
		gc= new GenomicCoords("chr1:400-500", null, null);
		bm.setGc(gc); // Move track to this position w/o adding this pos to bookmarks.

		gc= new GenomicCoords("chr1:200-300", null, null);
		bm.setGc(gc);
		bm.addBookmark("book2");
		
		gc= new GenomicCoords("chr1:1-1000", null, null);
		bm.setGc(gc);
		
		assertTrue(bm.getIntervalFeatureList().size() == 2);
		
		// Can set name
		assertEquals("book2", bm.getIntervalFeatureList().get(1).getName());
		System.out.println(bm.printToScreen());
	}
	
	@Test
	public void canRemoveBookmark() throws InvalidGenomicCoordsException, IOException, ClassNotFoundException, InvalidRecordException, SQLException{

		GenomicCoords gc= new GenomicCoords("chr1:1-100", null, null);
		TrackBookmark bm= new TrackBookmark(gc, "book1");
		bm.setNoFormat(true);
		
		gc= new GenomicCoords("chr1:200-300", null, null);
		bm.setGc(gc);
		bm.addBookmark("book2");
		
		gc= new GenomicCoords("chr1:1-1000", null, null); // Zoom out to span both bookmarks.  
		bm.setGc(gc);
		// Just check we created a correct test object with two intervals
		assertEquals(2, bm.getIntervalFeatureList().size());

		// Partial overlap not touched 
		bm.setGc(new GenomicCoords("chr1:1-1000", null, null));
		bm.removeBookmark();
		assertEquals(2, bm.getIntervalFeatureList().size());

		// Go to first bookmark and remove it:
		bm.setGc(new GenomicCoords("chr1:1-100", null, null));
		bm.removeBookmark();
		bm.setGc(new GenomicCoords("chr1:1-1000", null, null));
		assertEquals(1, bm.getIntervalFeatureList().size());

		// Go to second bookmark and remove it. No bookmarks left:
		bm.setGc(new GenomicCoords("chr1:200-300", null, null));
		bm.removeBookmark();
		assertTrue(bm.getIntervalFeatureList().size() == 0);
	}

	@Test
	public void canExportSetting() throws InvalidGenomicCoordsException, IOException, ClassNotFoundException, InvalidRecordException, SQLException{

		GenomicCoords gc= new GenomicCoords("chr1:1-100", null, null);
		TrackBookmark bm= new TrackBookmark(gc, "book1");
		bm.setNoFormat(true);
		
		gc= new GenomicCoords("chr1:200-300", null, null);
		bm.setGc(gc);
		bm.addBookmark("book2");

		gc= new GenomicCoords("chr2:2000-3000", null, null);
		bm.setGc(gc);
		bm.addBookmark(".");
		
		assertTrue(bm.settingsToString().contains("goto chr1:1-100"));
		assertTrue(bm.settingsToString().contains("goto chr2:2000-3000"));
		
		// System.out.println(bm.settingsToString());
	}
	
	@Test
	public void canPrintBookmarksAsList() throws InvalidGenomicCoordsException, IOException, ClassNotFoundException, InvalidRecordException, SQLException{

		GenomicCoords gc= new GenomicCoords("chr1:1-100", null, null);
		TrackBookmark bm= new TrackBookmark(gc, "book1");
		bm.setNoFormat(true);
		
		gc= new GenomicCoords("chr1:200-300", null, null);
		bm.setGc(gc);
		bm.addBookmark("book2");

		gc= new GenomicCoords("chr2:2000-3000", null, null);
		bm.setGc(gc);
		bm.addBookmark("foo bar");
		
		assertEquals(3, bm.asList().size());
		assertTrue(bm.asList().get(1).contains("chr1:200-300"));
		assertTrue(bm.asList().get(2).contains("foo bar"));
	}
	
}
