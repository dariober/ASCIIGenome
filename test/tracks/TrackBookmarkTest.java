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
		
		gc= new GenomicCoords("chr1:200-300", null, null);
		bm.setGc(gc);
		bm.add("book2");
		
		gc= new GenomicCoords("chr1:400-500", null, null);
		bm.setGc(gc); // Set it but do not add to bookamrks
		
		gc= new GenomicCoords("chr1:1-1000", null, null);
		bm.setGc(gc);
		
		assertTrue(bm.getIntervalFeatureList().size() == 2);
		
		// Can set name
		assertEquals("book2", bm.getIntervalFeatureList().get(1).getName());
		System.out.println(bm.printToScreen());
	}

}
