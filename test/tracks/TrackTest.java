package tracks;

import static org.junit.Assert.*;

import java.io.IOException;
import java.sql.SQLException;

import org.junit.Test;

import exceptions.InvalidGenomicCoordsException;
import exceptions.InvalidRecordException;
import samTextViewer.GenomicCoords;

public class TrackTest {

	@Test
	public void canExportSettings() throws InvalidGenomicCoordsException, IOException, ClassNotFoundException, InvalidRecordException, SQLException {
		
		String bgzFn= "test_data/refSeq.hg19.short.sort.bed.gz"; // "test_data/refSeq.hg19.short.sort.bed.gz";
		GenomicCoords gc= new GenomicCoords("chr1:16000000-20000000", null, null);
		TrackIntervalFeature tif= new TrackIntervalFeature(bgzFn, gc);
		tif.setTrackTag("name.bed#12");
		System.out.println(tif.settingsToString());
		
	}

}
