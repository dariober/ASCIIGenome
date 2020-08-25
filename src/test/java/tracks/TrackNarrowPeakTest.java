package tracks;

import static org.junit.Assert.*;

import java.io.IOException;
import java.sql.SQLException;

import org.junit.Before;
import org.junit.Test;

import coloring.Config;
import exceptions.InvalidColourException;
import exceptions.InvalidConfigException;
import exceptions.InvalidGenomicCoordsException;
import exceptions.InvalidRecordException;
import samTextViewer.GenomicCoords;

public class TrackNarrowPeakTest {

	@Before
	public void prepareConfig() throws IOException, InvalidConfigException{
		new Config(null);
	}

	@Test
	public void testCanInitNarrorPeakTrack() throws InvalidGenomicCoordsException, IOException, ClassNotFoundException, InvalidRecordException, SQLException, InvalidColourException {

		GenomicCoords gc= new GenomicCoords("chr1:4389687-26338119", 80, null, null);
		Track np= new TrackNarrowPeak("test_data/wgEncodeCaltechTfbsC2c12CebpbFCntrl50bE2p60hPcr1xPkRep2.narrowPeak.gz", gc);
		np.setNoFormat(true);
		np.setYLimitMin(0);
		np.setYLimitMin(10);
		np.setyMaxLines(5);
		String x = np.printToScreen();
		System.err.println(x);
		assertTrue(x.length() > 10);
	}

	// Test colour according to fold_change value [col #7]
	
	// Test set colour using featureColor
}
