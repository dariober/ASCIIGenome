package tracks;

import static org.junit.Assert.*;

import java.io.IOException;
import java.sql.SQLException;

import org.junit.Test;

import exceptions.InvalidGenomicCoordsException;
import exceptions.InvalidRecordException;
import samTextViewer.GenomicCoords;

public class TrackSeqRegexTest {

	@Test
	public void canInitializeTrack() throws InvalidGenomicCoordsException, IOException, ClassNotFoundException, InvalidRecordException, SQLException {
		
		GenomicCoords gc= new GenomicCoords("chr7:8000000-8001000", null, "test_data/chr7.fa");
		
		TrackSeqRegex trackSeqRegex= new TrackSeqRegex(gc); 
		trackSeqRegex.setSeqRegex("(?i)CC");
		trackSeqRegex.setNoFormat(true);
		assertTrue(trackSeqRegex.getIntervalFeatureList().size() > 30);
		assertTrue(trackSeqRegex.getIntervalFeatureList().size() < 300);
		
		gc= new GenomicCoords("chr7:8000000-8000050", null, "test_data/chr7.fa");
		trackSeqRegex.setGc(gc);
		assertTrue(trackSeqRegex.getIntervalFeatureList().size() < 30);
		assertTrue(trackSeqRegex.getIntervalFeatureList().size() > 5);
		
	}

	@Test
	public void methodsWork() throws InvalidGenomicCoordsException, IOException, ClassNotFoundException, InvalidRecordException, SQLException{
		
		GenomicCoords gc= new GenomicCoords("chr7:8000000-8001000", null, "test_data/chr7.fa");
		
		TrackSeqRegex trackSeqRegex= new TrackSeqRegex(gc);
		trackSeqRegex.setSeqRegex("(?i)CC");		
		trackSeqRegex.setNoFormat(true);
		trackSeqRegex.setGap(0); // Set gap mode
		
		// Change track height
		trackSeqRegex.setyMaxLines(3);
		assertTrue(trackSeqRegex.printToScreen().split("\n").length == 3);
		
		trackSeqRegex.setyMaxLines(10);
		assertTrue(trackSeqRegex.printToScreen().split("\n").length > 3);
		
		// Test title
		assertTrue(!trackSeqRegex.getTitle().isEmpty());
		
		// FIXME: Filters
//		trackSeqRegex.setHideRegex("chr7");
//		trackSeqRegex.setPrintMode(PrintRawLine.CLIP);
//		System.out.println(trackSeqRegex.printFeatures(100));
		
	}
	
	@Test
	public void canReturnTrackFromMatchingRegex() throws InvalidGenomicCoordsException, IOException, ClassNotFoundException, InvalidRecordException, SQLException{

		GenomicCoords gc= new GenomicCoords("seq:1-100", null, "test_data/seq_cg.fa");
		TrackSeqRegex trackSeqRegex= new TrackSeqRegex(gc);
		trackSeqRegex.setNoFormat(true);
		trackSeqRegex.setSeqRegex("(?i)atc");
		
		assertTrue(trackSeqRegex.printToScreen().startsWith(">>>"));
		
		
		// Match not found
		trackSeqRegex.setSeqRegex("FOOBAR");
		assertTrue(trackSeqRegex.printToScreen().isEmpty());
		
		// Regex not given
		trackSeqRegex.setSeqRegex("");
		assertTrue(trackSeqRegex.printToScreen().isEmpty());
		
		// Regex unset at init
		trackSeqRegex= new TrackSeqRegex(gc);
		trackSeqRegex.setNoFormat(true);
		
		// Palindromic
		gc= new GenomicCoords("seq:1-100", null, "test_data/seq_cg.fa");
		trackSeqRegex.setNoFormat(true);
		trackSeqRegex.setSeqRegex("cG"); // Case insensitive
		assertTrue(trackSeqRegex.printToScreen().trim().startsWith("|"));
		
		// Case sensitive
		trackSeqRegex.setCaseSensitive(true);
		trackSeqRegex.setSeqRegex("cG");
		assertTrue(trackSeqRegex.getIntervalFeatureList().size() == 0);
	}
	
	@Test 
	public void handlingMissingFastaFile() throws ClassNotFoundException, IOException, InvalidGenomicCoordsException, InvalidRecordException, SQLException{
		// Ref Sequence not given 
		boolean passed= false;
		try{
			new TrackSeqRegex(new GenomicCoords("seq:1-100", null, null));
		} catch (NullPointerException e) {
			passed= true;
		}
		assertTrue(passed);
	}

	
}