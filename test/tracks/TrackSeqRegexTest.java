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
		
		GenomicCoords gc= new GenomicCoords("chr7:8000000-8001000", null, 80, "test_data/chr7.fa");
		
		TrackSeqRegex trackSeqRegex= new TrackSeqRegex(gc); 
		trackSeqRegex.setSeqRegex("(?i)CC");
		trackSeqRegex.update();
		trackSeqRegex.setNoFormat(true);
		assertTrue(trackSeqRegex.getIntervalFeatureList().size() > 30);
		assertTrue(trackSeqRegex.getIntervalFeatureList().size() < 300);
		
		gc= new GenomicCoords("chr7:8000000-8000050", null, 80, "test_data/chr7.fa");
		trackSeqRegex.setGc(gc);
		trackSeqRegex.update();
		assertTrue(trackSeqRegex.getIntervalFeatureList().size() < 30);
		assertTrue(trackSeqRegex.getIntervalFeatureList().size() > 5);
		
	}

	@Test
	public void methodsWork() throws InvalidGenomicCoordsException, IOException, ClassNotFoundException, InvalidRecordException, SQLException{
		
		GenomicCoords gc= new GenomicCoords("chr7:8000000-8001000", null, 80, "test_data/chr7.fa");
		
		TrackSeqRegex trackSeqRegex= new TrackSeqRegex(gc);
		trackSeqRegex.setSeqRegex("(?i)CC");		
		trackSeqRegex.update();
		trackSeqRegex.setNoFormat(true);
		trackSeqRegex.setGap(0); // Set gap mode
		
		// Change track height
		trackSeqRegex.setyMaxLines(3);
		assertTrue(trackSeqRegex.printToScreen().split("\n").length == 3);
		
		trackSeqRegex.setyMaxLines(10);
		assertTrue(trackSeqRegex.printToScreen().split("\n").length > 3);
		
		// Test title
		assertTrue(!trackSeqRegex.getTitle().isEmpty());
	}
	
	@Test
	public void canReturnTrackFromMatchingRegex() throws InvalidGenomicCoordsException, IOException, ClassNotFoundException, InvalidRecordException, SQLException{

		GenomicCoords gc= new GenomicCoords("seq", 1, 100, null, 100, "test_data/seq_cg.fa");
		TrackSeqRegex trackSeqRegex= new TrackSeqRegex(gc);
		trackSeqRegex.setSeqRegex("(?i)atc");
		trackSeqRegex.setNoFormat(true);
		trackSeqRegex.update();
		
		assertTrue(trackSeqRegex.printToScreen().startsWith(">>>"));
		
		
		// Match not found
		trackSeqRegex.setSeqRegex("FOOBAR");
		trackSeqRegex.update();
		assertTrue(trackSeqRegex.printToScreen().isEmpty());
		
		// Regex not given
		trackSeqRegex.setSeqRegex("");
		trackSeqRegex.update();
		assertTrue(trackSeqRegex.printToScreen().isEmpty());
		
		// Regex unset at init
		trackSeqRegex= new TrackSeqRegex(gc);
		trackSeqRegex.update();
		trackSeqRegex.setNoFormat(true);
		
		// Palindromic
		gc= new GenomicCoords("seq", 1, 100, null, 100, "test_data/seq_cg.fa");
		trackSeqRegex.setSeqRegex("CG");
		trackSeqRegex.setNoFormat(true);
		trackSeqRegex.update();
		assertTrue(trackSeqRegex.printToScreen().trim().startsWith("|"));
//		System.out.println(trackSeqRegex.printToScreen());
//		System.out.println(gc.printableRefSeq(true));
	}
	
	@Test 
	public void handlingMissingFastaFile() throws ClassNotFoundException, IOException, InvalidGenomicCoordsException, InvalidRecordException, SQLException{
		// Ref Sequence not given 
		boolean passed= false;
		try{
			new TrackSeqRegex(new GenomicCoords("seq", 1, 100, null, 100, null));
		} catch (NullPointerException e) {
			passed= true;
		}
		assertTrue(passed);
	}

	
}
