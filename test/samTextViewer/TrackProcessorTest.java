package samTextViewer;

import static org.junit.Assert.assertTrue;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.List;

import org.junit.Test;

import exceptions.BamIndexNotFoundException;
import exceptions.InvalidGenomicCoordsException;
import exceptions.InvalidRecordException;
import tracks.TrackSet;

public class TrackProcessorTest {
		
	@Test
	public void canSaveTrackSettings() throws InvalidGenomicCoordsException, IOException, ClassNotFoundException, BamIndexNotFoundException, InvalidRecordException, SQLException{
		
		GenomicCoords gc= new GenomicCoords("chr7:1-100", null, null);
		List<String> genome= new ArrayList<String>();
		genome.add("test_data/ear045.oxBS.actb.bam");
		gc.setGenome(genome, true);

		GenomicCoordsHistory gch= new GenomicCoordsHistory();
		gch.add(gc);
		
		TrackSet trackSet= new TrackSet();
		trackSet.addTrackFromSource("test_data/hg19_genes.gtf.gz", gc, null);
		trackSet.addTrackFromSource("test_data/ear045.oxBS.actb.bam", gc, null);
		trackSet.addTrackFromSource("test_data/ear045.oxBS.actb.tdf", gc, null);
		
		TrackProcessor tp= new TrackProcessor(trackSet, gch); 
		
		File x= new File("deleteme.txt");
		x.deleteOnExit();
		tp.exportTrackSetSettings(x.getAbsolutePath());
		assertTrue(x.length() > 200); // File not empty.
		
	}

}
