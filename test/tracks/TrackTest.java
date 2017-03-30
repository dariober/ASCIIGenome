package tracks;

import static org.junit.Assert.assertTrue;

import java.io.IOException;
import java.sql.SQLException;

import org.junit.Test;

import coloring.Config;
import exceptions.InvalidColourException;
import exceptions.InvalidCommandLineException;
import exceptions.InvalidConfigException;
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
	
	@Test
	public void canCutLinesWhenPrinting() throws InvalidGenomicCoordsException, IOException, ClassNotFoundException, InvalidRecordException, SQLException, InvalidColourException, InvalidConfigException, InvalidCommandLineException{

		new Config(null);
		
		GenomicCoords gc= new GenomicCoords("chr1:1-100000", null, null);
		TrackIntervalFeature tif= new TrackIntervalFeature("test_data/hg19_genes_head.gtf", gc);
		tif.setNoFormat(true);
		tif.setPrintMode(PrintRawLine.FULL);
		tif.setCutScriptForPrinting("-d \\t|; -f 3-5,,11-9"); // Note reversing 11,10,9
		String printed= tif.printLines();
		assertTrue(printed.length() > 200);
		assertTrue(printed.startsWith("exon"));
		assertTrue(printed.trim().endsWith("\"OR4F5\""));
		
		// Columns can be duplicated!
		tif.setCutScriptForPrinting("-d \\t|; -f 1,1,2,2");
		printed= tif.printLines();
		assertTrue(printed.startsWith("chr1 chr1"));
		System.err.println(tif.printLines());
	}

}
