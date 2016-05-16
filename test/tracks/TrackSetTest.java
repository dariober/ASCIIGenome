package tracks;

import static org.junit.Assert.*;

import java.io.IOException;
import org.junit.Test;

import exceptions.InvalidCommandLineException;
import exceptions.InvalidGenomicCoordsException;
import samTextViewer.GenomicCoords;

public class TrackSetTest {

	@Test
	public void canSetVisibilityForTrackIntervalFeature() throws InvalidCommandLineException, IOException, InvalidGenomicCoordsException{
				
		TrackSet ts= new TrackSet();
		GenomicCoords gc= new GenomicCoords("chr1", 1, 100, null, 100, null);
		Track t1= new TrackIntervalFeature("test_data/hg19_genes_head.gtf", gc); t1.setFileTag("#10"); ts.addOrReplace(t1);
		Track t2= new TrackIntervalFeature("test_data/hg19_genes_head.gtf.gz", gc); t2.setFileTag("#11"); ts.addOrReplace(t2);
		Track t3= new TrackIntervalFeature("test_data/refSeq.bed", gc); t3.setFileTag("#30"); ts.addOrReplace(t3);
		
		String cmdInput= "visible exon intron .*#1.*"; // Set for #1...
		ts.setVisibilityForTrackIntervalFeature(cmdInput);

		assertEquals("exon", ts.getTrackSet().get("#10").getShowRegex());
		assertEquals("intron", ts.getTrackSet().get("#11").getHideRegex());

		assertEquals(".*", ts.getTrackSet().get("#30").getShowRegex()); // As default
		assertEquals("", ts.getTrackSet().get("#30").getHideRegex());

		// assertEquals(0, ts.getTrackSet().get("#20").getYmin(), 0.001);
		// assertEquals(10, ts.getTrackSet().get("#1").getYmax(), 0.001);
	}
		
	@Test
	public void canSetTrackHeight() throws InvalidCommandLineException, IOException, InvalidGenomicCoordsException{
				
		TrackSet ts= new TrackSet();
		Track t1= new Track(); t1.setFilename("foo.gz"); t1.setFileTag("#1"); ts.addOrReplace(t1);
		Track t2= new Track(); t2.setFilename("foo.txt"); t2.setFileTag("#20"); ts.addOrReplace(t2);
		Track t3= new Track(); t3.setFilename("bla.gz"); t3.setFileTag("#3"); ts.addOrReplace(t3);

		String cmdInput= "trackHeight 2 #\\d";
		ts.setTrackHeightForRegex(cmdInput);
				
		assertEquals(2, ts.getTrackSet().get("#1").getyMaxLines());
		assertEquals(10, ts.getTrackSet().get("#20").getyMaxLines()); 
		assertEquals(2, ts.getTrackSet().get("#1").getyMaxLines());

	}
	
	@Test
	public void canSetYlimits() throws InvalidCommandLineException{
				
		String cmdInput= "ylim 0 10 #\\d+";
		// List<Track> tracks= new ArrayList<Track>();
		
		TrackSet ts= new TrackSet();
		
		Track t1= new Track(); t1.setFilename("foo.gz"); t1.setFileTag("#1"); ts.addOrReplace(t1);
		Track t2= new Track(); t2.setFilename("foo.txt"); t2.setFileTag("#20"); ts.addOrReplace(t2);
		Track t3= new Track(); t3.setFilename("bla.gz"); t3.setFileTag("#3"); ts.addOrReplace(t3);
		
		ts.setTrackYlimitsForRegex(cmdInput);
				
		assertEquals(0, ts.getTrackSet().get("#1").getYLimitMin(), 0.001);
		assertEquals(0, ts.getTrackSet().get("#20").getYLimitMin(), 0.001);
		assertEquals(10, ts.getTrackSet().get("#1").getYLimitMax(), 0.001);
	}
	
}
