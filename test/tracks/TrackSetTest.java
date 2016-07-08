package tracks;

import static org.junit.Assert.*;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.junit.Test;

import exceptions.InvalidCommandLineException;
import exceptions.InvalidGenomicCoordsException;
import samTextViewer.GenomicCoords;

public class TrackSetTest {

	@Test
	public void canReorderTracks() throws InvalidGenomicCoordsException, IOException{
		TrackSet ts= new TrackSet();
		GenomicCoords gc= new GenomicCoords("chr1", 1, 100, null, 100, null);
	
		Track t1= new TrackIntervalFeature("test_data/refSeq.bed", gc); t1.setFileTag("#1"); ts.getTrackSet().put(t1.getFileTag(), t1);
		Track t2= new TrackIntervalFeature("test_data/refSeq.bed", gc); t2.setFileTag("#2"); ts.getTrackSet().put(t2.getFileTag(), t2);
		Track t3= new TrackIntervalFeature("test_data/refSeq.bed", gc); t3.setFileTag("#3"); ts.getTrackSet().put(t3.getFileTag(), t3);
				
		List<String> newOrder= new ArrayList<String>();
		newOrder.add("#1");
		newOrder.add("#2");
		newOrder.add("#3");
		ts.orderTracks(newOrder);
		assertEquals(newOrder, new ArrayList<String>(ts.getTrackSet().keySet()));

		// Handle missing tracks in new order
		newOrder= new ArrayList<String>();
		newOrder.add("#3");
		newOrder.add("#1");
		ts.orderTracks(newOrder);
		assertEquals(3, ts.getTrackSet().keySet().size());
		
		// Handle non existing tracks
		newOrder= new ArrayList<String>();
		ts.orderTracks(newOrder);
		assertEquals(3, ts.getTrackSet().keySet().size());

		newOrder= new ArrayList<String>();
		newOrder.add("1");
		newOrder.add("1");
		newOrder.add("2");
		newOrder.add("3");
		newOrder.add("3");
		newOrder.add("foo!");
		ts.orderTracks(newOrder);
		assertEquals(3, ts.getTrackSet().keySet().size());

		// Partial matches
		ts= new TrackSet();
		t1= new TrackIntervalFeature("test_data/refSeq.bed", gc); t1.setFileTag("#1"); ts.getTrackSet().put(t1.getFileTag(), t1);
		t2= new TrackIntervalFeature("test_data/refSeq.bed", gc); t2.setFileTag("#2"); ts.getTrackSet().put(t2.getFileTag(), t2);
		t3= new TrackIntervalFeature("test_data/refSeq.bed", gc); t3.setFileTag("#3"); ts.getTrackSet().put(t3.getFileTag(), t3);
				
		newOrder= new ArrayList<String>();
		newOrder.add("2");
		newOrder.add("3");
		newOrder.add("1");
		ts.orderTracks(newOrder);
		assertEquals("#2", new ArrayList<String>(ts.getTrackSet().keySet()).get(0));
		assertEquals("#3", new ArrayList<String>(ts.getTrackSet().keySet()).get(1));
	}
	
	@Test
	public void canSetVisibilityForTrackIntervalFeature() throws InvalidCommandLineException, IOException, InvalidGenomicCoordsException{
				
		TrackSet ts= new TrackSet();
		GenomicCoords gc= new GenomicCoords("chr1", 1, 100, null, 100, null);
		Track t1= new TrackIntervalFeature("test_data/hg19_genes_head.gtf", gc); t1.setFileTag("#10"); ts.getTrackSet().put(t1.getFileTag(), t1);
		Track t2= new TrackIntervalFeature("test_data/hg19_genes_head.gtf.gz", gc); t2.setFileTag("#11"); ts.getTrackSet().put(t2.getFileTag(), t2);
		Track t3= new TrackIntervalFeature("test_data/refSeq.bed", gc); t3.setFileTag("#30"); ts.getTrackSet().put(t3.getFileTag(), t3);
		
		String cmdInput= "visible exon intron .*#1.*"; // Set for #1...
		ts.setVisibilityForTrackIntervalFeature(cmdInput);

		assertEquals("exon", ts.getTrackSet().get("#10").getShowRegex());
		assertEquals("intron", ts.getTrackSet().get("#11").getHideRegex());

		assertEquals(".*", ts.getTrackSet().get("#30").getShowRegex()); // As default
		assertEquals("", ts.getTrackSet().get("#30").getHideRegex());

		// assertEquals(0, ts.getTrackSet().get("#20").getYmin(), 0.001);
		// assertEquals(10, ts.getTrackSet().get("#1").getYmax(), 0.001);
	}

	/*
	@Test
	public void canSetPrintFeaturesForTrackIntervalFeature() throws InvalidCommandLineException, IOException, InvalidGenomicCoordsException{
				
		TrackSet ts= new TrackSet();
		GenomicCoords gc= new GenomicCoords("chr1", 1, 20000, null, 100, null);
		Track t1= new TrackIntervalFeature("test_data/hg19_genes_head.gtf", gc); t1.setFileTag("#10"); ts.getTrackSet().put(t1.getFileTag(), t1);
		Track t2= new TrackIntervalFeature("test_data/hg19_genes_head.gtf.gz", gc); t2.setFileTag("#11"); ts.getTrackSet().put(t2.getFileTag(), t2);
		Track t3= new TrackIntervalFeature("test_data/refSeq.bed", gc); t3.setFileTag("#30"); ts.getTrackSet().put(t3.getFileTag(), t3);
		
		String cmdInput= "visible exon intron .*#1.*"; // Set for #1...
		ts.setPrintFeatureForTrackIntervalFeature(cmdInput, 100);

		assertTrue(ts.getTrackSet().get("#10").printFeatures(100).startsWith("chr1"));
	}
	*/

	@Test
	public void canSetBSMode() throws InvalidCommandLineException, IOException, InvalidGenomicCoordsException{
				
		TrackSet ts= new TrackSet();
		Track t1= new Track(); t1.setFilename("foo.gz"); t1.setFileTag("#1"); ts.getTrackSet().put(t1.getFileTag(), t1);
		Track t2= new Track(); t2.setFilename("foo.txt"); t2.setFileTag("#20"); ts.getTrackSet().put(t2.getFileTag(), t2);
		Track t3= new Track(); t3.setFilename("bla.gz"); t3.setFileTag("#3"); ts.getTrackSet().put(t3.getFileTag(), t3);

		String cmdInput= "bsMode #\\d$";
		ts.setBisulfiteModeForRegex(cmdInput);
		assertTrue(ts.getTrackSet().get("#1").isBisulf());
		assertTrue(! ts.getTrackSet().get("#20").isBisulf());
		
		ts.setBisulfiteModeForRegex(cmdInput); // Was set true, now becomes false
		assertTrue(! ts.getTrackSet().get("#1").isBisulf());
	}

	
	@Test
	public void canSetTrackColour() throws InvalidCommandLineException, IOException, InvalidGenomicCoordsException{
				
		TrackSet ts= new TrackSet();
		Track t1= new Track(); t1.setFilename("foo.gz"); t1.setFileTag("#1"); ts.getTrackSet().put(t1.getFileTag(), t1);
		Track t2= new Track(); t2.setFilename("foo.txt"); t2.setFileTag("#20"); ts.getTrackSet().put(t2.getFileTag(), t2);
		Track t3= new Track(); t3.setFilename("bla.gz"); t3.setFileTag("#3"); ts.getTrackSet().put(t3.getFileTag(), t3);

		String defaultColour= (new Track()).getTitleColour();
		
		String cmdInput= "trackColour RED #\\d$";
		ts.setTrackColourForRegex(cmdInput);
		assertEquals("red", ts.getTrackSet().get("#1").getTitleColour());
		assertEquals(defaultColour, ts.getTrackSet().get("#20").getTitleColour());
		assertEquals("red", ts.getTrackSet().get("#3").getTitleColour());
		
		// Non-existant colour: Reset to default
		cmdInput= "trackColour foo .*";
		ts.setTrackColourForRegex(cmdInput);
		assertEquals(defaultColour, ts.getTrackSet().get("#1").getTitleColour());
		
		// All reset to red
		cmdInput= "trackColour red";
		ts.setTrackColourForRegex(cmdInput);
		assertEquals("red", ts.getTrackSet().get("#1").getTitleColour());
		
		// All reset to default
		cmdInput= "trackColour";
		ts.setTrackColourForRegex(cmdInput);
		assertEquals(defaultColour, ts.getTrackSet().get("#1").getTitleColour());
	}
	
	@Test
	public void canSetTrackHeight() throws InvalidCommandLineException, IOException, InvalidGenomicCoordsException{
				
		TrackSet ts= new TrackSet();
		Track t1= new Track(); t1.setFilename("foo.gz"); t1.setFileTag("#1"); ts.getTrackSet().put(t1.getFileTag(), t1);
		Track t2= new Track(); t2.setFilename("foo.txt"); t2.setFileTag("#20"); ts.getTrackSet().put(t2.getFileTag(), t2);
		Track t3= new Track(); t3.setFilename("bla.gz"); t3.setFileTag("#3"); ts.getTrackSet().put(t3.getFileTag(), t3);

		String cmdInput= "trackHeight 2 #\\d$";
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
		
		Track t1= new Track(); t1.setFilename("foo.gz"); t1.setFileTag("#1"); ts.getTrackSet().put(t1.getFileTag(), t1);
		Track t2= new Track(); t2.setFilename("foo.txt"); t2.setFileTag("#20"); ts.getTrackSet().put(t2.getFileTag(), t2);
		Track t3= new Track(); t3.setFilename("bla.gz"); t3.setFileTag("#3"); ts.getTrackSet().put(t3.getFileTag(), t3);
		
		ts.setTrackYlimitsForRegex(cmdInput);
				
		assertEquals(0, ts.getTrackSet().get("#1").getYLimitMin(), 0.001);
		assertEquals(0, ts.getTrackSet().get("#20").getYLimitMin(), 0.001);
		assertEquals(10, ts.getTrackSet().get("#1").getYLimitMax(), 0.001);
	}
	
}
