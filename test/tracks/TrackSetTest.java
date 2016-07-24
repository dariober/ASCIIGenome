package tracks;

import static org.junit.Assert.*;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.junit.Test;

import exceptions.InvalidCommandLineException;
import exceptions.InvalidGenomicCoordsException;
import samTextViewer.GenomicCoords;
import samTextViewer.Utils;

public class TrackSetTest {

	@Test
	public void canAddBookmarkTrack() throws InvalidGenomicCoordsException, IOException{

		TrackSet ts= new TrackSet();
		GenomicCoords gc= new GenomicCoords("chr1", 1, 100, null, 100, null);
		ts.addBookmark_IN_PREP(gc, "bookmark_n1");
		assertTrue(ts.getTrackSet().containsKey(TrackSet.BOOKMARK_TAG));
		TrackIntervalFeature tif= (TrackIntervalFeature) ts.getTrackSet().get(TrackSet.BOOKMARK_TAG);
		assertEquals(1, tif.intervalFeatureSet.getIntervalMap().get("chr1").get(0).getFrom());
		
		GenomicCoords gc2= new GenomicCoords("chr2", 90, 100, null, 100, null);
		ts.addBookmark_IN_PREP(gc2, "bookmark_n2");
		assertEquals(90, tif.intervalFeatureSet.getIntervalMap().get("chr2").get(0).getFrom());
		
		GenomicCoords gc3= new GenomicCoords("chr2", 2, 100, null, 100, null);
		ts.addBookmark_IN_PREP(gc3, "bookmark_n3");
		assertEquals(2, tif.intervalFeatureSet.getIntervalMap().get("chr2").get(0).getFrom());
				
		// System.out.println(tif.intervalFeatureSet.getIntervalMap().get("chr2"));
		// System.out.println(tif.intervalFeatureSet.getIntervalMap().get("chr2").get(0).getRaw());
	}
	
	@Test
	public void canReorderTracks() throws InvalidGenomicCoordsException, IOException, InvalidCommandLineException{
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
		
		String cmdInput= "visible exon intron #1"; // Set for #1...
		ts.setVisibilityForTrackIntervalFeature(Utils.tokenize(cmdInput, " "));

		assertEquals("exon", ts.getTrackSet().get("#10").getShowRegex());
		assertEquals("intron", ts.getTrackSet().get("#11").getHideRegex());
		assertEquals(".*", ts.getTrackSet().get("#30").getShowRegex()); // As default
		assertEquals("^$", ts.getTrackSet().get("#30").getHideRegex());

		cmdInput= "visible exon intron #10 #3"; // Set for #1...
		ts.setVisibilityForTrackIntervalFeature(Utils.tokenize(cmdInput, " "));
		assertEquals("exon", ts.getTrackSet().get("#30").getShowRegex()); // As default
		assertEquals("intron", ts.getTrackSet().get("#30").getHideRegex());

		cmdInput= "visible"; // Set for #1...
		ts.setVisibilityForTrackIntervalFeature(Utils.tokenize(cmdInput, " "));
		assertEquals(".*", ts.getTrackSet().get("#30").getShowRegex()); // As default
		assertEquals("^$", ts.getTrackSet().get("#30").getHideRegex());
	}


	@Test
	public void canSetBSMode() throws InvalidCommandLineException, IOException, InvalidGenomicCoordsException{
				
		TrackSet ts= new TrackSet();
		Track t1= new Track(); t1.setFilename("foo.bam"); t1.setFileTag("#1"); ts.getTrackSet().put(t1.getFileTag(), t1);
		Track t2= new Track(); t2.setFilename("foo.bam"); t2.setFileTag("#20"); ts.getTrackSet().put(t2.getFileTag(), t2);
		Track t3= new Track(); t3.setFilename("bla.gz"); t3.setFileTag("#3"); ts.getTrackSet().put(t3.getFileTag(), t3);

		String cmdInput= "bsMode #\\d$";
		ts.setBisulfiteModeForRegex(Utils.tokenize(cmdInput, " "));
		assertTrue(ts.getTrackSet().get("#1").isBisulf());
		assertTrue(! ts.getTrackSet().get("#20").isBisulf());
		
		ts.setBisulfiteModeForRegex(Utils.tokenize(cmdInput, " ")); // Was set true, now becomes false
		assertTrue(! ts.getTrackSet().get("#1").isBisulf());
	}

	@Test
	public void canSetRpmForRegex() throws InvalidCommandLineException, IOException, InvalidGenomicCoordsException{
				
		TrackSet ts= new TrackSet();
		Track t1= new Track(); t1.setFilename("foo.bam"); t1.setFileTag("#1"); ts.getTrackSet().put(t1.getFileTag(), t1);
		Track t2= new Track(); t2.setFilename("foo.bam"); t2.setFileTag("#20"); ts.getTrackSet().put(t2.getFileTag(), t2);
		Track t3= new Track(); t3.setFilename("bla.txt"); t3.setFileTag("#3"); ts.getTrackSet().put(t3.getFileTag(), t3);

		String cmdInput= "rpm #1 #3";
		ts.setRpmForRegex(Utils.tokenize(cmdInput, " "));
		assertTrue(ts.getTrackSet().get("#1").isRpm());
		assertTrue(! ts.getTrackSet().get("#20").isRpm());
		
		ts.setRpmForRegex(Utils.tokenize(cmdInput, " ")); // Was set true, now becomes false
		assertTrue(! ts.getTrackSet().get("#1").isRpm());
	}

	@Test
	public void canSetFilterFlagForRegex() throws InvalidCommandLineException, IOException, InvalidGenomicCoordsException{
				
		TrackSet ts= new TrackSet();
		Track t1= new Track(); t1.setFilename("foo.bam"); t1.setFileTag("#1"); ts.getTrackSet().put(t1.getFileTag(), t1);
		Track t2= new Track(); t2.setFilename("foo.bam"); t2.setFileTag("#20"); ts.getTrackSet().put(t2.getFileTag(), t2);
		Track t3= new Track(); t3.setFilename("bla.txt"); t3.setFileTag("#3"); ts.getTrackSet().put(t3.getFileTag(), t3);

		String cmdInput= "-F 1024 #1 #3";
		ts.setFilterFlagForRegex(Utils.tokenize(cmdInput, " "));
		assertEquals(1024+4, ts.getTrackSet().get("#1").get_F_flag());
		assertEquals(4, ts.getTrackSet().get("#20").get_F_flag());
		
		cmdInput= "-F 16";
		ts.setFilterFlagForRegex(Utils.tokenize(cmdInput, " "));
		cmdInput= "-f 1024 #1 #3";
		ts.setFilterFlagForRegex(Utils.tokenize(cmdInput, " "));		
		cmdInput= "mapq 30 #1 #3";
		ts.setFilterFlagForRegex(Utils.tokenize(cmdInput, " "));		
		
		assertEquals(16+4, ts.getTrackSet().get("#1").get_F_flag());
		assertEquals(16+4, ts.getTrackSet().get("#20").get_F_flag());
		assertEquals(16+4, ts.getTrackSet().get("#3").get_F_flag());

		assertEquals(1024, ts.getTrackSet().get("#1").get_f_flag());
		assertEquals(0, ts.getTrackSet().get("#20").get_f_flag());
		assertEquals(30, ts.getTrackSet().get("#3").getMapq());
		
	}
	
	@Test
	public void canSetFeatureDisplayModeForRegex() throws InvalidCommandLineException, IOException, InvalidGenomicCoordsException{
				
		TrackSet ts= new TrackSet();
		Track t1= new Track(); t1.setFilename("foo.gz"); t1.setFileTag("#1"); ts.getTrackSet().put(t1.getFileTag(), t1);
		Track t2= new Track(); t2.setFilename("foo.txt"); t2.setFileTag("#20"); ts.getTrackSet().put(t2.getFileTag(), t2);
		Track t3= new Track(); t3.setFilename("bla.gz"); t3.setFileTag("#3"); ts.getTrackSet().put(t3.getFileTag(), t3);

		ts.setFeatureDisplayModeForRegex(Utils.tokenize("squash #1 #3", " "));
		assertEquals(FeatureDisplayMode.SQUASHED, t1.getFeatureDisplayMode());
		assertEquals(FeatureDisplayMode.SQUASHED, t3.getFeatureDisplayMode());
		
		ts.setFeatureDisplayModeForRegex(Utils.tokenize("squash #1 #3", " "));
		assertEquals(FeatureDisplayMode.EXPANDED, t1.getFeatureDisplayMode());
		assertEquals(FeatureDisplayMode.EXPANDED, t3.getFeatureDisplayMode());

		ts.setFeatureDisplayModeForRegex(Utils.tokenize("merge #1 #3", " "));
		assertEquals(FeatureDisplayMode.MERGED, t1.getFeatureDisplayMode());
		assertEquals(FeatureDisplayMode.MERGED, t3.getFeatureDisplayMode());

		ts.setFeatureDisplayModeForRegex(Utils.tokenize("squash", " "));
		assertEquals(FeatureDisplayMode.SQUASHED, t1.getFeatureDisplayMode());
		assertEquals(FeatureDisplayMode.SQUASHED, t3.getFeatureDisplayMode());
		
	}
	
	
	@Test
	public void canSetPrintMode() throws InvalidCommandLineException, IOException, InvalidGenomicCoordsException{
				
		TrackSet ts= new TrackSet();
		Track t1= new Track(); t1.setFilename("foo.gz"); t1.setFileTag("#1"); ts.getTrackSet().put(t1.getFileTag(), t1);
		Track t2= new Track(); t2.setFilename("foo.txt"); t2.setFileTag("#20"); ts.getTrackSet().put(t2.getFileTag(), t2);
		Track t3= new Track(); t3.setFilename("bla.gz"); t3.setFileTag("#3"); ts.getTrackSet().put(t3.getFileTag(), t3);

		ts.setPrintModeForRegex(Utils.tokenize("print #1 #3", " "));
		assertEquals(PrintRawLine.CLIP, t1.getPrintMode());
		assertEquals(PrintRawLine.CLIP, t3.getPrintMode());
		
		ts.setPrintModeForRegex(Utils.tokenize("print #1", " "));
		assertEquals(PrintRawLine.OFF, t1.getPrintMode());
		
		ts.setPrintModeForRegex(Utils.tokenize("print #1", " "));
		assertEquals(PrintRawLine.CLIP, t1.getPrintMode());
		
		ts.setPrintModeForRegex(Utils.tokenize("printFull #1", " "));
		assertEquals(PrintRawLine.FULL, t1.getPrintMode());
		
		ts.setPrintModeForRegex(Utils.tokenize("printFull #1", " "));
		assertEquals(PrintRawLine.OFF, t1.getPrintMode());
	}
	
	
	@Test
	public void canSetTrackColour() throws InvalidCommandLineException, IOException, InvalidGenomicCoordsException{
				
		TrackSet ts= new TrackSet();
		Track t1= new Track(); t1.setFilename("foo.gz"); t1.setFileTag("#1"); ts.getTrackSet().put(t1.getFileTag(), t1);
		Track t2= new Track(); t2.setFilename("foo.txt"); t2.setFileTag("#20"); ts.getTrackSet().put(t2.getFileTag(), t2);
		Track t3= new Track(); t3.setFilename("bla.gz"); t3.setFileTag("#3"); ts.getTrackSet().put(t3.getFileTag(), t3);

		String defaultColour= (new Track()).getTitleColour();
		
		String cmdInput= "trackColour RED #\\d$";
		ts.setTrackColourForRegex(Utils.tokenize(cmdInput, " "));
		assertEquals("red", ts.getTrackSet().get("#1").getTitleColour());
		assertEquals(defaultColour, ts.getTrackSet().get("#20").getTitleColour());
		assertEquals("red", ts.getTrackSet().get("#3").getTitleColour());
		
		// Non-existant colour: Reset to default
		cmdInput= "trackColour foo .*";
		ts.setTrackColourForRegex(Utils.tokenize(cmdInput, " "));
		assertEquals(defaultColour, ts.getTrackSet().get("#1").getTitleColour());
		
		// All reset to red
		cmdInput= "trackColour red";
		ts.setTrackColourForRegex(Utils.tokenize(cmdInput, " "));
		assertEquals("red", ts.getTrackSet().get("#1").getTitleColour());
		
		// All reset to default
		cmdInput= "trackColour";
		ts.setTrackColourForRegex(Utils.tokenize(cmdInput, " "));
		assertEquals(defaultColour, ts.getTrackSet().get("#1").getTitleColour());
		
		// Multiple regexes
		ts= new TrackSet();
		t1= new Track(); t1.setFilename("foo.gz"); t1.setFileTag("#1"); ts.getTrackSet().put(t1.getFileTag(), t1);
		t2= new Track(); t2.setFilename("foo.txt"); t2.setFileTag("#20"); ts.getTrackSet().put(t2.getFileTag(), t2);
		t3= new Track(); t3.setFilename("bla.gz"); t3.setFileTag("#3"); ts.getTrackSet().put(t3.getFileTag(), t3);
		
		cmdInput= "trackColour red #1 #3 #1";
		ts.setTrackColourForRegex(Utils.tokenize(cmdInput, " "));
		assertEquals("red", ts.getTrackSet().get("#1").getTitleColour());
		assertEquals(defaultColour, ts.getTrackSet().get("#20").getTitleColour());
		assertEquals("red", ts.getTrackSet().get("#3").getTitleColour());
	}
	
	@Test
	public void canSetTrackHeight() throws InvalidCommandLineException, IOException, InvalidGenomicCoordsException{
				
		TrackSet ts= new TrackSet();
		Track t1= new Track(); t1.setFilename("foo.gz"); t1.setFileTag("#1"); ts.getTrackSet().put(t1.getFileTag(), t1);
		Track t2= new Track(); t2.setFilename("foo.txt"); t2.setFileTag("#20"); ts.getTrackSet().put(t2.getFileTag(), t2);
		Track t3= new Track(); t3.setFilename("bla.gz"); t3.setFileTag("#3"); ts.getTrackSet().put(t3.getFileTag(), t3);

		String cmdInput= "trackHeight 2 #1 #3";
		ts.setTrackHeightForRegex(Utils.tokenize(cmdInput, " "));
				
		assertEquals(2, ts.getTrackSet().get("#1").getyMaxLines());
		assertEquals(10, ts.getTrackSet().get("#20").getyMaxLines()); 
		assertEquals(2, ts.getTrackSet().get("#3").getyMaxLines());
		
		cmdInput= "trackHeight 99"; // Defsult regex: All tracks captured 
		ts.setTrackHeightForRegex(Utils.tokenize(cmdInput, " "));
		assertEquals(99, ts.getTrackSet().get("#1").getyMaxLines());
		assertEquals(99, ts.getTrackSet().get("#20").getyMaxLines()); 
		assertEquals(99, ts.getTrackSet().get("#3").getyMaxLines());
		
	}
	
	@Test
	public void canSetYlimits() throws InvalidCommandLineException{
				
		TrackSet ts= new TrackSet();	
		Track t1= new Track(); t1.setFilename("foo.gz"); t1.setFileTag("#1"); ts.getTrackSet().put(t1.getFileTag(), t1);
		Track t2= new Track(); t2.setFilename("foo.txt"); t2.setFileTag("#20"); ts.getTrackSet().put(t2.getFileTag(), t2);
		Track t3= new Track(); t3.setFilename("bla.gz"); t3.setFileTag("#3"); ts.getTrackSet().put(t3.getFileTag(), t3);
		
		String cmdInput= "ylim 10 20 #1 #2";
		ts.setTrackYlimitsForRegex(Utils.tokenize(cmdInput, " "));
		assertEquals(10, ts.getTrackSet().get("#1").getYLimitMin(), 0.001);
		assertEquals(20, ts.getTrackSet().get("#1").getYLimitMax(), 0.001);
		assertEquals(10, ts.getTrackSet().get("#20").getYLimitMin(), 0.001);
		assertEquals(20, ts.getTrackSet().get("#20").getYLimitMax(), 0.001);
		
		cmdInput= "ylim 90 99";
		ts.setTrackYlimitsForRegex(Utils.tokenize(cmdInput, " "));
		assertEquals(90, ts.getTrackSet().get("#1").getYLimitMin(), 0.001);
		assertEquals(99, ts.getTrackSet().get("#1").getYLimitMax(), 0.001);
		assertEquals(90, ts.getTrackSet().get("#20").getYLimitMin(), 0.001);
		assertEquals(99, ts.getTrackSet().get("#20").getYLimitMax(), 0.001);
		assertEquals(90, ts.getTrackSet().get("#3").getYLimitMin(), 0.001);
		assertEquals(99, ts.getTrackSet().get("#3").getYLimitMax(), 0.001);

	}
	
	@Test
	public void canPrintTrackInfo(){
		
		TrackSet ts= new TrackSet();
		assertEquals("", ts.showTrackInfo());
		
		Track t1= new Track(); t1.setFilename("/path/to/foo.gz"); t1.setFileTag("foo#1"); ts.getTrackSet().put(t1.getFileTag(), t1);
		Track t2= new Track(); t2.setFilename("/path/to/foo.vcf"); t2.setFileTag("bar#20"); ts.getTrackSet().put(t2.getFileTag(), t2);
		Track t3= new Track(); t3.setFilename("/path/to/bla.gz"); t3.setFileTag("baz#3"); ts.getTrackSet().put(t3.getFileTag(), t3);

		assertTrue(ts.showTrackInfo().startsWith("foo"));
		assertTrue(ts.showTrackInfo().endsWith("BED"));
		
	}
	
}
