package tracks;

import static org.junit.Assert.*;

import java.io.File;
import java.io.IOException;
import java.net.MalformedURLException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.io.FileUtils;
import org.junit.Test;

import exceptions.BamIndexNotFoundException;
import exceptions.InvalidCommandLineException;
import exceptions.InvalidGenomicCoordsException;
import exceptions.InvalidRecordException;
import samTextViewer.GenomicCoords;
import samTextViewer.Utils;

public class TrackSetTest {

	@Test 
	public void canTrimGenomicCoordinatesForTrack() throws InvalidGenomicCoordsException, IOException, ClassNotFoundException, BamIndexNotFoundException, InvalidRecordException, SQLException{
		
		GenomicCoords gc= new GenomicCoords("chr7:5565052-5571960", null, null);
		
		TrackSet trackSet= new TrackSet();
		trackSet.addTrackFromSource("test_data/refSeq.hg19.bed.gz", gc, null);
		trackSet.addTrackFromSource("test_data/hg19_genes.gtf.gz", gc, null);
		
		GenomicCoords cropped= trackSet.trimCoordsForTrack(Utils.tokenize("trim refSeq.hg19", " "));
		assertEquals(5566778+1, (int)cropped.getFrom());
		assertEquals(5567378, (int)cropped.getTo());
		
		// No feature in #1: What happens if we try to trim
		gc= new GenomicCoords("chr7:5568506-5575414", null, null);
		trackSet= new TrackSet();
		trackSet.addTrackFromSource("test_data/refSeq.hg19.bed.gz", gc, null);
		trackSet.addTrackFromSource("test_data/hg19_genes.gtf.gz", gc, null);
		// No change:
		cropped= trackSet.trimCoordsForTrack(Utils.tokenize("trim refSeq.hg19", " "));
		assertEquals(5568506, (int)cropped.getFrom());
		assertEquals(5575414, (int)cropped.getTo());
		
		// Trim when feature(s) extend beyond the current window: No cropping
		gc= new GenomicCoords("chr7:5566843-5567275", null, null);
		trackSet= new TrackSet();
		trackSet.addTrackFromSource("test_data/refSeq.hg19.bed.gz", gc, null);
		trackSet.addTrackFromSource("test_data/hg19_genes.gtf.gz", gc, null);

		cropped= trackSet.trimCoordsForTrack(Utils.tokenize("trim refSeq.hg19", " "));
		assertEquals(5566843, (int)cropped.getFrom());
		assertEquals(5567275, (int)cropped.getTo());

		// Trim w/o tracks
		trackSet= new TrackSet();
		cropped= trackSet.trimCoordsForTrack(Utils.tokenize("trim refSeq.hg19", " "));
		assertEquals(null, cropped);
	}
	
	@Test
	public void canPrintFeaturesToFile() throws InvalidGenomicCoordsException, IOException, ClassNotFoundException, BamIndexNotFoundException, InvalidRecordException, SQLException, InvalidCommandLineException{
		
		// --------------------------------------------------------------------
		// Prepare coords and trackSet
		GenomicCoords gc= new GenomicCoords("chr7:5565052-5571960", null, null);
		
		TrackSet trackSet= new TrackSet();
		trackSet.addTrackFromSource("test_data/hg19_genes.gtf.gz", gc, null);
		trackSet.addTrackFromSource("test_data/hg19_genes.gtf.gz", gc, null);
		// --------------------------------------------------------------------

		// No redirection file: Nothing done
		List<String>cmdInput= Utils.tokenize("print #1 >", " ");
		trackSet.setPrintModeAndPrintFeaturesForRegex(cmdInput);

		// Invalid output: Non existent dir:
		File ff= new File("test_data/foobar/delete.me");
		ff.deleteOnExit();
		cmdInput= Utils.tokenize("print #1 > test_data/foobar/delete.me", " ");
		boolean failed= false;
		try{
			trackSet.setPrintModeAndPrintFeaturesForRegex(cmdInput);
		} catch (IOException e){
			failed= true;
		}
		assertTrue(failed);		
		
		// Now give an output file.
		ff= new File("deleteme.gtf");
		ff.deleteOnExit();
		cmdInput= Utils.tokenize("print #1 > deleteme.gtf", " ");

		trackSet.setPrintModeAndPrintFeaturesForRegex(cmdInput);
		assertTrue(ff.exists());
		assertTrue(ff.length() > 200);
		assertEquals(13, FileUtils.readLines(ff).size());
		
		// Append to file
		cmdInput= Utils.tokenize("print #1 >> deleteme.gtf", " ");
		trackSet.setPrintModeAndPrintFeaturesForRegex(cmdInput);
		assertEquals(26, FileUtils.readLines(ff).size());
		
	}
	
	@Test
	public void canDropTracks() throws InvalidGenomicCoordsException, IOException, ClassNotFoundException, BamIndexNotFoundException, InvalidRecordException, SQLException, InvalidCommandLineException{

		GenomicCoords gc= new GenomicCoords("chr7:1-100", null, null);
		
		TrackSet trackSet= new TrackSet();
		trackSet.addTrackFromSource("test_data/hg19_genes.gtf.gz", gc, null);
		trackSet.addTrackFromSource("test_data/hg18_var_sample.wig.v2.1.30.tdf", gc, null);
		trackSet.addTrackFromSource("test_data/hg18_var_sample.wig.v2.1.30.tdf", gc, null);
		
		// Test only
		ArrayList<String> cmdInput = Utils.tokenize("dropTracks -t .*", " ");
		trackSet.dropTracksWithRegex(cmdInput);
		assertEquals(3, trackSet.getTrackList().size());
		
		cmdInput = Utils.tokenize("dropTracks gtf.gz", " ");
		trackSet.dropTracksWithRegex(cmdInput);
		assertEquals(2, trackSet.getTrackList().size());
		
		// Drop all
		cmdInput = Utils.tokenize("dropTracks .*", " ");
		String msg= trackSet.dropTracksWithRegex(cmdInput);
		assertEquals(0, trackSet.getTrackList().size());
		assertTrue(msg.length() > 10);
	}
	
	@Test // Disable to save time
	public void canAddTrackFromSourcename() throws InvalidGenomicCoordsException, IOException, BamIndexNotFoundException, InvalidRecordException, ClassNotFoundException, SQLException{
		
		GenomicCoords gc= new GenomicCoords("chr7:1-100", null, null);
		
		TrackSet trackSet= new TrackSet();
		trackSet.addTrackFromSource("test_data/hg19_genes.gtf.gz", gc, null);
		trackSet.addTrackFromSource("test_data/hg18_var_sample.wig.v2.1.30.tdf", gc, null);
		assertEquals(2, trackSet.getTrackList().size());
		
		trackSet.addTrackFromSource("test_data/ds051.actb.bam", gc, null);
		assertEquals(4, trackSet.getTrackList().size());
		
	}
		
	@Test // Disable to save time
	public void canInitializeFromListOfFileNames() throws InvalidGenomicCoordsException, IOException, InvalidRecordException, ClassNotFoundException, SQLException{
		
		List<String> inputFileList= new ArrayList<String>();
		inputFileList.add("test_data/ear045.oxBS.actb.bam");
		inputFileList.add("test_data/hg19_genes.gtf.gz");
		inputFileList.add("test_data/posNeg.bedGraph.gz");
		
		GenomicCoords gc= new GenomicCoords("chr1:1-100", null, null);
		TrackSet trackSet= new TrackSet(inputFileList, gc);

		assertEquals(4, trackSet.getTrackList().size()); // MEMO: BAM files add 2 tracks. 
		
	}
	
	@Test
	public void canAddBookmarkTrack() throws InvalidGenomicCoordsException, IOException, ClassNotFoundException, InvalidRecordException, SQLException{

		List<String>cmdInput= new ArrayList<String>();
		cmdInput.add("bookmark");
		cmdInput.add("bookmark_1");
		TrackSet ts= new TrackSet();
		GenomicCoords gc= new GenomicCoords("chr1:1-100", null, null);
		
		ts.bookmark(gc, cmdInput);
		assertTrue(ts.getTrackList().size() == 1);
		
		ts.getTrackList().get(0).setNoFormat(true);

		GenomicCoords gc2 = new GenomicCoords("chr1:1-1000", null, null);
		ts.bookmark(gc2, cmdInput);
		
		TrackBookmark bm = (TrackBookmark) ts.getTrackList().get(0);
		assertEquals(2, bm.getIntervalFeatureList().size());

		bm.setPrintMode(PrintRawLine.CLIP);
		System.out.println(bm.printFeaturesToFile()); // NB: it prints twice the same gc becouse the position is nt changed
		
	}
	
	@Test
	public void canReorderTracks() throws InvalidGenomicCoordsException, IOException, InvalidCommandLineException{
		TrackSet ts= new TrackSet();
		Track t1= new Track(); t1.setFilename("foo.gz");  ts.addTrack(t1, "foo.gz");
		Track t2= new Track(); t2.setFilename("foo.txt"); ts.addTrack(t2, "foo.txt");
		Track t3= new Track(); t3.setFilename("bla.gz"); ts.addTrack(t3, "bla.gz");

		List<String> newOrder= new ArrayList<String>();
		newOrder.add("foo.gz#1");
		newOrder.add("foo.txt#2");
		newOrder.add("bla.gz#3");
		ts.orderTracks(newOrder);
		assertEquals(newOrder, new ArrayList<String>(ts.getTrackTags()));

		// Handle missing tracks in new order
		newOrder= new ArrayList<String>();
		newOrder.add("bla.gz#3");
		newOrder.add("foo.txt#2");
		ts.orderTracks(newOrder);
		assertEquals(3, ts.getTrackList().size());
		
		// Handle non existing tracks
		newOrder= new ArrayList<String>();
		ts.orderTracks(newOrder);
		assertEquals(3, ts.getTrackList().size());

		newOrder= new ArrayList<String>();
		newOrder.add("1");
		newOrder.add("1");
		newOrder.add("2");
		newOrder.add("3");
		newOrder.add("3");
		newOrder.add("foo!");
		ts.orderTracks(newOrder);
		assertEquals(3, ts.getTrackList().size());

		// Partial matches
		newOrder= new ArrayList<String>();
		newOrder.add("2");
		newOrder.add("3");
		newOrder.add("1");
		ts.orderTracks(newOrder);
		assertEquals("foo.txt#2", ts.getTrackTags().get(0));
		assertEquals("bla.gz#3", ts.getTrackTags().get(1));
		
		ts.orderTracks(new ArrayList<String>());
		assertEquals("bla.gz#3", ts.getTrackTags().get(0));
		assertEquals("foo.gz#1", ts.getTrackTags().get(1));

	}
	
	@Test
	public void canSetFilterForTrackIntervalFeature() throws InvalidCommandLineException, IOException, InvalidGenomicCoordsException, ClassNotFoundException, InvalidRecordException, SQLException{
				
		TrackSet ts= new TrackSet();
		GenomicCoords gc= new GenomicCoords("chr1:1-100", null, null);
		Track t1= new TrackIntervalFeature("test_data/hg19_genes_head.gtf", gc); ts.addTrack(t1, "x");
		Track t2= new TrackIntervalFeature("test_data/hg19_genes_head.gtf.gz", gc); ts.addTrack(t2, "x");
		Track t3= new TrackIntervalFeature("test_data/refSeq.bed", gc); ts.addTrack(t3, "x");
		
		// MEMO: Track tags:
		// [hg19_genes_head.gtf#1, hg19_genes_head.gtf.gz#2, refSeq.bed#3]
		
		String cmdInput= "grep -i exon -e intron #1"; // Set for #1...
		//String cmdInput= "filter exon intron #1"; // Set for #1...
		ts.setFilterForTrackIntervalFeature(Utils.tokenize(cmdInput, " "));

		assertEquals("exon", ts.getTrack(t1).getShowRegex());
		assertEquals("intron", ts.getTrack(t1).getHideRegex());
		assertEquals(".*", ts.getTrack(t3).getShowRegex()); // As default
		assertEquals("^$", ts.getTrack(t3).getHideRegex());

		// cmdInput= "filter exon intron #1 #3"; // Set for #1...
		cmdInput= "grep -i exon -e intron #1 #3"; // Set for #1...
		ts.setFilterForTrackIntervalFeature(Utils.tokenize(cmdInput, " "));
		assertEquals("exon", ts.getTrack(t3).getShowRegex());
		assertEquals("intron", ts.getTrack(t3).getHideRegex());

		cmdInput= "grep"; // Reset all to default
		ts.setFilterForTrackIntervalFeature(Utils.tokenize(cmdInput, " "));
		assertEquals(".*", ts.getTrack(t3).getShowRegex()); // As default
		assertEquals("^$", ts.getTrack(t3).getHideRegex());
	}


	@Test
	public void canSetBSMode() throws InvalidCommandLineException, IOException, InvalidGenomicCoordsException{

		TrackSet ts= new TrackSet();
		Track t1= new Track(); t1.setFilename("foo.bam");  ts.addTrack(t1, "foo.bam");
		Track t2= new Track(); t2.setFilename("bar.bam"); ts.addTrack(t2, "bar.bam");
		Track t3= new Track(); t3.setFilename("foo.bam"); ts.addTrack(t3, "foo.bam");
		
		String cmdInput= "bsMode foo.*#\\d";
		ts.setBisulfiteModeForRegex(Utils.tokenize(cmdInput, " "));
		assertTrue(ts.getTrack(t1).isBisulf());
		assertTrue(! ts.getTrack(t2).isBisulf());
		
		ts.setBisulfiteModeForRegex(Utils.tokenize(cmdInput, " ")); // Was set true, now becomes false
		assertTrue(! ts.getTrack(t1).isBisulf());
	}

	@Test
	public void canSetRpmForRegex() throws InvalidCommandLineException, IOException, InvalidGenomicCoordsException, ClassNotFoundException, InvalidRecordException, SQLException{
		
		TrackSet ts= new TrackSet();
		Track t1= new Track(); ts.addTrack(t1, "x");
		Track t2= new Track(); ts.addTrack(t2, "x");
		Track t3= new Track(); ts.addTrack(t3, "x");

		
		String cmdInput= "bsMode #1 #3";
		ts.setRpmForRegex(Utils.tokenize(cmdInput, " "));
		assertTrue(ts.getTrack(t1).isRpm());
		assertTrue(! ts.getTrack(t2).isRpm());
		
		ts.setRpmForRegex(Utils.tokenize(cmdInput, " ")); // Was set true, now becomes false
		assertTrue(! ts.getTrack(t1).isRpm());
	}

	@Test
	public void canEditTrackNamesForRegex() throws InvalidCommandLineException, IOException, InvalidGenomicCoordsException, SQLException, InvalidRecordException, ClassNotFoundException{
				
		TrackSet ts= new TrackSet();
		Track t1= new Track(); ts.addTrack(t1, "foo.gff");
		Track t2= new Track(); ts.addTrack(t2, "foo.bed");
		Track t3= new Track(); ts.addTrack(t3, "baz.narrowPeak");
		
		ts.editNamesForRegex(Utils.tokenize("editNames foo FOO", " "));
		
		assertTrue(ts.getTrackList().get(0).getTrackTag().startsWith("FOO.gff"));
		assertTrue(ts.getTrackList().get(1).getTrackTag().startsWith("FOO.bed"));
		assertTrue(ts.getTrackList().get(2).getTrackTag().startsWith("baz"));

		// Replace with nothing: ""
		ts.editNamesForRegex(Utils.tokenize("editNames .bed \"\"", " "));
		assertTrue(ts.getTrackList().get(1).getTrackTag().startsWith("FOO#"));
	}
	
	@Test
	public void canSetFilterFlagForRegex() throws InvalidCommandLineException, IOException, InvalidGenomicCoordsException, SQLException, InvalidRecordException, ClassNotFoundException{
				
		TrackSet ts= new TrackSet();
		Track t1= new Track(); ts.addTrack(t1, "x");
		Track t2= new Track(); ts.addTrack(t2, "x");
		Track t3= new Track(); ts.addTrack(t3, "x");

		// String cmdInput= "-F 1024 #1 #3";
		// Reset all three filters
		String cmdInput= "samtools -q 10 -F 1024 -f 16 #1 #3";
		ts.setSamFilterForRegex(Utils.tokenize(cmdInput, " "));
		assertEquals(1024+4, ts.getTrack(t1).get_F_flag());
		assertEquals(16, ts.getTrack(t1).get_f_flag());
		assertEquals(10, ts.getTrack(t1).getMapq());
		// Not changed
		assertEquals(4, ts.getTrack(t2).get_F_flag());
		assertEquals(0, ts.getTrack(t2).getMapq());
		// assertEquals(4, ts.getTrack(t2).get_F_flag());
		
		// Set one filter for all tracks, the others return to zero:
		cmdInput= "samtools -f 16";
		ts.setSamFilterForRegex(Utils.tokenize(cmdInput, " "));
		assertEquals(16, ts.getTrack(t1).get_f_flag());
		assertEquals(16, ts.getTrack(t2).get_f_flag());
		assertEquals(16, ts.getTrack(t3).get_f_flag());

		assertEquals(4, ts.getTrack(t1).get_F_flag()); // Set to 4 (unmapped)
		assertEquals(4, ts.getTrack(t2).get_F_flag());
		assertEquals(4, ts.getTrack(t3).get_F_flag());

		assertEquals(0, ts.getTrack(t1).getMapq());
		assertEquals(0, ts.getTrack(t2).getMapq());
		assertEquals(0, ts.getTrack(t3).getMapq());
		
	}
	
	@Test
	public void canSetFeatureDisplayModeForRegex() throws InvalidCommandLineException, IOException, InvalidGenomicCoordsException{
				
		TrackSet ts= new TrackSet();
		Track t1= new Track(); ts.addTrack(t1, "x");
		Track t2= new Track(); ts.addTrack(t2, "x");
		Track t3= new Track(); ts.addTrack(t3, "x");

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
		Track t1= new Track(); ts.addTrack(t1, "x");
		Track t2= new Track(); ts.addTrack(t2, "x");
		Track t3= new Track(); ts.addTrack(t3, "x");

		ts.setPrintModeAndPrintFeaturesForRegex(Utils.tokenize("print #1 #3", " "));
		assertEquals(PrintRawLine.CLIP, t1.getPrintMode());
		assertEquals(PrintRawLine.CLIP, t3.getPrintMode());
		
		ts.setPrintModeAndPrintFeaturesForRegex(Utils.tokenize("print #1", " "));
		assertEquals(PrintRawLine.OFF, t1.getPrintMode());
		
		ts.setPrintModeAndPrintFeaturesForRegex(Utils.tokenize("print #1", " "));
		assertEquals(PrintRawLine.CLIP, t1.getPrintMode());
		
		ts.setPrintModeAndPrintFeaturesForRegex(Utils.tokenize("print -full #1", " "));
		assertEquals(PrintRawLine.FULL, t1.getPrintMode());
		
		ts.setPrintModeAndPrintFeaturesForRegex(Utils.tokenize("print -full #1", " "));
		assertEquals(PrintRawLine.OFF, t1.getPrintMode());
	}
	
	
	@Test
	public void canSetTrackColour() throws InvalidCommandLineException, IOException, InvalidGenomicCoordsException{
				
		TrackSet ts= new TrackSet();
		Track t1= new Track(); t1.setFilename("foo.gz");  ts.addTrack(t1, "foo.gz");
		Track t2= new Track(); t2.setFilename("foo.txt"); ts.addTrack(t2, "foo.txt");
		Track t3= new Track(); t3.setFilename("bla.gz"); ts.addTrack(t3, "bla.gz");

		String defaultColour= (new Track()).getTitleColour();
			
		String cmdInput= "trackColour RED gz#\\d$";
		ts.setTrackColourForRegex(Utils.tokenize(cmdInput, " "));
		assertEquals("red", ts.getTrack(t1).getTitleColour());
		assertEquals(defaultColour, ts.getTrack(t2).getTitleColour());
		assertEquals("red", ts.getTrack(t3).getTitleColour());
		
		// Non-existant colour: Throw exception
		cmdInput= "trackColour foo .*";
		boolean passed= false;
		try{
			ts.setTrackColourForRegex(Utils.tokenize(cmdInput, " "));
		} catch(InvalidCommandLineException e){
			passed= true;
		}
		assertTrue(passed);
		
		// All reset to red
		cmdInput= "trackColour red";
		ts.setTrackColourForRegex(Utils.tokenize(cmdInput, " "));
		assertEquals("red", ts.getTrack(t1).getTitleColour());
		
		// All reset to default
		cmdInput= "trackColour";
		ts.setTrackColourForRegex(Utils.tokenize(cmdInput, " "));
		assertEquals(defaultColour, ts.getTrack(t1).getTitleColour());
		
		cmdInput= "trackColour red #1 #3 #1";
		ts.setTrackColourForRegex(Utils.tokenize(cmdInput, " "));
		
		assertEquals("red", ts.getTrack(t1).getTitleColour());
		assertEquals(defaultColour, ts.getTrack(t2).getTitleColour());
		assertEquals("red", ts.getTrack(t3).getTitleColour());
	}
	
	@Test
	public void canSetTrackHeight() throws InvalidCommandLineException, IOException, InvalidGenomicCoordsException, ClassNotFoundException, InvalidRecordException, SQLException{
				
		TrackSet ts= new TrackSet();
		Track t1= new Track(); ts.addTrack( t1, "x");
		Track t2= new Track(); ts.addTrack(t2, "x");
		Track t3= new Track(); ts.addTrack(t3, "x");

		String cmdInput= "trackHeight 2 #1 #3";
		ts.setTrackHeightForRegex(Utils.tokenize(cmdInput, " "));

		assertEquals(2, ts.getTrack(t1).getyMaxLines());
		assertEquals(10, ts.getTrack(t2).getyMaxLines()); 
		assertEquals(2, ts.getTrack(t3).getyMaxLines());
		
		cmdInput= "trackHeight 99"; // Defsult regex: All tracks captured 
		ts.setTrackHeightForRegex(Utils.tokenize(cmdInput, " "));
		assertEquals(99, ts.getTrack(t1).getyMaxLines());
		assertEquals(99, ts.getTrack(t2).getyMaxLines()); 
		assertEquals(99, ts.getTrack(t3).getyMaxLines());
		
	}
	
	@Test
	public void canSetYlimits() throws InvalidCommandLineException, MalformedURLException, ClassNotFoundException, IOException, InvalidGenomicCoordsException, InvalidRecordException, SQLException{
				
		TrackSet ts= new TrackSet();
		Track t1= new Track(); ts.addTrack(t1, "x");
		Track t2= new Track(); ts.addTrack(t2, "x");
		Track t3= new Track(); ts.addTrack(t3, "x");
		
		String cmdInput= "ylim 10 20 #1 #2";
		ts.setTrackYlimitsForRegex(Utils.tokenize(cmdInput, " "));
		assertEquals(10, ts.getTrack(t1).getYLimitMin(), 0.001);
		assertEquals(20, ts.getTrack(t1).getYLimitMax(), 0.001);
		assertEquals(10, ts.getTrack(t2).getYLimitMin(), 0.001);
		assertEquals(20, ts.getTrack(t2).getYLimitMax(), 0.001);
		
		cmdInput= "ylim 90 99";
		ts.setTrackYlimitsForRegex(Utils.tokenize(cmdInput, " "));
		assertEquals(90, ts.getTrack(t1).getYLimitMin(), 0.001);
		assertEquals(99, ts.getTrack(t1).getYLimitMax(), 0.001);
		assertEquals(90, ts.getTrack(t2).getYLimitMin(), 0.001);
		assertEquals(99, ts.getTrack(t2).getYLimitMax(), 0.001);
		assertEquals(90, ts.getTrack(t3).getYLimitMin(), 0.001);
		assertEquals(99, ts.getTrack(t3).getYLimitMax(), 0.001);

	}
	
	@Test
	public void canPrintTrackInfo(){
		
		TrackSet ts= new TrackSet();
		assertEquals("", ts.showTrackInfo());

		Track t1= new Track(); t1.setFilename("/path/to/foo.gz"); ts.addTrack(t1, "foo.gz");
		Track t2= new Track(); t2.setFilename("/path/to/foo.vcf"); ts.addTrack(t2, "foo.vcf");
		Track t3= new Track(); t3.setFilename("/path/to/bla.gz"); ts.addTrack(t3, "bla.gz");

		assertTrue(ts.showTrackInfo().contains("foo.gz"));
		assertTrue(ts.showTrackInfo().contains("/path/to/foo.gz"));
		assertTrue(ts.showTrackInfo().contains("BED"));

	}

}
