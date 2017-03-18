package tracks;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.io.IOException;
import java.net.MalformedURLException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.io.FileUtils;
import org.junit.Test;

import exceptions.BamIndexNotFoundException;
import exceptions.InvalidColourException;
import exceptions.InvalidCommandLineException;
import exceptions.InvalidGenomicCoordsException;
import exceptions.InvalidRecordException;
import samTextViewer.GenomicCoords;
import samTextViewer.Utils;

public class TrackSetTest {

	@Test
	public void canConstructTrackSetFromURL() throws InvalidGenomicCoordsException, IOException, ClassNotFoundException, BamIndexNotFoundException, InvalidRecordException, SQLException{

		GenomicCoords gc= new GenomicCoords("chr7:5565052-5571960", null, null);
		
		TrackSet trackSet= new TrackSet();
		trackSet.addTrackFromSource("http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeCaltechRnaSeq/wgEncodeCaltechRnaSeqGm12878R2x75Il400SigRep2V2.bigWig", gc, null);
		trackSet.addTrackFromSource("http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeCaltechRnaSeq/wgEncodeCaltechRnaSeqGm12878R2x75Il400JunctionsRep2V3.bigBed", gc, null);
		trackSet.addTrackFromSource("http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeCaltechRnaSeq/wgEncodeCaltechRnaSeqGm12878R2x75Il400GeneGencV3cRep2V3.gtf.gz", gc, null);
		trackSet.addTrackFromSource("http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeCaltechRnaSeq/wgEncodeCaltechRnaSeqGm12878R2x75Il400SigRep2V2.bigWig", gc, null);
		
	}
	
	@Test
	public void canGetListOfOpenedFiles() throws ClassNotFoundException, IOException, BamIndexNotFoundException, InvalidGenomicCoordsException, InvalidRecordException, SQLException{

		GenomicCoords gc= new GenomicCoords("chr7:5565052-5571960", null, null);
		
		TrackSet trackSet= new TrackSet();
		trackSet.addTrackFromSource("test_data/refSeq.hg19.bed.gz", gc, null);
		trackSet.addTrackFromSource("test_data/hg19_genes.gtf.gz", gc, null);
		trackSet.addTrackFromSource("test_data/refSeq.hg19.bed.gz", gc, null);
		trackSet.addTrackFromSource("ftp://ftp.ensembl.org/pub/release-86/gff3/homo_sapiens/Homo_sapiens.GRCh38.86.chromosome.18.gff3.gz", gc, null);

		assertEquals(3, trackSet.getOpenedFiles().size());
		assertTrue(trackSet.getOpenedFiles().iterator().next().length() > 30); // Check we are getting full path
	}
	
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
	public void canPrintFeaturesToFile() throws InvalidGenomicCoordsException, IOException, ClassNotFoundException, BamIndexNotFoundException, InvalidRecordException, SQLException, InvalidCommandLineException, InvalidColourException{
		
		// --------------------------------------------------------------------
		// Prepare coords and trackSet
		GenomicCoords gc= new GenomicCoords("chr7:5565052-5571960", null, null);
		
		TrackSet trackSet= new TrackSet();
		trackSet.addTrackFromSource("test_data/hg19_genes.gtf.gz", gc, null);
		trackSet.addTrackFromSource("test_data/hg19_genes.gtf.gz", gc, null);
		// --------------------------------------------------------------------

		// No redirection file: Nothing done
		List<String>cmdInput= Utils.tokenize("print #1 >", " ");
		boolean failed= false;
		try{
			trackSet.setPrintModeAndPrintFeaturesForRegex(cmdInput);
		} catch (InvalidCommandLineException e){
			failed= true;
		}
		assertTrue(failed);		

		// Invalid output: Non existent dir:
		File ff= new File("test_data/foobar/delete.me");
		ff.deleteOnExit();
		cmdInput= Utils.tokenize("print #1 > test_data/foobar/delete.me", " ");
		failed= false;
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
	public void canAddBookmarkTrack() throws InvalidGenomicCoordsException, IOException, ClassNotFoundException, InvalidRecordException, SQLException, InvalidColourException{

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
	public void canSetAwkForTrackIntervalFeature() throws InvalidCommandLineException, IOException, InvalidGenomicCoordsException, ClassNotFoundException, InvalidRecordException, SQLException{
				
		TrackSet ts= new TrackSet();
		GenomicCoords gc= new GenomicCoords("chr1:1-100", null, null);
		Track t1= new TrackIntervalFeature("test_data/hg19_genes_head.gtf", gc); ts.addTrack(t1, "x");
		Track t2= new TrackIntervalFeature("test_data/hg19_genes_head.gtf.gz", gc); ts.addTrack(t2, "x");
		Track t3= new TrackIntervalFeature("test_data/refSeq.bed", gc); ts.addTrack(t3, "x");
		
		// Set for one track
		String cmdInput= "awk   '$3 == \"exon\"' #1"; 		
		ts.setAwkForTrackIntervalFeature(Utils.tokenize(cmdInput, " "));
		assertEquals("-F '\\t' '$3 == \"exon\"'", ts.getTrack(t1).getAwk());
		assertEquals("", ts.getTrack(t3).getAwk()); // As default
		
		// Use custom delim, some tracks
		cmdInput= "awk -F _ '$3 == 10' #1 #3"; 	
		ts.setAwkForTrackIntervalFeature(Utils.tokenize(cmdInput, " "));
		assertEquals("-F _ '$3 == 10'", ts.getTrack(t1).getAwk());
		assertEquals("-F _ '$3 == 10'", ts.getTrack(t3).getAwk());
			
		// Use custom delim: All tracks
		cmdInput= "awk -v FOO=foo -F _ '$3 == 20'"; 	
		ts.setAwkForTrackIntervalFeature(Utils.tokenize(cmdInput, " "));
		assertEquals("-v FOO=foo -F _ '$3 == 20'", ts.getTrack(t1).getAwk());
		assertEquals("-v FOO=foo -F _ '$3 == 20'", ts.getTrack(t3).getAwk());
		
		// Turn off one track
		cmdInput= "awk -off #2";
		ts.setAwkForTrackIntervalFeature(Utils.tokenize(cmdInput, " "));
		assertEquals("", ts.getTrack(t2).getAwk());

		// Turn off all tracks
		cmdInput= "awk";
		ts.setAwkForTrackIntervalFeature(Utils.tokenize(cmdInput, " "));
		assertEquals("", ts.getTrack(t1).getAwk());
		assertEquals("", ts.getTrack(t2).getAwk());
		assertEquals("", ts.getTrack(t3).getAwk());
	}
	
	@Test
	public void canSetBSMode() throws InvalidCommandLineException, IOException, InvalidGenomicCoordsException{

		TrackSet ts= new TrackSet();
		Track t1= new Track(); t1.setFilename("foo.bam");  ts.addTrack(t1, "foo.bam");
		Track t2= new Track(); t2.setFilename("bar.bam"); ts.addTrack(t2, "bar.bam");
		Track t3= new Track(); t3.setFilename("foo.bam"); ts.addTrack(t3, "foo.bam");
		
		String cmdInput= "BSseq foo.*#\\d";
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

		
		String cmdInput= "rpm #1 #3";
		ts.setRpmForRegex(Utils.tokenize(cmdInput, " "));
		assertTrue(ts.getTrack(t1).isRpm());
		assertTrue(! ts.getTrack(t2).isRpm());
		
		ts.setRpmForRegex(Utils.tokenize(cmdInput, " ")); // Was set true, now becomes false
		assertTrue(! ts.getTrack(t1).isRpm());
		
		ts.setRpmForRegex(Utils.tokenize("rpm -on", " "));
		for(Track tr : ts.getTrackList()){
			assertTrue(tr.isRpm());
		}
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

		ts.setFeatureDisplayModeForRegex(Utils.tokenize("featureDisplayMode #1 #3", " "));
		assertEquals(FeatureDisplayMode.COLLAPSED, t1.getFeatureDisplayMode());
		assertEquals(FeatureDisplayMode.COLLAPSED, t3.getFeatureDisplayMode());

		// Toggle back to expanded
		ts.setFeatureDisplayModeForRegex(Utils.tokenize("featureDisplayMode #1 #3", " "));
		assertEquals(FeatureDisplayMode.EXPANDED, t1.getFeatureDisplayMode());
		assertEquals(FeatureDisplayMode.EXPANDED, t3.getFeatureDisplayMode());
		
		ts.setFeatureDisplayModeForRegex(Utils.tokenize("featureDisplayMode -collapsed", " "));
		assertEquals(FeatureDisplayMode.COLLAPSED, t1.getFeatureDisplayMode());
		assertEquals(FeatureDisplayMode.COLLAPSED, t2.getFeatureDisplayMode());
		assertEquals(FeatureDisplayMode.COLLAPSED, t3.getFeatureDisplayMode());

	}
	
	@Test
	public void canSetPrintMode() throws InvalidCommandLineException, IOException, InvalidGenomicCoordsException, InvalidColourException{
				
		TrackSet ts= new TrackSet();
		Track t1= new Track(); ts.addTrack(t1, "x");
		Track t2= new Track(); ts.addTrack(t2, "x");
		Track t3= new Track(); ts.addTrack(t3, "x");

		ts.setPrintModeAndPrintFeaturesForRegex(Utils.tokenize("print #1 #3", " "));
		assertEquals(PrintRawLine.CLIP, t1.getPrintMode());
		assertEquals(PrintRawLine.CLIP, t3.getPrintMode());

		ts.setPrintModeAndPrintFeaturesForRegex(Utils.tokenize("print -off #1", " "));
		assertEquals(PrintRawLine.OFF, t1.getPrintMode());
		assertEquals(PrintRawLine.CLIP, t3.getPrintMode());

		ts.setPrintModeAndPrintFeaturesForRegex(Utils.tokenize("print -full", " "));
		assertEquals(PrintRawLine.FULL, t1.getPrintMode());
		assertEquals(PrintRawLine.FULL, t2.getPrintMode());
		assertEquals(PrintRawLine.FULL, t3.getPrintMode());
		
		ts.setPrintModeAndPrintFeaturesForRegex(Utils.tokenize("print #1", " "));
		assertEquals(PrintRawLine.OFF, t1.getPrintMode());
		
		ts.setPrintModeAndPrintFeaturesForRegex(Utils.tokenize("print #1", " "));
		assertEquals(PrintRawLine.CLIP, t1.getPrintMode());
		
		ts.setPrintModeAndPrintFeaturesForRegex(Utils.tokenize("print -full #1", " "));
		assertEquals(PrintRawLine.FULL, t1.getPrintMode());
		
		// Change number of lines but leave printmode as it is
		ts.setPrintModeAndPrintFeaturesForRegex(Utils.tokenize("print -n 10 #1", " "));
		assertEquals(PrintRawLine.FULL, t1.getPrintMode());
		
		// Set -n with track in mode OFF: -n turns the printing ON 
		ts.getTrackList().get(0).setPrintMode(PrintRawLine.OFF);
		ts.setPrintModeAndPrintFeaturesForRegex(Utils.tokenize("print -n 10", " "));
		assertTrue( ! t1.getPrintMode().equals(PrintRawLine.OFF));
	}
	
	
	@Test
	public void canSetTrackColour() throws InvalidCommandLineException, IOException, InvalidGenomicCoordsException, InvalidColourException{
				
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
		} catch(InvalidColourException e){
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
		
		// Invert selection
		cmdInput= "trackColour -v blue #1";
		ts.setTrackColourForRegex(Utils.tokenize(cmdInput, " "));
//		 assertTrue( ! ts.getTrack(t1).getTitleColour().equals("blue"));
//		 assertEquals("blue", ts.getTrack(t2).getTitleColour());
//		 assertEquals("blue", ts.getTrack(t3).getTitleColour());

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
	
		// Invert selection
		cmdInput= "trackHeight -v 10 #1";  
		ts.setTrackHeightForRegex(Utils.tokenize(cmdInput, " "));
		assertEquals(99, ts.getTrack(t1).getyMaxLines());
		assertEquals(10, ts.getTrack(t2).getyMaxLines());
		assertEquals(10, ts.getTrack(t3).getyMaxLines());
		
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

		// First reset all
		cmdInput= "ylim 0 10";
		ts.setTrackYlimitsForRegex(Utils.tokenize(cmdInput, " "));
		for(Track tr : ts.getTrackList()){
			assertEquals(0, tr.getYLimitMin(), 0.001);
			assertEquals(10, tr.getYLimitMax(), 0.001);
		}
		
		// These limits may be confused for regexes. Make sure the sublisiung is correct:
		cmdInput= "ylim -1 2 #1";
		ts.setTrackYlimitsForRegex(Utils.tokenize(cmdInput, " "));
		assertEquals(-1, ts.getTrack(t1).getYLimitMin(), 0.001);
		assertEquals(2, ts.getTrack(t1).getYLimitMax(), 0.001);
		assertEquals(0, ts.getTrack(t2).getYLimitMin(), 0.001);
		assertEquals(10, ts.getTrack(t2).getYLimitMax(), 0.001);
		assertEquals(0, ts.getTrack(t3).getYLimitMin(), 0.001);
		assertEquals(10, ts.getTrack(t3).getYLimitMax(), 0.001);

		// First reset all
		cmdInput= "ylim 0 10";
		ts.setTrackYlimitsForRegex(Utils.tokenize(cmdInput, " "));
		for(Track tr : ts.getTrackList()){
			assertEquals(0, tr.getYLimitMin(), 0.001);
			assertEquals(10, tr.getYLimitMax(), 0.001);
		}
		
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
