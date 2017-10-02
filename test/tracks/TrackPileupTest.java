package tracks;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import org.apache.commons.io.FileUtils;
import org.junit.BeforeClass;
import org.junit.Test;

import com.google.common.base.Splitter;
import com.google.common.base.Stopwatch;

import coloring.Config;
import exceptions.InvalidColourException;
import exceptions.InvalidConfigException;
import exceptions.InvalidGenomicCoordsException;
import exceptions.InvalidRecordException;
import htsjdk.samtools.filter.MappingQualityFilter;
import htsjdk.samtools.filter.SamRecordFilter;
import samTextViewer.GenomicCoords;

public class TrackPileupTest {

    @BeforeClass
    public static void init() throws IOException, InvalidConfigException {
        new Config(null);
    }
	
	@Test
	public void canPrintConsensusSequence() throws InvalidGenomicCoordsException, IOException, ClassNotFoundException, InvalidRecordException, SQLException, InvalidColourException{

		GenomicCoords gc= new GenomicCoords("chr7:5566779-5566799", 80, null, "test_data/chr7.fa");
		TrackPileup tc= new TrackPileup("test_data/ds051.short.bam", gc);
		tc.setNoFormat(true);
		
		assertTrue(tc.getPrintableConsensusSequence().startsWith("=TT========="));
		
		// Advance coordinates and check consensus is updated:
		gc= new GenomicCoords("chr7:5566780-5566800", 80, null, "test_data/chr7.fa");
		tc.setGc(gc);
		assertTrue(tc.getPrintableConsensusSequence().startsWith("TT========="));
		
		// Large window doesn't show consensus 
		gc= new GenomicCoords("chr7:5566779-5566879", 80, null, "test_data/chr7.fa");
		tc= new TrackPileup("test_data/ds051.short.bam", gc);
		tc.setNoFormat(true);
		assertEquals("", tc.getPrintableConsensusSequence());
		
		// Region with no coverage
		gc= new GenomicCoords("chr7:1-100", 80, null, "test_data/chr7.fa");
		tc= new TrackPileup("test_data/ds051.short.bam", gc);
		tc.setNoFormat(true);
		assertEquals("", tc.getPrintableConsensusSequence());

	}

	@Test
	public void canFilterReadsContainingVariants() throws InvalidGenomicCoordsException, IOException, ClassNotFoundException, InvalidRecordException, SQLException, InvalidColourException{
		GenomicCoords gc= new GenomicCoords("chr7:1000001-1000081",100, null, "test_data/chr7.fa");
		TrackPileup tp= new TrackPileup("test_data/variant_reads.sam", gc);
		tp.setNoFormat(true);
		tp.setyMaxLines(10);
		assertEquals(10, tp.printToScreen().split("\n").length); // N. reads stacked in this interval before filtering

		tp.setVariantReadInInterval("chr7", 1000001, 1000001);
		System.err.println(tp.printToScreen());
		assertEquals(1, tp.printToScreen().trim().split("\n").length);
	}
	
	@Test
	public void canProcessReadsWithMissingSequence() throws ClassNotFoundException, IOException, InvalidGenomicCoordsException, InvalidRecordException, SQLException, InvalidColourException, InvalidConfigException{
	
		new Config(null);

		GenomicCoords gc= new GenomicCoords("chr7:1-1000", 80, null, null);
		TrackPileup tr= new TrackPileup("test_data/missingReadSeq.bam", gc);
		tr.setNoFormat(true);
		tr.printToScreen();
		assertTrue(tr.printToScreen().trim().startsWith("_"));

		assertEquals(1, (int)tr.getDepth(gc.getChrom(), gc.getFrom(), gc.getTo()).entrySet().iterator().next().getValue());

	}
	
	@Test
	public void canPrintProfile() throws ClassNotFoundException, IOException, InvalidGenomicCoordsException, InvalidRecordException, SQLException, InvalidColourException, InvalidConfigException{
	
		new Config(null);
	
		GenomicCoords gc= new GenomicCoords("chr7:5566776-5566796", 80, null, null);
		TrackPileup tr= new TrackPileup("test_data/ds051.short.bam", gc);
		System.err.println(tr.getScreenScores());
		tr.setNoFormat(true);
		System.err.println(tr.printToScreen());
		assertTrue(tr.getScreenScores().size() > 1); // Here we only test the method doesn't crash
		
		gc= new GenomicCoords("chr7:5,554,740-5,554,780", 80, null, null);
		tr= new TrackPileup("test_data/ear045.oxBS.actb.bam", gc);
		assertTrue(tr.getScreenScores().size() > 1);
		
		assertTrue(tr.printToScreen().length() > 50);
		assertTrue(tr.getTitle().length() > 50);
		System.err.println(tr.printToScreen());
		System.err.println(tr.getTitle());
		
		tr.setRpm(true);
		tr.getScreenScores();
		
	}

	@Test
	public void canCorrectRangeByRpm() throws ClassNotFoundException, IOException, InvalidGenomicCoordsException, InvalidRecordException, SQLException, InvalidColourException, InvalidConfigException{
	
		GenomicCoords gc= new GenomicCoords("chr7:5566776-5566796", 80, null, null);
		TrackPileup tr= new TrackPileup("test_data/ds051.short.bam", gc);
		System.err.println(tr.getScreenScores());
		assertTrue(tr.getTitle().contains("range[1.0 22.0]"));
		tr.setRpm(true);
		System.err.println(tr.getScreenScores());
		assertTrue( tr.getTitle().contains("1000000")); // The range should contain 1,000,000 because this is the entire size of the file
		                                                // rpm= 22/22*1,000,000
	}

	
	@Test
	public void canCollectCoverage() throws ClassNotFoundException, IOException, InvalidGenomicCoordsException, InvalidRecordException, SQLException{
		GenomicCoords gc= new GenomicCoords("chr7:5566736-5566856", 80, null, null);
		TrackPileup tr= new TrackPileup("test_data/ds051.short.bam", gc);
		
		assertEquals(79, tr.getDepth(gc.getChrom(), gc.getFrom(), gc.getTo()).size());
		assertEquals(1, (int)tr.getDepth(gc.getChrom(), gc.getFrom(), gc.getTo()).get(5566778)); // Depths checked against mpileup
		assertEquals(5, (int)tr.getDepth(gc.getChrom(), gc.getFrom(), gc.getTo()).get(5566782));
		assertEquals(18, (int)tr.getDepth(gc.getChrom(), gc.getFrom(), gc.getTo()).get(5566856));
		
		gc= new GenomicCoords("chr7:5522059-5612125", 80, null, null);
		long t0= System.currentTimeMillis();
		tr= new TrackPileup("test_data/ear045.oxBS.actb.bam", gc);
		Map<Integer, Integer> depth = tr.getDepth(gc.getChrom(), gc.getFrom(), gc.getTo());
		long t1= System.currentTimeMillis();
		assertTrue(t1-t0 < 10000); // Processing time (in ms) is acceptably small
		assertTrue(t1-t0 > 100); // But not suspiciously small
		
		System.err.println("Time to parse " + depth.size() + " positions: " + (t1-t0) + " ms");
	}
	
	public static void sameAsMpileup() throws ClassNotFoundException, IOException, InvalidGenomicCoordsException, InvalidRecordException, SQLException{

		GenomicCoords gc= new GenomicCoords("chr7:5522059-5612125", 80, null, null);
		TrackPileup tr= new TrackPileup("test_data/ear045.oxBS.actb.bam", gc);
		Map<Integer, Integer> depth = tr.getDepth(gc.getChrom(), gc.getFrom(), gc.getTo());

		// See test_data/README.md for obtaining this test file (samtools mpileup ...)
		String expPileup= FileUtils.readFileToString(new File("test_data/ear045.oxBS.actb.pileup"));
		List<String> expList = Splitter.on("\n").omitEmptyStrings().splitToList(expPileup);
	
		// mpileup and TrackPileup hit the same positions
		assertEquals(expList.size(), depth.size());

		// Same depth at same positions
		int i= 0;
		for(int obsPos : depth.keySet()){
			int obsDepth=depth.get(obsPos);
			int expPos= Integer.parseInt(Splitter.on("\t").splitToList(expList.get(i)).get(1));
			int expDepth= Integer.parseInt(Splitter.on("\t").splitToList(expList.get(i)).get(3));
			try{
				assertEquals(expPos, obsPos);
				assertTrue(Math.abs((expDepth - obsDepth)) < 6);
			} catch(AssertionError e){
				System.err.println("At iteration: " + i);
				System.err.println(expList.get(i));
				System.err.println("Observed depth: " + obsDepth);
				throw e; 
			}
			i++;
		}				
	}
	
	@Test
	public void canFilterReadsWithGrep() throws InvalidGenomicCoordsException, IOException, ClassNotFoundException, InvalidRecordException, SQLException, InvalidColourException{
		GenomicCoords gc= new GenomicCoords("chr7:5566000-5567000",80, null, null);
		TrackPileup tr= new TrackPileup("test_data/ds051.short.bam", gc);
		tr.setNoFormat(true);
		tr.setYLimitMin(0);
		tr.setYLimitMax(10);
		tr.setShowHideRegex("NCNNNCCC", "\\t5566779\\t");
		assertEquals(4, tr.printToScreen().trim().split("\n").length);
	}

	@Test
	public void canFilterReadsSamtools() throws InvalidGenomicCoordsException, IOException, ClassNotFoundException, InvalidRecordException, SQLException, InvalidColourException{
		
		GenomicCoords gc= new GenomicCoords("chr7:5566778",80, null, null);
		TrackPileup tr= new TrackPileup("test_data/ds051.short.bam", gc);
		tr.setNoFormat(true);
		assertTrue(tr.getTitle().contains("22.0")); // Before filtering
		List<SamRecordFilter> samRecordFilter= new ArrayList<SamRecordFilter>();
		samRecordFilter.add(new MappingQualityFilter(51));
		tr.setSamRecordFilter(samRecordFilter);
		
		assertEquals("", tr.printToScreen().trim());

		samRecordFilter= new ArrayList<SamRecordFilter>();
		samRecordFilter.add(new MappingQualityFilter(11));
		tr.setSamRecordFilter(samRecordFilter);
		assertTrue(tr.getTitle().contains("17.0"));
		
	}

	@Test
	public void canFilterReadsWithGrepAndAwk() throws InvalidGenomicCoordsException, IOException, ClassNotFoundException, InvalidRecordException, SQLException, InvalidColourException{
		GenomicCoords gc= new GenomicCoords("chr7:5566000-5567000",80, null, null);
		TrackPileup tr= new TrackPileup("test_data/ds051.short.bam", gc);
		tr.setNoFormat(true);
		tr.setyMaxLines(1000);
		assertTrue(tr.getTitle().contains("22.0")); // N. reads before filtering
		// assertTrue(tr.getTitle().contains("22/22")); // N. reads before filtering
		tr.setShowHideRegex("NCNNNCCC", FeatureFilter.DEFAULT_HIDE_REGEX);
		tr.setAwk("'$4 != 5566779'");
		System.err.println(tr.getTitle());
		assertTrue(tr.getTitle().contains("4.0"));
		// assertTrue(tr.getTitle().contains("4/22"));
	}

	@Test
	public void canMergePositionsIntoIntervals(){
		
		List<Integer> x= new ArrayList<Integer>();
		assertEquals(0, TrackPileup.mergePositionsInIntervals(x).size());
		
		x.clear();
		x.add(1); x.add(3); x.add(9);
		assertEquals("[[1, 1], [3, 3], [9, 9]]", TrackPileup.mergePositionsInIntervals(x).toString());

		x.clear(); // Handle duplicates
		x.add(1); x.add(1); x.add(2); x.add(3); x.add(3);
		assertEquals("[[1, 3]]", TrackPileup.mergePositionsInIntervals(x).toString());
		
		x.clear();
		x.add(1); x.add(2); x.add(3); x.add(4);
		assertEquals("[[1, 4]]", TrackPileup.mergePositionsInIntervals(x).toString());
		
		x.clear();
		x.add(1); x.add(2); x.add(4); x.add(5); x.add(6); x.add(10); x.add(11);
		assertEquals("[[1, 2], [4, 6], [10, 11]]", TrackPileup.mergePositionsInIntervals(x).toString());

		x.clear();
		x.add(1); x.add(2); x.add(4); x.add(5); x.add(6); x.add(11);
		assertEquals("[[1, 2], [4, 6], [11, 11]]", TrackPileup.mergePositionsInIntervals(x).toString());
		
		x.clear();
		x.add(2); x.add(3); x.add(1);
		boolean pass= false;
		try{
			TrackPileup.mergePositionsInIntervals(x);
		} catch(RuntimeException e){
			pass= true;
		}
		assertTrue(pass);
		
		Stopwatch sw= Stopwatch.createStarted();
		List<Integer> bigList= new ArrayList<Integer>();
		for(int i= 0; i < 10000000; i++){
			bigList.add(i);
		}
		System.err.println(sw);
		sw.reset(); sw.start();
		TrackPileup.mergePositionsInIntervals(bigList);
	}
	
	@Test
	public void canCollectCoverageAtOnePos() throws ClassNotFoundException, IOException, InvalidGenomicCoordsException, InvalidRecordException, SQLException{
		GenomicCoords gc= new GenomicCoords("chr7:5588536-5588536", 80, null, null);
		TrackPileup tr= new TrackPileup("test_data/ear045.oxBS.actb.bam", gc);
		assertEquals(1, tr.getDepth(gc.getChrom(), gc.getFrom(), gc.getTo()).size());
	}
	
	@Test
	public void canHandleZeroReads() throws ClassNotFoundException, IOException, InvalidGenomicCoordsException, InvalidRecordException, SQLException{
		GenomicCoords gc= new GenomicCoords("chr1:1-1000", 80, null, null);
		TrackPileup tr= new TrackPileup("test_data/ear045.oxBS.actb.bam", gc);
		assertEquals(0, tr.getDepth(gc.getChrom(), gc.getFrom(), gc.getTo()).size());
	}

	@Test
	public void canConstructFromUnsortedInput() throws ClassNotFoundException, IOException, InvalidGenomicCoordsException, InvalidRecordException, SQLException{
		GenomicCoords gc= new GenomicCoords("chr1:1-1000", 80, null, null);
		new TrackPileup("/Users/db291g/Tritume/reads.sam", gc);
	}	
}
