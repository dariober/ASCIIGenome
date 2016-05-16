package tracks;

import static org.junit.Assert.*;

import java.io.IOException;
import java.util.List;

import org.junit.Test;

import exceptions.InvalidGenomicCoordsException;
import samTextViewer.GenomicCoords;

public class TrackIntervalFeatureTest {
	
	@Test
	public void canReadTabix() throws IOException, InvalidGenomicCoordsException{
		
		String bgzFn= "test_data/refSeq.hg19.short.sort.bed.gz"; // "test_data/refSeq.hg19.short.sort.bed.gz";
		GenomicCoords gc= new GenomicCoords("chr1:16000000-20000000", null, 70, null);
		
		IntervalFeatureSet ifs= new IntervalFeatureSet(bgzFn);
		List<IntervalFeature> xset = ifs.getFeaturesInInterval(gc.getChrom(), gc.getFrom(), gc.getTo());
		assertEquals(3, xset.size());
	}
	
	@Test
	public void canPrintGtfFeatures() throws IOException, InvalidGenomicCoordsException{

		String intervalFileName= "test_data/refSeq.bed";
		//IntervalFeatureSet ifs= new IntervalFeatureSet(new File(intervalFileName)); 
		GenomicCoords gc= new GenomicCoords("chr1:1-70", null, 70, null);
		TrackIntervalFeature tif= new TrackIntervalFeature(intervalFileName, gc);
		String gtfString= tif.toString();
		System.out.println(gtfString);
	}
	
	@Test
	public void canConstructTrack() throws InvalidGenomicCoordsException, IOException {
		String intervalFileName= "test_data/refSeq.bed";
		//IntervalFeatureSet ifs= new IntervalFeatureSet(new File(intervalFileName)); 
		GenomicCoords gc= new GenomicCoords("chr1:1-70", null, 70, null);
		TrackIntervalFeature tif= new TrackIntervalFeature(intervalFileName, gc);
		tif.setNoFormat(true);           		 
		String exp= "||||||||||          ||||||||||                                        ";
		assertEquals(exp, tif.printToScreen());

		gc= new GenomicCoords("chr1:1-70", null, 35, null);
		tif= new TrackIntervalFeature(intervalFileName, gc);
		tif.setNoFormat(true);
		exp= "|||||     |||||                    ";
		assertEquals(exp, tif.printToScreen());	
	}

	@Test
	public void canAssignFeatureText() throws InvalidGenomicCoordsException, IOException {
		String intervalFileName= "test_data/hg19_genes_head.gtf";
		GenomicCoords gc= new GenomicCoords("chr1:11874-12227", null, 5, null);
		TrackIntervalFeature tif= new TrackIntervalFeature(intervalFileName, gc);
		tif.setNoFormat(true);
		assertEquals("EEEEE", tif.printToScreen());
	}
	
	@Test
	public void canStackFeatures() throws InvalidGenomicCoordsException, IOException{
		String intervalFileName= "test_data/overlapped.bed";
		int windowSize= 70;
		GenomicCoords gc= new GenomicCoords("chr1:1-70", null, windowSize, null);
		TrackIntervalFeature tif= new TrackIntervalFeature(intervalFileName, gc);
		tif.setNoFormat(true);           		 
		assertEquals(windowSize * 2 + 1, tif.printToScreen().length());
		
		String exp= "" +
"||||||||||     |||||||||||||||          ||||||||||                    \n" +
"     |||||||||||||||                                                  ";
		assertEquals(exp, tif.printToScreen());
	
	}
	
	//@Test
	public void canConstructTrackfromGtf() throws IOException, InvalidGenomicCoordsException{
		
		String intervalFileName= "test_data/hg19_genes.gtf.gz";
		//IntervalFeatureSet ifs= new IntervalFeatureSet(new File(intervalFileName)); 
		// System.out.println(ifs.getIntervalMap().size());
		GenomicCoords gc= new GenomicCoords("chr1:1-13000", null, 70, null);
		TrackIntervalFeature tif= new TrackIntervalFeature(intervalFileName, gc);
		
		gc= new GenomicCoords("chr7:5566000-5571000", null, 70, null);
		tif= new TrackIntervalFeature(intervalFileName, gc);
		System.out.println(tif);
	}
}
