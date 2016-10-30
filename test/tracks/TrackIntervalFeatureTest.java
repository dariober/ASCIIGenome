package tracks;

import static org.junit.Assert.*;

import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.List;
import org.junit.Test;

import exceptions.InvalidGenomicCoordsException;
import exceptions.InvalidRecordException;
import samTextViewer.GenomicCoords;

public class TrackIntervalFeatureTest {
	
	
	@Test
	public void canReadTabix() throws IOException, InvalidGenomicCoordsException, ClassNotFoundException, InvalidRecordException, SQLException{
		
		String bgzFn= "test_data/refSeq.hg19.short.sort.bed.gz"; // "test_data/refSeq.hg19.short.sort.bed.gz";
		GenomicCoords gc= new GenomicCoords("chr1:16000000-20000000", null, null);
		TrackIntervalFeature tif= new TrackIntervalFeature(bgzFn, gc);
		List<IntervalFeature> xset = tif.getFeaturesInInterval(gc.getChrom(), gc.getFrom(), gc.getTo());
		assertEquals(3, xset.size());
		
	}
	
	@Test
	public void canPrintGtfFeatures() throws IOException, InvalidGenomicCoordsException, ClassNotFoundException, InvalidRecordException, SQLException{

		String intervalFileName= "test_data/refSeq.bed";
		GenomicCoords gc= new GenomicCoords("chr1:1-70", null, null);
		TrackIntervalFeature tif= new TrackIntervalFeature(intervalFileName, gc);
		String gtfString= tif.toString();
		System.out.println(gtfString);
	}
		
	@Test
	public void canConstructTrack() throws InvalidGenomicCoordsException, IOException, ClassNotFoundException, InvalidRecordException, SQLException {
		String intervalFileName= "test_data/refSeq.bed";
		//IntervalFeatureSet ifs= new IntervalFeatureSet(new File(intervalFileName)); 
		GenomicCoords gc= new GenomicCoords("chr1:1-70", null, null);
		TrackIntervalFeature tif= new TrackIntervalFeature(intervalFileName, gc);
		tif.setNoFormat(true);           		 		
		String exp= "||||||||||          ||||||||||                                        ";
		assertEquals(exp, tif.printToScreen());

		gc= new GenomicCoords("chr1:1-70", null, null);
		tif= new TrackIntervalFeature(intervalFileName, gc);
		tif.setNoFormat(true);
		assertTrue(tif.printToScreen().startsWith("||||"));	
		assertTrue(tif.printToScreen().endsWith("    "));	

	}

	@Test
	public void canAssignFeatureText() throws InvalidGenomicCoordsException, IOException, ClassNotFoundException, InvalidRecordException, SQLException {
		String intervalFileName= "test_data/hg19_genes_head.gtf";
		GenomicCoords gc= new GenomicCoords("chr1:11874-12227", null, null);
		TrackIntervalFeature tif= new TrackIntervalFeature(intervalFileName, gc);
		tif.setNoFormat(true);
		tif.setGtfAttributeForName("N/A");
		assertTrue(tif.printToScreen().startsWith("EEEEE"));
	}

	@Test
	public void canPrintFeatureWithName() throws InvalidGenomicCoordsException, IOException, ClassNotFoundException, InvalidRecordException, SQLException{
		
		String intervalFileName= "test_data/hg19_genes_head.gtf";
		GenomicCoords gc= new GenomicCoords("chr1:11874-12052", null, null);
		TrackIntervalFeature tif= new TrackIntervalFeature(intervalFileName, gc);
		tif.setNoFormat(true);
		tif.setGtfAttributeForName(null);
		assertTrue(tif.printToScreen().startsWith("NR_046018_1_EEEEEEEE"));
		
	}
	
	@Test
	public void canStackFeatures() throws InvalidGenomicCoordsException, IOException, ClassNotFoundException, InvalidRecordException, SQLException{
		String intervalFileName= "test_data/overlapped.bed";
		GenomicCoords gc= new GenomicCoords("chr1:1-70", null, null);
		TrackIntervalFeature tif= new TrackIntervalFeature(intervalFileName, gc);
		tif.setNoFormat(true);           		 
		assertEquals(70 * 2 + 1, tif.printToScreen().length());
		
		String exp= "" +
"||||||||||     |||||||||||||||          ||||||||||                    \n" +
"     |||||||||||||||                                                  ";
		assertEquals(exp, tif.printToScreen());
	
	}
	
	@Test
	public void canConstructTrackfromGtf() throws IOException, InvalidGenomicCoordsException, ClassNotFoundException, InvalidRecordException, SQLException{
		
		String intervalFileName= "test_data/hg19_genes.gtf.gz";
		GenomicCoords gc= new GenomicCoords("chr1:1-13000", null, null);
		TrackIntervalFeature tif= new TrackIntervalFeature(intervalFileName, gc);
		
		gc= new GenomicCoords("chr7:5566000-5571000", null, null);
		tif= new TrackIntervalFeature(intervalFileName, gc);
		System.out.println(tif);
	}
	
	@Test
	public void canShowHideFeatureByRegexFilters() throws InvalidGenomicCoordsException, IOException, ClassNotFoundException, InvalidRecordException, SQLException{
		
		String intervalFileName= "test_data/hg19_genes_head.gtf.gz";
		GenomicCoords gc= new GenomicCoords("chr1:10000-100000", null, null);
		TrackIntervalFeature tif= new TrackIntervalFeature(intervalFileName, gc);
		tif.setNoFormat(true);
		tif.setHideRegex("\texon\t");
		assertTrue(tif.getIntervalFeatureList().size() == 3);

		tif.setHideRegex("^$");
		tif.setShowRegex("WASH7P");
		assertTrue(tif.getIntervalFeatureList().size() == 11);
	}
	
	@Test
	public void canShowAndHide_getFeaturesInInterval() throws IOException, InvalidGenomicCoordsException, ClassNotFoundException, InvalidRecordException, SQLException{
		
		String intervalFileName= "test_data/hg19_genes_head.gtf.gz";
		GenomicCoords gc= new GenomicCoords("chr1:10000-100000", null, null);
		TrackIntervalFeature tif= new TrackIntervalFeature(intervalFileName, gc);

		tif.setShowRegex("start_codon");
		tif.setHideRegex("OR4F");
		List<IntervalFeature> subset = tif.getFeaturesInInterval("chr1", 1, 500000000);
		assertEquals(40, subset.size());

	}


	@Test
	public void canShowAndHide_coordsOfNextFeature() throws IOException, InvalidGenomicCoordsException, ClassNotFoundException, InvalidRecordException, SQLException{
		
		GenomicCoords gc= new GenomicCoords("chr1:1", null, null);
		String intervalFileName= "test_data/hg19_genes_head.gtf.gz";
		TrackIntervalFeature tif= new TrackIntervalFeature(intervalFileName, gc);

		
		tif.setShowRegex(".*exon.*");
		tif.setHideRegex(".*DDX11L1.*");
		GenomicCoords curr = tif.coordsOfNextFeature(gc);
		assertEquals(14362, (int) curr.getFrom());
		
	}

	@Test
	public void canShowAndHide_findNextRegex() throws IOException, InvalidGenomicCoordsException, ClassNotFoundException, InvalidRecordException, SQLException{
		
		GenomicCoords gc= new GenomicCoords("chr1:1", null, null);
		String intervalFileName= "test_data/hg19_genes_head.gtf.gz";
		TrackIntervalFeature tif= new TrackIntervalFeature(intervalFileName, gc);
		
		tif.setShowRegex(".*exon.*");
		tif.setHideRegex(".*DDX11L1.*");
		GenomicCoords curr = tif.findNextMatch(gc, ".*gene_id.*");
		assertEquals(14362, (int) curr.getFrom());
	}

	@Test
	public void canReadFileWithComments() throws IOException, InvalidGenomicCoordsException, ClassNotFoundException, InvalidRecordException, SQLException{

		GenomicCoords gc= new GenomicCoords("chr1:1-1000000", null, null);
		String intervalFileName= "test_data/refSeq.bed";
		TrackIntervalFeature tif= new TrackIntervalFeature(intervalFileName, gc);
		assertEquals(2, tif.getFeaturesInInterval("chr1", 0, 100000).size());
	
	}

	@Test
	public void canFindNextFeatureOnChromGivenRegex() throws IOException, InvalidGenomicCoordsException, ClassNotFoundException, InvalidRecordException, SQLException{

		GenomicCoords gc= new GenomicCoords("chr1:1-1000000", null, null);
		String intervalFileName= "test_data/refSeq.hg19.short.sort-2.bed";
		TrackIntervalFeature tif= new TrackIntervalFeature(intervalFileName, gc);
		
		IntervalFeature x = tif.findNextRegexInGenome(".*NM_.*", "chr1", 20000000);
		assertTrue(x.getRaw().contains("NM_013943_utr3_5_0_chr1_25167429_f"));
		x = tif.findNextRegexInGenome(".*NM_.*", "chr1", 80000000);
		assertTrue(x.getRaw().contains("NM_001080397_utr3_8_0_chr1_8404074_f"));
		
		x = tif.findNextRegexInGenome("NotPresent", "chr1", 1);
		assertEquals(null, x);
		
	}

	@Test
	public void canGetCoordsOfNextFeature() throws IOException, InvalidGenomicCoordsException, ClassNotFoundException, InvalidRecordException, SQLException{

		GenomicCoords gc= new GenomicCoords("chr1:8000000-20000000", null, null);
		TrackIntervalFeature tif= new TrackIntervalFeature("test_data/refSeq.hg19.short.sort-2.bed", gc);
		
		GenomicCoords newGc= tif.coordsOfNextFeature(gc);
		assertEquals(25167428+1, (int)newGc.getFrom());
		assertEquals(25167428+gc.getGenomicWindowSize(), (int)newGc.getTo());
	
		// Handling no next feature
		gc= new GenomicCoords("chr2:8000000-20000000", null, null);
		newGc= tif.coordsOfNextFeature(gc);
		assertEquals(gc, newGc);
		
		gc= new GenomicCoords("chr1:100000000-101000000", null, null);
		newGc= tif.coordsOfNextFeature(gc);
		assertEquals(gc, newGc);
		
	}

	@Test
	public void canFindAllRegex() throws IOException, InvalidGenomicCoordsException, ClassNotFoundException, InvalidRecordException, SQLException{

		GenomicCoords gc= new GenomicCoords("chr18:1-10000", null, null);
		TrackIntervalFeature tif= new TrackIntervalFeature("test_data/hg19_genes_head.gtf", gc);
		
		GenomicCoords matched = tif.genomicCoordsAllChromMatchInGenome(".*\"WASH7P\".*", gc);
		assertEquals("chr1", matched.getChrom());
		assertEquals(14362, (int)matched.getFrom());
		assertEquals(29370, (int)matched.getTo());
				
		// No match
		gc= new GenomicCoords("chr18:1-10000", null, null);
		tif= new TrackIntervalFeature("test_data/hg19_genes_head.gtf", gc);
		matched = tif.genomicCoordsAllChromMatchInGenome(".*\"FOOBAR\".*", gc);
		assertEquals("chr18", matched.getChrom());
		assertEquals(1, (int)matched.getFrom());
		assertEquals(10000, (int)matched.getTo());
	}

	@Test
	public void canFetchInterval() throws IOException, InvalidGenomicCoordsException, ClassNotFoundException, InvalidRecordException, SQLException{

		GenomicCoords gc= new GenomicCoords("chr18:1-10000", null, null);
		TrackIntervalFeature tif= new TrackIntervalFeature("test_data/refSeq.hg19.short.bed", gc);

		List<IntervalFeature> interval= tif.getFeaturesInInterval("chr1", 20000000, 40000000);
		assertEquals(25167428+1, interval.get(0).getFrom()); // Note adding 1 because bed is 0-based
		assertEquals(33586132, interval.get(interval.size()-1).getTo());

		// Nothing to fetch: Range not in bed
		tif= new TrackIntervalFeature("test_data/refSeq.hg19.short.bed", gc);
		interval= tif.getFeaturesInInterval("chr1", 500000000, 600000000);
		assertEquals(0, interval.size());
		// Nothing to fetch: chrom not in bed:
		interval= tif.getFeaturesInInterval("chrNonSense", 1, 10);
		assertEquals(0, interval.size());		
		
	}

	@Test
	public void canPrintMappingOfFeaturesToScreen() throws IOException, InvalidGenomicCoordsException, ClassNotFoundException, InvalidRecordException, SQLException{
		List<Double> rulerMap= new ArrayList<Double>();
		for(int i= 14000; i < 14400; i += 10){
			rulerMap.add((double)i);
		}
		System.out.println(rulerMap);

		GenomicCoords gc= new GenomicCoords("chr18:1-10000", null, null);
		TrackIntervalFeature tif= new TrackIntervalFeature("test_data/refSeq.hg19.short.bed", gc);
		
		System.out.println(tif.getFeaturesInInterval("chr1", 0, 1000000000).get(0));
	}

	@Test
	public void canReadVCFTabix() throws IOException, InvalidGenomicCoordsException, ClassNotFoundException, InvalidRecordException, SQLException{		

		GenomicCoords gc= new GenomicCoords("chr18:1-10000", null, null);
		TrackIntervalFeature tif= new TrackIntervalFeature("test_data/CHD.exon.2010_03.sites.vcf.gz", gc);

		List<IntervalFeature> xset = tif.getFeaturesInInterval("1", 1, 10000000);
		assertEquals(9, xset.size());
		IntervalFeature x = xset.get(1);
		assertEquals("1", x.getChrom());
		assertEquals(1108138, x.getFrom());		
	}

	
	@Test
	public void canReadUnsortedVCF() throws IOException, InvalidGenomicCoordsException, ClassNotFoundException, InvalidRecordException, SQLException{		

		GenomicCoords gc= new GenomicCoords("chr18:1-10000", null, null);
		TrackIntervalFeature tif= new TrackIntervalFeature("test_data/CHD.exon.2010_03.sites.unsorted.vcf", gc);
		List<IntervalFeature> xset = tif.getFeaturesInInterval("1", 1, 10000000);
		assertEquals(9, xset.size());
		IntervalFeature x = xset.get(1);
		assertEquals("1", x.getChrom());
		assertEquals(1108138, x.getFrom());
	}


	@Test
	public void canReadFeaturesOfLengthOne() throws IOException, InvalidGenomicCoordsException, ClassNotFoundException, InvalidRecordException, SQLException{		

		GenomicCoords gc= new GenomicCoords("chr18:1-10000", null, null);
		TrackIntervalFeature tif= new TrackIntervalFeature("test_data/refSeqZero.bed", gc);		
		
		List<IntervalFeature> xset = tif.getFeaturesInInterval("chr1", 1, 100);
		assertEquals(1, xset.size());
		
	}

	@Test
	public void canReadFromURL() throws IOException, InvalidGenomicCoordsException, ClassNotFoundException, InvalidRecordException, SQLException{

		System.err.println("canReadFromURL: This can take  a while...");
		String urlStr= "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeHaibTfbs/wgEncodeHaibTfbsA549Atf3V0422111Etoh02PkRep1.broadPeak.gz";

		GenomicCoords gc= new GenomicCoords("chr18:1-10000", null, null);
		TrackIntervalFeature tif= new TrackIntervalFeature(urlStr, gc);		
		
		List<IntervalFeature> xset = tif.getFeaturesInInterval("chr1", 1, 1000000);
		assertEquals(2, xset.size());
		assertEquals(878407+1, xset.get(0).getFrom());
		
	}


	
}
