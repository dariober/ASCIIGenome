	package tracks;

import static org.junit.Assert.*;

import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import org.apache.commons.lang3.text.StrTokenizer;
import org.junit.Test;

import exceptions.InvalidGenomicCoordsException;
import exceptions.InvalidRecordException;
import samTextViewer.GenomicCoords;
import tracks.IntervalFeature;
import tracks.IntervalFeatureSet;

public class IntervalFeatureSetTest {

	@Test
	public void canInitialize() throws ClassNotFoundException, IOException, InvalidGenomicCoordsException, InvalidRecordException, SQLException{
		
		IntervalFeatureSet set= new IntervalFeatureSet("test_data/hg19_genes_head.gtf");
		List<IntervalFeature> x = set.getFeaturesInInterval("chr1", 1000000, 5000000);
		assertTrue(x.get(0).getRaw().startsWith("chr1"));
		assertTrue(x.get(0).getFrom() >= 1000000);
		
	}
	
	@Test
	public void canShowAndHide_getFeaturesInInterval() throws IOException, InvalidGenomicCoordsException, ClassNotFoundException, InvalidRecordException, SQLException{
		// For Map:
		IntervalFeatureSet set= new IntervalFeatureSet("test_data/hg19_genes_head.gtf");
		set.setShowRegex("start_codon");
		set.setHideRegex("OR4F");
		List<IntervalFeature> subset = set.getFeaturesInInterval("chr1", 1, 500000000);
		assertEquals(40, subset.size());

		// Same for tabix
		set= new IntervalFeatureSet("test_data/hg19_genes_head.gtf.gz");
		set.setShowRegex("start_codon");
		set.setHideRegex("OR4F");
		subset = set.getFeaturesInInterval("chr1", 1, 500000000);
		assertEquals(40, subset.size());
		
		// TESTME: coordsOfNextFeature
		// findNextString
	}

	@Test
	public void canShowAndHide_coordsOfNextFeature() throws IOException, InvalidGenomicCoordsException, ClassNotFoundException, InvalidRecordException, SQLException{
		GenomicCoords gc= new GenomicCoords("chr1:1", null, 100, null);

		IntervalFeatureSet set= new IntervalFeatureSet("test_data/hg19_genes_head.gtf");
		set.setShowRegex(".*exon.*");
		set.setHideRegex(".*DDX11L1.*");
		GenomicCoords curr = set.coordsOfNextFeature(gc);
		assertEquals(14362, (int) curr.getFrom());
	}

	@Test
	public void canShowAndHide_findNextRegex() throws IOException, InvalidGenomicCoordsException, ClassNotFoundException, InvalidRecordException, SQLException{
		GenomicCoords gc= new GenomicCoords("chr1:1", null, 100, null);

		IntervalFeatureSet set= new IntervalFeatureSet("test_data/hg19_genes_head.gtf");
		set.setShowRegex(".*exon.*");
		set.setHideRegex(".*DDX11L1.*");
		GenomicCoords curr = set.findNextMatch(gc, ".*gene_id.*");
		assertEquals(14362, (int) curr.getFrom());
	}
	
	@Test
	public void canReadFileWithComments() throws IOException, InvalidGenomicCoordsException, ClassNotFoundException, InvalidRecordException, SQLException{

		IntervalFeatureSet set= new IntervalFeatureSet("test_data/refSeq.bed");
		assertEquals(2, set.getFeaturesInInterval("chr1", 0, 100000).size());
	
	}
	
	@Test
	public void test2() throws IOException, InvalidGenomicCoordsException, ClassNotFoundException, InvalidRecordException, SQLException{
		Set<String> chroms= new HashSet<String>();
		chroms.add("c1");
		chroms.add("c2");
		chroms.add("c3");
		chroms.add("c4");
		chroms.add("c5");
		chroms.add("c6");
		
		IntervalFeatureSet set= new IntervalFeatureSet("test_data/refSeq.hg19.short.sort-2.bed");
		System.out.println(set.getChromListStartingAt(chroms, "c3"));
		System.out.println(set.getChromListStartingAt(chroms, "c1"));
		System.out.println(set.getChromListStartingAt(chroms, "c6"));
	}

	
	@Test
	public void canFindNextFeatureOnChromGivenRegex() throws IOException, InvalidGenomicCoordsException, ClassNotFoundException, InvalidRecordException, SQLException{

		StrTokenizer str= new StrTokenizer("one\n   two '   three four'");
		str.setQuoteChar('\'');
		List<String> xs = str.getTokenList();		
		
		System.out.println(xs);
		
		IntervalFeatureSet set= new IntervalFeatureSet("test_data/refSeq.hg19.short.sort-2.bed");
		
		IntervalFeature x = set.findNextRegexInGenome(".*NM_.*", "chr1", 20000000);
		assertTrue(x.getRaw().contains("NM_013943_utr3_5_0_chr1_25167429_f"));
		x = set.findNextRegexInGenome(".*NM_.*", "chr1", 80000000);
		assertTrue(x.getRaw().contains("NM_001080397_utr3_8_0_chr1_8404074_f"));
		
		x = set.findNextRegexInGenome("NotPresent", "chr1", 1);
		assertEquals(null, x);
		
		// Tabix
		set= new IntervalFeatureSet("test_data/refSeq.hg19.short.sort.bed.gz");
		x = set.findNextRegexInGenome(".*NM_.*", "chr1", 20000000);
		assertTrue(x.getRaw().contains("NM_013943_utr3_5_0_chr1_25167429_f"));
		
		x = set.findNextRegexInGenome(".*NM_.*", "chr1", 80000000);
		assertTrue(x.getRaw().contains("NM_001080397_utr3_8_0_chr1_8404074_f"));
	
		set= new IntervalFeatureSet("test_data/refSeq.hg19.short.sort-2.bed");
		x = set.findNextRegexInGenome(".*NM_.*", "chr1", 20000000);
		int i= 0;
		while(i < 20){
			x = set.findNextRegexInGenome(".*NM_.*", x.getChrom(), x.getFrom());
			assertTrue(x.getRaw().contains("NM_"));
			i++;
		}

	}

	
	@Test
	public void canGetCoordsOfNextFeature() throws IOException, InvalidGenomicCoordsException, ClassNotFoundException, InvalidRecordException, SQLException{
		IntervalFeatureSet set= new IntervalFeatureSet("test_data/refSeq.hg19.short.sort-2.bed");
		
		GenomicCoords gc= new GenomicCoords("chr1:8000000-20000000", null, 100, null);
		GenomicCoords newGc= set.coordsOfNextFeature(gc);
		assertEquals(25167428+1, (int)newGc.getFrom());
		assertEquals(25167428+gc.getGenomicWindowSize(), (int)newGc.getTo());
	
		// Handling no next feature
		gc= new GenomicCoords("chr2:8000000-20000000", null, 100, null);
		newGc= set.coordsOfNextFeature(gc);
		assertEquals(gc, newGc);
		
		gc= new GenomicCoords("chr1:100000000-101000000", null, 100, null);
		newGc= set.coordsOfNextFeature(gc);
		assertEquals(gc, newGc);
		
	}

	@Test
	public void canGetCoordsOfNextFeatureTabix() throws IOException, InvalidGenomicCoordsException, ClassNotFoundException, InvalidRecordException, SQLException{
		IntervalFeatureSet set= new IntervalFeatureSet("test_data/refSeq.hg19.short.sort.bed.gz");
		
		GenomicCoords gc= new GenomicCoords("chr1:8000000-20000000", null, 100, null);
		GenomicCoords newGc= set.coordsOfNextFeature(gc);
		assertEquals(25167428+1, (int)newGc.getFrom());
		assertEquals(25167428+gc.getGenomicWindowSize(), (int)newGc.getTo());		

		gc= new GenomicCoords("chr2:8000000-20000000", null, 100, null);
		newGc= set.coordsOfNextFeature(gc);
		assertEquals(gc, newGc);
		
		gc= new GenomicCoords("chr1:100000000-101000000", null, 100, null);
		newGc= set.coordsOfNextFeature(gc);
		assertEquals(gc, newGc);
	}
	
	@Test
	public void canFindAllRegex() throws IOException, InvalidGenomicCoordsException, ClassNotFoundException, InvalidRecordException, SQLException{
		
		IntervalFeatureSet set= new IntervalFeatureSet("test_data/hg19_genes_head.gtf");
		GenomicCoords gc= new GenomicCoords("chr18:1-10000", null, 100, null);
		GenomicCoords matched = set.genomicCoordsAllChromMatchInGenome(".*\"WASH7P\".*", gc);
		assertEquals("chr1", matched.getChrom());
		assertEquals(14362, (int)matched.getFrom());
		assertEquals(29370, (int)matched.getTo());
		
		// With tabix
		set= new IntervalFeatureSet("test_data/hg19_genes_head.gtf.gz");
		gc= new GenomicCoords("chr18:1-10000", null, 100, null);
		matched = set.genomicCoordsAllChromMatchInGenome(".*\"WASH7P\".*", gc);
		assertEquals("chr1", matched.getChrom());
		assertEquals(14362, (int)matched.getFrom());
		assertEquals(29370, (int)matched.getTo());
		
		// No match
		set= new IntervalFeatureSet("test_data/hg19_genes_head.gtf.gz");
		gc= new GenomicCoords("chr18:1-10000", null, 100, null);
		matched = set.genomicCoordsAllChromMatchInGenome(".*\"FOOBAR\".*", gc);
		assertEquals("chr18", matched.getChrom());
		assertEquals(1, (int)matched.getFrom());
		assertEquals(10000, (int)matched.getTo());
	}
	
	@Test
	public void test() throws IOException, InvalidGenomicCoordsException, ClassNotFoundException, InvalidRecordException, SQLException {
		// Note that refSeq.hg19.short.bed is not sorted by pos. 
		IntervalFeatureSet set= new IntervalFeatureSet("test_data/refSeq.hg19.short.bed");
		assertEquals("NM_001080397_utr3_8_0_chr1_8404074_f", set.getFeaturesInInterval("chr1", 0, 1000000000).get(0).getName());
		set= new IntervalFeatureSet("test_data/refSeq.hg19.bed.gz");		
	}

	@Test
	public void canFetchInterval() throws IOException, InvalidGenomicCoordsException, ClassNotFoundException, InvalidRecordException, SQLException{

		IntervalFeatureSet set= new IntervalFeatureSet("test_data/refSeq.hg19.short.bed");
		List<IntervalFeature> interval= set.getFeaturesInInterval("chr1", 20000000, 40000000);
		assertEquals(25167428+1, interval.get(0).getFrom()); // Note adding 1 because bed is 0-based
		assertEquals(33586132, interval.get(interval.size()-1).getTo());

		// Nothing to fetch: Range not in bed
		set= new IntervalFeatureSet("test_data/refSeq.hg19.bed.gz");
		interval= set.getFeaturesInInterval("chr1", 500000000, 600000000);
		assertEquals(0, interval.size());
		// Nothing to fetch: chrom not in bed:
		interval= set.getFeaturesInInterval("chrNonSense", 1, 10);
		assertEquals(0, interval.size());		
		
	}
	
	@Test
	public void canPrintMappingOfFeaturesToScreen() throws IOException, InvalidGenomicCoordsException, ClassNotFoundException, InvalidRecordException, SQLException{
		List<Double> rulerMap= new ArrayList<Double>();
		for(int i= 14000; i < 14400; i += 10){
			rulerMap.add((double)i);
		}
		System.out.println(rulerMap);
		IntervalFeatureSet set= new IntervalFeatureSet("test_data/refSeq.hg19.bed.gz");
		System.out.println(set.getFeaturesInInterval("chr1", 0, 1000000000).get(0));
	}
	
	@Test
	public void canReadVCFTabix() throws IOException, InvalidGenomicCoordsException, ClassNotFoundException, InvalidRecordException, SQLException{		
		IntervalFeatureSet set= new IntervalFeatureSet("test_data/CHD.exon.2010_03.sites.vcf.gz");
		List<IntervalFeature> xset = set.getFeaturesInInterval("1", 1, 10000000);
		assertEquals(9, xset.size());
		IntervalFeature x = xset.get(1);
		assertEquals("1", x.getChrom());
		assertEquals(1108138, x.getFrom());		
	}
	
	@Test
	public void canReadUnsortedVCF() throws IOException, InvalidGenomicCoordsException, ClassNotFoundException, InvalidRecordException, SQLException{		
		IntervalFeatureSet set= new IntervalFeatureSet("test_data/CHD.exon.2010_03.sites.unsorted.vcf");
		List<IntervalFeature> xset = set.getFeaturesInInterval("1", 1, 10000000);
		assertEquals(9, xset.size());
		IntervalFeature x = xset.get(1);
		assertEquals("1", x.getChrom());
		assertEquals(1108138, x.getFrom());
	}
	
	@Test
	public void canReadFeaturesOfLengthOne() throws IOException, InvalidGenomicCoordsException, ClassNotFoundException, InvalidRecordException, SQLException{		
		IntervalFeatureSet set= new IntervalFeatureSet("test_data/refSeqZero.bed");
		List<IntervalFeature> xset = set.getFeaturesInInterval("chr1", 1, 100);
		assertEquals(1, xset.size());
		
		// Tabix
		set= new IntervalFeatureSet("test_data/refSeqZero.bed.gz");
		xset = set.getFeaturesInInterval("chr1", 1, 100);
		assertEquals(1, xset.size());
	}
	
	@Test
	public void canReadFromURL() throws IOException, InvalidGenomicCoordsException, ClassNotFoundException, InvalidRecordException, SQLException{
		System.err.println("canReadFromURL: This can take  a while...");
		String urlStr= "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeHaibTfbs/wgEncodeHaibTfbsA549Atf3V0422111Etoh02PkRep1.broadPeak.gz";
		
		IntervalFeatureSet set= new IntervalFeatureSet(urlStr);
		List<IntervalFeature> xset = set.getFeaturesInInterval("chr1", 1, 1000000);
		assertEquals(2, xset.size());
		assertEquals(878407+1, xset.get(0).getFrom());
		
	}

}
