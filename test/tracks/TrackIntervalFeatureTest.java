package tracks;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.List;

import org.junit.Test;

import exceptions.InvalidColourException;
import exceptions.InvalidGenomicCoordsException;
import exceptions.InvalidRecordException;
import samTextViewer.GenomicCoords;

public class TrackIntervalFeatureTest {
	
	@Test
	public void canPrintChromsomeNames() throws ClassNotFoundException, IOException, InvalidGenomicCoordsException, InvalidRecordException, SQLException{
		
		String intervalFileName= "test_data/hg19_genes.gtf.gz";
		GenomicCoords gc= new GenomicCoords("7:5527151-5530709", null, null);
		TrackIntervalFeature tif= new TrackIntervalFeature(intervalFileName, gc);
		
		assertTrue(tif.getChromosomeNames().size() > 10);
		
		tif= new TrackIntervalFeature("test_data/wgEncodeDukeDnase8988T.fdr01peaks.hg19.bb", gc);
		assertTrue(tif.getChromosomeNames().size() > 10);
	}

	@Test
	public void transcriptGTFToOneLine() throws InvalidGenomicCoordsException, IOException, ClassNotFoundException, InvalidRecordException, SQLException{
		
		// Expect:
		// ccccccccccc-----cccccae------------------------------------eeeeeeeee
		
		String intervalFileName= "test_data/hg19.gencode_genes_v19.gtf.gz";
		GenomicCoords gc= new GenomicCoords("chr7:5568562-5572120", null, null);
		TrackIntervalFeature tif= new TrackIntervalFeature(intervalFileName, gc);
		tif.setNoFormat(true);
		assertTrue(tif.printToScreen().trim().startsWith("ccc"));
		assertTrue(tif.printToScreen().trim().endsWith("eee"));
		assertEquals(6, tif.getIntervalFeatureList().size());

	}
	
	@Test
	public void canHideTitle() throws ClassNotFoundException, InvalidColourException, IOException, InvalidGenomicCoordsException, InvalidRecordException, SQLException{
		// See issue #42
		String intervalFileName= "test_data/hg19.gencode_genes_v19.gtf.gz";
		GenomicCoords gc= new GenomicCoords("chr7", null, null);
		TrackIntervalFeature tif= new TrackIntervalFeature(intervalFileName, gc);
		tif.setHideTitle(true);
		assertEquals("", tif.getTitle());
	}
	
	@Test
	public void transcriptGFFToOneLine() throws InvalidGenomicCoordsException, IOException, ClassNotFoundException, InvalidRecordException, SQLException{
		
		String intervalFileName= "test_data/Homo_sapiens.GRCh38.86.ENST00000331789.gff3";
		GenomicCoords gc= new GenomicCoords("7:5527151-5530709", null, null);
		TrackIntervalFeature tif= new TrackIntervalFeature(intervalFileName, gc);
		tif.setNoFormat(true);
		assertTrue(tif.printToScreen().startsWith("uuuuu"));
		assertTrue(tif.printToScreen().endsWith("www"));
		System.out.println("PRINTING:" + tif.printToScreen());
		
		tif.setNoFormat(false);
		assertTrue(tif.printToScreen().trim().startsWith("["));
	}
	
	@Test
	public void canAddNameToGFFTranscript() throws InvalidGenomicCoordsException, IOException, ClassNotFoundException, InvalidRecordException, SQLException{
		String intervalFileName= "test_data/Homo_sapiens.GRCh38.86.ENST00000331789.gff3";

		GenomicCoords gc= new GenomicCoords("7:5527151-5530709", null, null);
		TrackIntervalFeature tif= new TrackIntervalFeature(intervalFileName, gc);
		tif.setNoFormat(true);
		
		assertTrue(tif.printToScreen().contains("ACTB-001"));
		
		tif.setNoFormat(false);
		assertTrue(tif.printToScreen().trim().startsWith("["));

		
	}
	
	@Test
	public void canReadBigBed() throws InvalidGenomicCoordsException, IOException, ClassNotFoundException, InvalidRecordException, SQLException{
		
		String filename= "test_data/wgEncodeDukeDnase8988T.fdr01peaks.hg19.bb";

		GenomicCoords gc= new GenomicCoords("chr1:1-800170", null, null);
		TrackIntervalFeature tif= new TrackIntervalFeature(filename, gc);
		tif.setNoFormat(true);
		assertEquals(12, tif.getIntervalFeatureList().size());
		assertEquals(564665+1, tif.getIntervalFeatureList().get(0).getFrom());
		
	}
	
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
		tif.setNoFormat(true);
		assertEquals(2, tif.getIntervalFeatureList().size());
		
		assertEquals("||||", tif.printToScreen().substring(0,  4));

		tif.setNoFormat(false);
		assertTrue(tif.printToScreen().trim().startsWith("[")); // trim is necessary to remove escape \033
		assertTrue(tif.printToScreen().length() > 100);
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
		tif.setGtfAttributeForName("-na");
		assertTrue(tif.printToScreen().startsWith("EEEEE"));
	}

	@Test
	public void canStackFeatures() throws InvalidGenomicCoordsException, IOException, ClassNotFoundException, InvalidRecordException, SQLException{
		String intervalFileName= "test_data/overlapped.bed";
		GenomicCoords gc= new GenomicCoords("chr1:1-70", null, null);
		TrackIntervalFeature tif= new TrackIntervalFeature(intervalFileName, gc);
		tif.setNoFormat(true);           		 
		assertEquals(70 * 2 + 1, tif.printToScreen().length());
		
		String exp= "" +
">>>>>>>>>>     >>>>>>>>>>>>>>>          <<<<<<<<<<                    \n" +
"     >>>>>>>>>>>>>>>                                                  ";
		assertEquals(exp, tif.printToScreen());
	
	}

	@Test
	public void canStackFeaturesInOneLine() throws InvalidGenomicCoordsException, IOException, ClassNotFoundException, InvalidRecordException, SQLException{
		String intervalFileName= "test_data/overlapped.bed";
		GenomicCoords gc= new GenomicCoords("chr1:1-70", null, null);
		TrackIntervalFeature tif= new TrackIntervalFeature(intervalFileName, gc);
		tif.setNoFormat(true);  
		tif.setFeatureDisplayMode(FeatureDisplayMode.ONELINE);
		
		String exp= ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>          <<<<<<<<<<";

		assertEquals(exp, tif.printToScreen().trim());
	
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
	public void canApplyAwk_getFeaturesInInterval() throws IOException, InvalidGenomicCoordsException, ClassNotFoundException, InvalidRecordException, SQLException{
		
		String intervalFileName= "test_data/hg19_genes_head.gtf.gz";
		GenomicCoords gc= new GenomicCoords("chr1:10000-100000", null, null);
		TrackIntervalFeature tif= new TrackIntervalFeature(intervalFileName, gc);
		
		String awk= "'$3 == \"start_codon\" && $9 !~ \"OR4F\"'";
		tif.setAwk(awk); // Note use single quotes
	
		assertEquals("-F '\\t' " + awk, tif.getAwk()); // Check -F arg has been prepended.
		
		List<IntervalFeature> subset = tif.getFeaturesInInterval("chr1", 1, 500000000);
		assertEquals(40, subset.size());

		tif.setAwk("-F \\t '($5 - $4) > 1000'"); // Filter for feature size > x
		subset = tif.getFeaturesInInterval("chr1", 1, 500000000);
		assertEquals(23, subset.size());	
		
		tif.setAwk("  "); // Remove filter w/o args.
		subset = tif.getFeaturesInInterval("chr1", 1, 500000000);
		assertEquals(1000, subset.size());	
		
		// Invalid script: Ugly stackTrace printed. All records returned
		tif.setAwk("$foo");
		subset = tif.getFeaturesInInterval("chr1", 1, 500000000);
		assertEquals(1000, subset.size());

		// awk output is neither empty nor equal to input
		// Exception expected.
		boolean pass= false;
		try{
			tif.setAwk("'{print 999}'");
		} catch(InvalidGenomicCoordsException e){
			pass= true;
		}
		assertTrue(pass);
		
	}

	@Test
	public void canApplyAwkAndGrep_getFeaturesInInterval() throws IOException, InvalidGenomicCoordsException, ClassNotFoundException, InvalidRecordException, SQLException{
		
		String intervalFileName= "test_data/hg19_genes_head.gtf.gz";
		GenomicCoords gc= new GenomicCoords("chr1:10000-100000", null, null);
		TrackIntervalFeature tif= new TrackIntervalFeature(intervalFileName, gc);

		tif.setAwk("'$3 == \"start_codon\"");
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
		GenomicCoords curr = tif.coordsOfNextFeature(gc, false);
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
		
		GenomicCoords newGc= tif.coordsOfNextFeature(gc, false);
		assertEquals(25167428+1, (int)newGc.getFrom());
		assertEquals(25167428+gc.getGenomicWindowSize(), (int)newGc.getTo());
	
		// Next feature is on next chrom, current chrom not in file at all.
		gc= new GenomicCoords("foo:1-10000", null, null);
		newGc= tif.coordsOfNextFeature(gc, false);
		assertEquals("chr1", newGc.getChrom());
		
		gc= new GenomicCoords("chr1:100000000-101000000", null, null);
		newGc= tif.coordsOfNextFeature(gc, false);
		assertEquals("chr3", newGc.getChrom());
		
		gc= new GenomicCoords("chr1:10000000-10001000", null, null);
		newGc= tif.coordsOfNextFeature(gc, true);
		assertEquals(8404074, (int)newGc.getFrom());

		newGc= tif.coordsOfNextFeature(newGc, true);
		assertEquals(67208779, (int)newGc.getFrom());
	}

	@Test
	public void canGetCoordsOfPreviousFeature() throws IOException, InvalidGenomicCoordsException, ClassNotFoundException, InvalidRecordException, SQLException{

		GenomicCoords gc= new GenomicCoords("chr1:8000000-20000000", null, null);
		TrackIntervalFeature tif= new TrackIntervalFeature("test_data/refSeq.hg19.short.sort-2.bed", gc);
		
		gc= new GenomicCoords("chr1:10000000-10001000", null, null);
		GenomicCoords newGc= tif.coordsOfNextFeature(gc, true);
		assertEquals(8404074, (int)newGc.getFrom());

		newGc= tif.coordsOfNextFeature(newGc, true);
		assertEquals(67208779, (int)newGc.getFrom());
		
		// Exactly at the start of a feature and move to previous one. This is the feature:
		// chrM hg19_wgEncodeGencodeBasicV19 exon 15957 16024 ...
		gc= new GenomicCoords("chrM:15957-17259", null, null);
		tif= new TrackIntervalFeature("test_data/hg19.gencode_genes_v19.gtf.gz", gc);
		newGc= tif.coordsOfNextFeature(gc, true);
		assertEquals(15889, (int)newGc.getFrom()); // chrM hg19_wgEncodeGencodeBasicV19 exon 15889 15954
		
		gc= new GenomicCoords("chrM:14672-14898", null, null);
		tif= new TrackIntervalFeature("test_data/hg19.gencode_genes_v19.gtf.gz", gc);
		newGc= tif.coordsOfNextFeature(gc, true);
		assertEquals(14150, (int)newGc.getFrom());

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
	public void canPrintRawLines() throws ClassNotFoundException, IOException, InvalidGenomicCoordsException, InvalidRecordException, SQLException{
		GenomicCoords gc= new GenomicCoords("chr1:1-40000", null, null);
		TrackIntervalFeature tif= new TrackIntervalFeature("test_data/hg19_genes_head.gtf", gc);
		tif.setPrintMode(PrintRawLine.CLIP);

		tif.setPrintRawLineCount(-1);
		assertEquals(20, tif.printFeaturesToFile().split("\n").length);

		tif.setPrintRawLineCount(5);
		assertEquals(5 + 1, tif.printFeaturesToFile().split("\n").length); // +1 for the string of omitted count.
		
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

	// @Test
	public void canConstructFromUcscGenePred() throws InvalidGenomicCoordsException, IOException, ClassNotFoundException, InvalidRecordException, SQLException{

		// From database
		TrackIntervalFeature tif= new TrackIntervalFeature("dm6:refGene", new GenomicCoords("chr3R", null, null));		
		List<IntervalFeature> xset = tif.getFeaturesInInterval("chr3R", 1, 6000000);
		assertTrue(xset.size() > 1000);		
		
		// From local file
		tif= new TrackIntervalFeature("test_data/refGene.hg19.chr7.txt.gz", new GenomicCoords("chr7", null, null));
		xset = tif.getFeaturesInInterval("chr7", 5000000, 6000000);
		assertTrue(xset.size() > 100);
		
		// From direct url connection
		tif= new TrackIntervalFeature("http://hgdownload.soe.ucsc.edu/goldenPath/dm6/database/refGene.txt.gz", new GenomicCoords("chr3R", null, null));		
		xset = tif.getFeaturesInInterval("chr3R", 1, 6000000);
		assertTrue(xset.size() > 1000);		
	}
	
}
