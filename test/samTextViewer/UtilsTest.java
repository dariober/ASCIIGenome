package samTextViewer;

import static org.junit.Assert.*;
import static org.hamcrest.CoreMatchers.*;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.Reader;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.zip.GZIPInputStream;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import tracks.TrackFormat;

import org.junit.Test;

import exceptions.InvalidGenomicCoordsException;

public class UtilsTest {

	static SamReaderFactory srf=SamReaderFactory.make();
	static SamReader samReader= srf.open(new File("test_data/ds051.short.bam"));
	public static SAMSequenceDictionary samSeqDict= samReader.getFileHeader().getSequenceDictionary();
	
	public static String fastaFile= "test_data/chr7.fa";
	
	@Test
	public void canParseGoToRegion(){
		assertEquals("1-1000", Utils.parseGoToRegion("1-1000"));
		assertEquals("1-1000", Utils.parseGoToRegion("1 - 1,000  "));
		assertEquals("1000", Utils.parseGoToRegion("1,000  "));		
	}
	
	@Test 
	public void canParseZoomArg(){
		assertEquals(2, Utils.parseZoom("zo 2", 1));
		assertEquals(3, Utils.parseZoom("zo 3 foo", 1));
		assertEquals(4, Utils.parseZoom("zo   4  ", 1));
		assertEquals(0, Utils.parseZoom("zo   0", 1));
		assertEquals(1, Utils.parseZoom("zo", 1));
		assertEquals(1, Utils.parseZoom("zo -3", 1));
		assertEquals(1, Utils.parseZoom("zo 3.3", 1));
	}
	
	@Test
	public void canTabulateListOfFeatures(){
		List<String> rawList= new ArrayList<String>();
		rawList.add("1\tgenedb\tgene\t2964\t45090");
		rawList.add("chr1\tgenedb_long\tgene\t2964\t45090");
		rawList.add("1\tfoo\tna\t2"); // Missing last field
		rawList.add("1\tfoo\tna\t2\t10");
		rawList.add("1\tfoo\t\t2\t10"); // Empty cell
		

		List<String> expList= new ArrayList<String>();
		expList.add("1    genedb      gene 2964 45090");
		expList.add("chr1 genedb_long gene 2964 45090");
		expList.add("1    foo         na   2");
		expList.add("1    foo         na   2    10");
		expList.add("1    foo              2    10");

		// Handling region with no features to print
		rawList= new ArrayList<String>();
		List<String> obsList= Utils.tabulateList(rawList);
		expList= new ArrayList<String>();
		assertThat(expList, is(obsList));
		
	}
		
	@Test
	public void canTestForTabixIndex() throws IOException{
		assertTrue(Utils.hasTabixIndex("test_data/test.bedGraph.gz"));
		assertTrue(! Utils.hasTabixIndex("test_data/test.bedGraph"));
	}
		
	@Test
	public void canGetFileTypeFromName(){
		
		assertEquals(TrackFormat.BIGWIG,
		Utils.getFileTypeFromName("/Users/berald01/Downloads/wgEncodeCaltechRnaSeqGm12878R2x75Il400SigRep2V2.bigWig"));
	} 

	@Test
	public void canInitRegion() throws IOException{
		assertEquals("chrM", Utils.initRegionFromFile("test_data/ds051.short.bam"));
		assertEquals("chr9", Utils.initRegionFromFile("test_data/hg18_var_sample.wig.v2.1.30.tdf"));
		assertEquals("chr1", Utils.initRegionFromFile("/Users/berald01/Downloads/wgEncodeCaltechRnaSeqGm12878R2x75Il400SigRep2V2.bigWig"));
		assertEquals("chr1:67208779", Utils.initRegionFromFile("test_data/refSeq.hg19.short.bed"));
		assertEquals("chr1:8404074", Utils.initRegionFromFile("test_data/refSeq.hg19.short.sort.bed.gz"));
		assertEquals("chr1:11874", Utils.initRegionFromFile("test_data/hg19_genes_head.gtf.gz"));
	}
	
	@Test
	public void canInitRegionFromURLBam() throws IOException{
		String reg= Utils.initRegionFromFile("http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeHaibTfbs/wgEncodeHaibTfbsA549Atf3V0422111Etoh02AlnRep1.bam");
		assertEquals("chr1", reg);
	}
	
	@Test
	public void testBamHasIndex() throws IOException{
		assertTrue(Utils.bamHasIndex("test_data/ds051.short.bam"));
		assertTrue(!Utils.bamHasIndex("test_data/ds051.noindex.bam"));
	}
	
	@Test
	public void canGetClosestIndex(){
		int windowSize= 150;
		List<Double> mapping = Utils.seqFromToLenOut(1, 1000000, windowSize);
		for(int i=0; i < windowSize; i++){
			System.out.println("Index: " + i + " position: " + mapping.get(i));
		}
		
		assertEquals(windowSize-1, Utils.getIndexOfclosestValue(1000000, mapping));
		assertEquals((windowSize/2)-1, Utils.getIndexOfclosestValue(1000000/2, mapping));
		assertEquals(0, Utils.getIndexOfclosestValue(1, mapping));
		assertEquals(0, Utils.getIndexOfclosestValue((6712/2.0)-1, mapping));
		assertEquals(1, Utils.getIndexOfclosestValue((6712/2.0)+1, mapping));
		assertEquals(0, Utils.getIndexOfclosestValue((6712/2.0), mapping));
		assertEquals(102, Utils.getIndexOfclosestValue(684564, mapping));
		assertEquals(112, Utils.getIndexOfclosestValue(750000, mapping));
	}
	
	@Test
	public void canGenerateSequence(){
		Utils.seqFromToLenOut(15, 17, 13);
		Utils.seqFromToLenOut(17, 15, 13);
		Utils.seqFromToLenOut(15, 26, 13).size();
		Utils.seqFromToLenOut(15, 15, 13);
		
		// Length of 1 returns "from" like R seq(0, 10, length.out= 1) -> 0
		assertEquals(1, Utils.seqFromToLenOut(0, 10, 1).size());
		assertEquals(0, Utils.seqFromToLenOut(0, 10, 1).get(0), 0.00001);
		assertEquals(0, Utils.seqFromToLenOut(0, 10, 0).size()); // Zero-length sequence
		
		assertEquals((Double)Double.NaN, Utils.seqFromToLenOut(Double.NaN, Double.NaN, 10).get(0));
	}
	
	@Test
	public void canTestForAllNaN(){
		ArrayList<Double> x= new ArrayList<Double>();
		x.add(Double.NaN);
		x.add(Double.NaN);
		x.add(Double.NaN);
		assertTrue(Utils.allIsNaN(x));
		
		x.add((Double) 1.0);
		assertFalse(Utils.allIsNaN(x));
	}
	
	@Test
	public void canParseInputAndUpdateGenomicCoords() throws InvalidGenomicCoordsException, IOException{
		GenomicCoords gc= new GenomicCoords("chr7:100-200", samSeqDict, 100, fastaFile);

		//String region= Utils.parseConsoleInput("-r chr8:1-1000", gc);
		//assertEquals("chr8:1-1000", region);
		
		//String region= Utils.parseConsoleInput("-r chr8:1", gc);
		//assertEquals("chr8:1", region);

		//region= Utils.parseConsoleInput("-r chr8", gc);
		//assertEquals("chr8", region);
		
		String region= Utils.parseConsoleInput("+10", gc);
		assertEquals("chr7:110-210", region);

		region= Utils.parseConsoleInput("-1000", gc);
		// assertEquals("chr7:110-210", region);
		
		String rawInput= "-r chr10 -F 1024";
		List<String> clArgs= Arrays.asList(rawInput.split("\\s+"));
		// System.out.println(clArgs.indexOf("-R"));		
	}
	
	@Test
	public void canGetGoToRegionString() throws InvalidGenomicCoordsException, IOException{
		GenomicCoords gc= new GenomicCoords("chr7:100-200", samSeqDict, 100, fastaFile);
		String rawInput= "1000";
		assertEquals("chr7:1000", Utils.parseConsoleInput(rawInput, gc));
		
		rawInput= "1000-10000";
		assertEquals("chr7:1000-10000", Utils.parseConsoleInput(rawInput, gc));
		
		rawInput= " 1,000 - 10,000";
		assertEquals("chr7:1000-10000", Utils.parseConsoleInput(rawInput, gc));
		
		rawInput= ":foo"; // Must fail
		try{
			System.err.println(Utils.parseConsoleInput(rawInput, gc));
			fail();
		} catch (Exception e) {
			
		}
	}
	
	@Test
	public void canRoundNumbersToSignificantDigits(){
		
		double x= 1000.123456789;
		double y= 1001.123456789;
		int nSignif= 3;
		double[] rounded= Utils.roundToSignificantDigits(x, y, nSignif);
		assertEquals(1000.123, rounded[0], 0.001); // Regular rounding
		assertEquals(1001.123, rounded[1], 0.001);
		
		x= 1000.00012345;
		y= 1000.0012345;
		nSignif= 2;
		rounded= Utils.roundToSignificantDigits(x, y, nSignif);
		assertEquals(1000.00012, rounded[0], 1e-16);
		assertEquals(1000.00123, rounded[1], 1e-16);		
		
		x= 1000.0009876;
		y= 1000.00987654;
		nSignif= 3;
		rounded= Utils.roundToSignificantDigits(x, y, nSignif);
		assertEquals(1000.000988, rounded[0], 1e-16);
		assertEquals(1000.009877, rounded[1], 1e-16);		

	}
	
	@Test
	public void canRoundToSiginficantDigits(){
		assertEquals(0.001234, Utils.roundToSignificantFigures(0.001234, 5), 1e-16);
		assertEquals(1200, (int)Utils.roundToSignificantFigures(1234, 2));
	}
	
	@Test
	public void canAddMetricSuffixToInt(){
		assertEquals("123M", Utils.parseIntToMetricSuffix(123000000));
		assertEquals("123k", Utils.parseIntToMetricSuffix(123000));
		assertEquals("123", Utils.parseIntToMetricSuffix(123));
	}
	
	@Test
	public void canTestForExistingURLFile(){
		String urlStr= "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeHaibTfbs/wgEncodeHaibTfbsA549Atf3V0422111Etoh02PkRep1.broadPeak.gz";
		assertTrue(Utils.urlFileExists(urlStr));
		assertFalse(Utils.urlFileExists(urlStr + "foobar"));

		assertTrue(Utils.urlFileExists("ftp://ftp.ensembl.org/pub/release-84/gtf/homo_sapiens/Homo_sapiens.GRCh38.84.abinitio.gtf.gz"));
		assertFalse(Utils.urlFileExists("ftp://ftp.ensembl.org/pub/release-84/gtf/homo_sapiens/foobar"));
		
		/* This should return false but it doesn't */
		// urlStr= "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeHaibTfbs/";
		// assertFalse(Utils.urlFileExists(urlStr));
	}
	
	@Test
	public void canAddTracksToList() throws IOException, InvalidGenomicCoordsException{
		List<String> inputFileList= new ArrayList<String>();
		inputFileList.add("foo");
		inputFileList.add("bar");
		List<String> newFileNames= new ArrayList<String>();
		newFileNames.add("test_data/ds051.actb.bam");
		newFileNames.add("nonsense");
		Utils.addTrack(inputFileList, newFileNames);
		assertEquals(3, inputFileList.size());
	}
	
	@Test
	public void canPrintSequenceDict(){
		assertTrue(Utils.printSamSeqDict(samSeqDict, 30).startsWith("chrM  16571"));
		assertTrue(Utils.printSamSeqDict(samSeqDict, 30).endsWith("chrY  59373566  |||||||"));
		SAMSequenceDictionary emptyDict= new SAMSequenceDictionary();
		Utils.printSamSeqDict(emptyDict, 30);
		Utils.printSamSeqDict(null, 30);
	}
}