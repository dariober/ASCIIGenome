package samTextViewer;

import static org.junit.Assert.*;
import static org.hamcrest.CoreMatchers.*;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Pattern;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.filter.AlignedFilter;
import htsjdk.samtools.filter.MappingQualityFilter;
import htsjdk.samtools.filter.SamRecordFilter;
import tracks.IntervalFeature;
import tracks.TrackCoverage;
import tracks.TrackFormat;

import org.junit.Test;

import com.google.common.base.Joiner;

import exceptions.InvalidCommandLineException;
import exceptions.InvalidGenomicCoordsException;
import exceptions.InvalidRecordException;
import filter.FirstOfPairFilter;
import filter.FlagToFilter;
import filter.ReadNegativeStrandFilter;

public class UtilsTest {

	static SamReaderFactory srf=SamReaderFactory.make();
	static SamReader samReader= srf.open(new File("test_data/ds051.short.bam"));
	public static SAMSequenceDictionary samSeqDict= samReader.getFileHeader().getSequenceDictionary();
	
	public static String fastaFile= "test_data/chr7.fa";

	@Test
	public void canGetRangeOfListOfValues(){
		List<Double> y= new ArrayList<Double>();
		y.add(1.0);
		y.add(10.0);
		y.add(3.0);
		y.add(Double.NaN);
		assertEquals(1.0, Utils.range(y)[0], 0.0001); // Min
		assertEquals(10.0, Utils.range(y)[1], 0.0001); // Max
		
		// Only NaN
		List<Double> nan= new ArrayList<Double>();
		nan.add(Double.NaN);
		nan.add(Double.NaN);
		nan.add(Double.NaN);
		assertTrue(Utils.range(nan)[0].isNaN());
		assertTrue(Utils.range(nan)[1].isNaN());
		
		// Length of one
		List<Double> y1= new ArrayList<Double>();
		y1.add(1.0);
		assertEquals(1.0, Utils.range(y1)[0], 0.0001);
		assertEquals(1.0, Utils.range(y1)[1], 0.0001);
		
		// Zero length
		List<Double> y0= new ArrayList<Double>();
		assertTrue(Utils.range(y0)[0].isNaN());
		assertTrue(Utils.range(y0)[1].isNaN());
		
	}
	
	@Test
	public void testSortByValueReverse(){
		
		Map<Character, Integer> baseCount= new LinkedHashMap<Character, Integer>();
		baseCount.put('A', 0);
		baseCount.put('C', 0);
		baseCount.put('G', 0);
		baseCount.put('T', 0);
		baseCount.put('N', 0);

		for(int i= 0; i < 10000; i++){
			int count= baseCount.get('G') + 1; 
			baseCount.put('G', count);
		}
		assertEquals('G', (char)Utils.sortByValue(baseCount).keySet().iterator().next());
	} 
	
	// @Test // Not run
	public void throwsGentleMessageOnMissingFaIndex(){
		String fastaFile= "test_data/noindex.fa";
		Utils.checkFasta(fastaFile);
	} 
	
	@Test
	public void canMergeIntervals() throws InvalidGenomicCoordsException{
		
		// Zero len list
		List<IntervalFeature> intv= new ArrayList<IntervalFeature>();
		assertEquals(0, Utils.mergeIntervalFeatures(intv).size());
		
		
		/* MEME: Start of bed features must be augmented by 1 */
		// One feature
		intv.add(new IntervalFeature("chr1 0 10 x1".replaceAll(" ", "\t"), TrackFormat.BED));
		assertEquals(1, Utils.mergeIntervalFeatures(intv).get(0).getFrom());
		// Test the name is taken from the original feature since only one interval is merged (i.e. no merging at all)
		assertEquals(intv.get(0).getName(), Utils.mergeIntervalFeatures(intv).get(0).getName());
		
		// One feature overalapping
		intv.add(new IntervalFeature("chr1 5 10".replaceAll(" ", "\t"), TrackFormat.BED));
		IntervalFeature expected= new IntervalFeature("chr1 0 10".replaceAll(" ", "\t"), TrackFormat.BED);
		
		assertEquals(expected.getFrom(), Utils.mergeIntervalFeatures(intv).get(0).getFrom());
		assertTrue(expected.equals(Utils.mergeIntervalFeatures(intv).get(0)));

		intv.add(new IntervalFeature("chr1 20 100".replaceAll(" ", "\t"), TrackFormat.BED));
		assertEquals(2, Utils.mergeIntervalFeatures(intv).size());
		assertEquals(21, Utils.mergeIntervalFeatures(intv).get(1).getFrom());
		assertEquals(100, Utils.mergeIntervalFeatures(intv).get(1).getTo());
		
		intv.add(new IntervalFeature("chr1 30 110".replaceAll(" ", "\t"), TrackFormat.BED));
		intv.add(new IntervalFeature("chr1 50 110".replaceAll(" ", "\t"), TrackFormat.BED));
		assertEquals(2, Utils.mergeIntervalFeatures(intv).size());
		assertEquals(21, Utils.mergeIntervalFeatures(intv).get(1).getFrom());
		assertEquals(110, Utils.mergeIntervalFeatures(intv).get(1).getTo());
		
		// Touching features get merged into a single one
		intv.clear();
		intv.add(new IntervalFeature("chr1 0 10".replaceAll(" ", "\t"), TrackFormat.BED));
		intv.add(new IntervalFeature("chr1 10 20".replaceAll(" ", "\t"), TrackFormat.BED));
		assertEquals(1, Utils.mergeIntervalFeatures(intv).size());
		assertEquals(1, Utils.mergeIntervalFeatures(intv).get(0).getFrom());
		assertEquals(20, Utils.mergeIntervalFeatures(intv).get(0).getTo());

		// Touching GFF feature 
		intv.clear();
		intv.add(new IntervalFeature("chr1 . . 1 10 . . .".replaceAll(" ", "\t"), TrackFormat.GFF));
		intv.add(new IntervalFeature("chr1 . . 11 20 . . .".replaceAll(" ", "\t"), TrackFormat.GFF));
		assertEquals(1, Utils.mergeIntervalFeatures(intv).size());
		assertEquals(1, Utils.mergeIntervalFeatures(intv).get(0).getFrom());
		assertEquals(20, Utils.mergeIntervalFeatures(intv).get(0).getTo());

		// Nothing to merge 
		intv.clear();
		intv.add(new IntervalFeature("chr1 . . 1 10 . . .".replaceAll(" ", "\t"), TrackFormat.GFF));
		intv.add(new IntervalFeature("chr1 . . 20 30 . . .".replaceAll(" ", "\t"), TrackFormat.GFF));
		intv.add(new IntervalFeature("chr1 . . 40 50 . . .".replaceAll(" ", "\t"), TrackFormat.GFF));
		assertEquals(3, Utils.mergeIntervalFeatures(intv).size());
		
		intv.add(new IntervalFeature("chr1 . . 40 50 . . .".replaceAll(" ", "\t"), TrackFormat.GFF));
		assertEquals(3, Utils.mergeIntervalFeatures(intv).size());
	
	}
	
	@Test
	public void testStringContainsRegex(){
		String x= "foobarbaz";
		String regex= "b.r";
		System.out.println("PATTERN:" + Pattern.compile(regex).matcher(x).find());
	}
	
	@Test
	public void canParseGoToRegion(){
		assertEquals("1-1000", Utils.parseGoToRegion("1-1000"));
		assertEquals("1-1000", Utils.parseGoToRegion("1 - 1,000  "));
		assertEquals("1000", Utils.parseGoToRegion("1,000  "));
		assertEquals("1000-1400", Utils.parseGoToRegion("1000 1200 1300 1400"));
		assertEquals("1000-1400", Utils.parseGoToRegion("1k 1200 1300 1.4k"));
		assertEquals("1000-1400000", Utils.parseGoToRegion(" **1k*******1.4M** "));
		assertEquals("1000-5000000", Utils.parseGoToRegion("--**1k*******1.4M**------5.0M--"));
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
	public void canGetBamReadCount(){
		assertEquals(15098, Utils.getAlignedReadCount(new File("test_data/ds051.actb.bam")));
	}

	@Test
	public void canCountReadsInWindow2() throws InvalidGenomicCoordsException, IOException{
		GenomicCoords gc= new GenomicCoords("chr7:5524838-5611878", samSeqDict, fastaFile);
		List<SamRecordFilter> filters= new ArrayList<SamRecordFilter>();
		
		assertEquals(100377, Utils.countReadsInWindow("test_data/ear045.oxBS.actb.bam", gc, filters));

		filters.add(new AlignedFilter(true));
		assertEquals(100265, Utils.countReadsInWindow("test_data/ear045.oxBS.actb.bam", gc, filters));

		filters= new ArrayList<SamRecordFilter>();
		filters.add(new AlignedFilter(true));
		filters.add(new ReadNegativeStrandFilter(false));
		assertEquals(50157, Utils.countReadsInWindow("test_data/ear045.oxBS.actb.bam", gc, filters));

		filters= new ArrayList<SamRecordFilter>();
		filters.add(new AlignedFilter(true));
		filters.add(new ReadNegativeStrandFilter(true));
		assertEquals(50108, Utils.countReadsInWindow("test_data/ear045.oxBS.actb.bam", gc, filters));

		filters= FlagToFilter.flagToFilterList(80, 1026);
		assertEquals(2729, Utils.countReadsInWindow("test_data/ear045.oxBS.actb.bam", gc, filters));

		filters= FlagToFilter.flagToFilterList(80, 1026);
		filters.add(new MappingQualityFilter(30));
		assertEquals(1592, Utils.countReadsInWindow("test_data/ear045.oxBS.actb.bam", gc, filters));

	}
	
	
	@Test
	public void canCountReadsInWindow() throws InvalidGenomicCoordsException, IOException{
		GenomicCoords gc= new GenomicCoords("chr7:5522436-5613572", samSeqDict, fastaFile);
		List<SamRecordFilter> filters= new ArrayList<SamRecordFilter>();
		
		filters.add(new MappingQualityFilter(30)); // Same as   
		filters.add(new FirstOfPairFilter(true));  // samtools view -q 30 -f 64
		
		long t0= System.currentTimeMillis();
		for(int i= 0; i < 10; i++){
			assertEquals(42770, Utils.countReadsInWindow("test_data/ear045.oxBS.actb.bam", gc, filters));
		}
		long t1= System.currentTimeMillis();
		System.out.println("TIME TO FILTER: " + (t1-t0));
		
		gc= new GenomicCoords("chr7:5524838-5611878", samSeqDict, fastaFile);
		
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
	public void canInitRegion() throws IOException, InvalidGenomicCoordsException{
		assertEquals("chrM", Utils.initRegionFromFile("test_data/ds051.short.bam"));
		assertEquals("chr9", Utils.initRegionFromFile("test_data/hg18_var_sample.wig.v2.1.30.tdf"));
		assertEquals("chr1", Utils.initRegionFromFile("/Users/berald01/Downloads/wgEncodeCaltechRnaSeqGm12878R2x75Il400SigRep2V2.bigWig"));
		assertEquals("chr1:67208779", Utils.initRegionFromFile("test_data/refSeq.hg19.short.bed"));
		assertEquals("chr1:8404074", Utils.initRegionFromFile("test_data/refSeq.hg19.short.sort.bed.gz"));
		assertEquals("chr1:11874", Utils.initRegionFromFile("test_data/hg19_genes_head.gtf.gz"));
	}
	
	@Test
	public void canInitRegionFromURLBam() throws IOException, InvalidGenomicCoordsException{
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
		GenomicCoords gc= new GenomicCoords("chr7:100-200", samSeqDict, fastaFile);

		//String region= Utils.parseConsoleInput("-r chr8:1-1000", gc);
		//assertEquals("chr8:1-1000", region);
		
		//String region= Utils.parseConsoleInput("-r chr8:1", gc);
		//assertEquals("chr8:1", region);

		//region= Utils.parseConsoleInput("-r chr8", gc);
		//assertEquals("chr8", region);
		List<String> tokens= new ArrayList<String>();
		tokens.add("+10");
		String region= Utils.parseConsoleInput(tokens, gc);
		assertEquals("chr7:110-210", region);

		tokens.set(0, "-1000");
		region= Utils.parseConsoleInput(tokens, gc);
		// assertEquals("chr7:110-210", region);
		
	}
	
	@Test
	public void canGetGoToRegionString() throws InvalidGenomicCoordsException, IOException{
		GenomicCoords gc= new GenomicCoords("chr7:100-200", samSeqDict, fastaFile);
		List<String> tokens= new ArrayList<String>();
		tokens.add("1000");
		assertEquals("chr7:1000", Utils.parseConsoleInput(tokens, gc));
		
		tokens.set(0, "1000-10000");
		assertEquals("chr7:1000-10000", Utils.parseConsoleInput(tokens, gc));

		tokens.set(0, "1,000 - 10,000");
		assertEquals("chr7:1000-10000", Utils.parseConsoleInput(tokens, gc));
		
		tokens.set(0, ":foo");
		try{
			System.err.println(Utils.parseConsoleInput(tokens, gc));
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
	public void canAddTracksToList() throws IOException, InvalidGenomicCoordsException, InvalidCommandLineException{
		List<String> inputFileList= new ArrayList<String>();
		inputFileList.add("foo");
		inputFileList.add("bar");
		List<String> newFileNames= new ArrayList<String>();
		newFileNames.add("test_data/ds051.actb.bam");
		newFileNames.add("nonsense");
		Utils.addSourceName(inputFileList, newFileNames);
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
	
	@Test
	public void canSplitStringInTokens(){

		String str= "foo    bar    baz   ";
		assertTrue(Utils.tokenize(str, " ").contains("bar"));
		assertTrue(Utils.tokenize(str, " ").contains("baz"));
		// assertEquals("bar", Utils.tokenize(str, " ").get(1));
		
		
		// Note use of quotes
		ArrayList<String> xx= Utils.tokenize("'foo && bar' "
				+ "&& bar"
				+ "&&baz "
				+ "&& 'foo && biz'"
				+ "&& 'foo && \' biz'", "&&");
		for (String token : xx) {
			System.out.println(token);
		}
		
		assertEquals("foo && bar", xx.get(0));
		assertEquals("bar", xx.get(1));
		assertEquals("baz", xx.get(2));
		assertEquals("foo && biz", xx.get(3));

		xx= Utils.tokenize("gene \"ACTB\"", "&&");
		assertEquals("gene \"ACTB\"", xx.get(0));
		
		// Unusual input:
		assertEquals(null, Utils.tokenize(null, " "));
		assertTrue(Utils.tokenize("", " ").size() == 0);
		assertTrue(Utils.tokenize("   ", " ").size() == 0);
		
		// Reverse token
		String cmdInput= "goto chr1 0 100";
		assertEquals(cmdInput, Joiner.on(" ").join(Utils.tokenize(cmdInput, " ")));
		assertEquals("goto chr1 0 100", Joiner.on(" ").join(Utils.tokenize("   goto   chr1   0   100  ", " ")));
		
	}
	
	@Test
	public void canPrintToStdoutOrFile() throws InvalidGenomicCoordsException, IOException, ClassNotFoundException, InvalidRecordException, SQLException{

		GenomicCoords gc= new GenomicCoords("chr7:5566770-5566870", samSeqDict, fastaFile);
		TrackCoverage tc= new TrackCoverage("test_data/ds051.short.bam", gc, false);
		
		File filename= new File("tmp.txt");
		filename.delete();
		Utils.printer(tc.getTitle(), "tmp.txt");
		Utils.printer(tc.printToScreen(), "tmp.txt");
		// Check tmp.txt looks ok.
	}
	
	@Test
	public void canGetWritableFileOrNull() throws InvalidGenomicCoordsException, IOException{
		GenomicCoords gc= new GenomicCoords("chr7:1-200", samSeqDict, fastaFile);
		String x= Utils.parseCmdinputToGetSnapshotFile("save", gc);
		assertEquals("chr7_1-200.txt", x);
		x= Utils.parseCmdinputToGetSnapshotFile("save /tmp/foo.txt", gc);
		assertEquals("/tmp/foo.txt", x);
	}

	@Test
	public void testPng() throws IOException{
		Utils.convertTextFileToGraphic(new File("test_data/chr7_5564857-5570489.txt"), new File("tmp.png"));
	}
	
	@Test
	public void canConvertCoordsToString(){
		assertEquals("chr1:1-100", Utils.coordinatesToString("chr1", 1, 100));
		assertEquals("chr1:1", Utils.coordinatesToString("chr1", 1, null));
		assertEquals("chr1:1", Utils.coordinatesToString("chr1", 1, null));
		assertEquals("chr1:1", Utils.coordinatesToString("chr1", null, null));
		assertEquals("chr1:1", Utils.coordinatesToString("chr1", null, -1));
		assertEquals("chr1:10", Utils.coordinatesToString("chr1", 10, 9));
	}
	
	
}