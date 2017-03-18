package samTextViewer;

import static org.hamcrest.CoreMatchers.is;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertThat;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Pattern;

import org.apache.commons.io.FileUtils;
import org.junit.Test;

import com.google.common.base.Joiner;

import exceptions.InvalidCommandLineException;
import exceptions.InvalidGenomicCoordsException;
import exceptions.InvalidRecordException;
import faidx.UnindexableFastaFileException;
import filter.FirstOfPairFilter;
import filter.FlagToFilter;
import filter.ReadNegativeStrandFilter;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.filter.AlignedFilter;
import htsjdk.samtools.filter.MappingQualityFilter;
import htsjdk.samtools.filter.SamRecordFilter;
import jline.console.ConsoleReader;
import jline.console.history.History;
import jline.console.history.MemoryHistory;
import tracks.IntervalFeature;
import tracks.TrackFormat;

public class UtilsTest {

	static SamReaderFactory srf=SamReaderFactory.make();
	static SamReader samReader= srf.open(new File("test_data/ds051.short.bam"));
	public static SAMSequenceDictionary samSeqDict= samReader.getFileHeader().getSequenceDictionary();
	
	public static String fastaFile= "test_data/chr7.fa";
	
	@Test
	public void canPadMultilineString(){
		
		assertEquals("foo  ", Utils.padEndMultiLine("foo", 5));
		
		// Empty string is returned as is.
		assertEquals("", Utils.padEndMultiLine("", 3));
		
		// One newline is expanded to TWO strings 
		assertEquals("   \n   ", Utils.padEndMultiLine("\n", 3));
		
		// Note starting with emtpy line, which is going to be padded.
		String x= "\nfoo\n1234567890";
		String padded= Utils.padEndMultiLine(x, 5);
		
		String[] p = padded.split("\n");
		assertEquals("     ", p[0]);
		assertEquals("foo  ", p[1]);
		assertEquals("1234567890", p[2]);
	}
	
	@Test
	public void canFilterUsingAwk() throws IOException{
		// Note single quotes around the awk script
		assertTrue(Utils.passAwkFilter("chr1\t10\t100", "-v VAR=5 '$2 > VAR && $1'"));
		
		// Using single quotes inside awk is tricky. 
		// See https://www.gnu.org/software/gawk/manual/html_node/Quoting.html 
		// and see http://stackoverflow.com/questions/9899001/how-to-escape-single-quote-in-awk-inside-printf
		// for using '\'' as a single quote
		assertTrue(Utils.passAwkFilter("chr'1", "'$1 == \"chr'\''1\"'"));
		
		assertTrue( ! Utils.passAwkFilter("'chr1\t10\t100", "-v VAR=50 '$2 > VAR'"));
		
		assertTrue( Utils.passAwkFilter("'chr1\t10\t100", "'($3 - $2) > 50'"));
		assertTrue( ! Utils.passAwkFilter("'chr1\t10\t100", "'($3 - $2) > 500'"));
		
		assertTrue(Utils.passAwkFilter("'chr1\t10\t100'", "")); // Empty script equals to no filter.
		assertTrue(Utils.passAwkFilter("'chr1\t10\t100'", "  "));
		
		// Valid awk script but output is not empty and not equal to input. Return NULL
		assertEquals(null, Utils.passAwkFilter("'chr1\t10\t100'", "'{print $1}'"));
		
		// Broken awk script:
		boolean pass= false;
		try{
			assertEquals(null, Utils.passAwkFilter("'chr1\t10\t100'", "'print {'"));
		} catch(IOException e){
			pass= true;
		}
		assertTrue(pass);
	}
	
	@Test
	public void canFilterSamTagWithAwk() throws IOException{
		String rec= "read\t0\tchr7\t5566778\t50\t5M\t*\t0\t0\tCTCAT\tIIIII\tMD:Z:75\tRG:Z:1\tXG:i:0\tNH:i:1\tNM:i:0\tXM:i:0\tXN:i:0\tXO:i:0\tAS:i:0\tYT:Z:UU";
		
		//Filter for NH tag value
		assertTrue(Utils.passAwkFilter(rec, "'getSamTag(\"NH\") > 0'"));
		assertFalse(Utils.passAwkFilter(rec, "'getSamTag(\"NH\") > 10'"));
		
		// Missing tag
		assertFalse(Utils.passAwkFilter(rec, "'getSamTag(\"ZZ\") > 0'"));
		// MIssing tag searched but not used
		assertTrue(Utils.passAwkFilter(rec, "'{getSamTag(\"ZZ\"); print $0'}"));
		
		long t0= System.currentTimeMillis();
		int i= 0;
		while(i < 1000){
			Utils.passAwkFilter(rec, "'getSamTag(\"NH\") > 0'");
			i++;
		}
		long t1= System.currentTimeMillis();
		assertTrue((t1-t0) < 3000); // It can filter reasonably fast (?) 
	}
	
	@Test
	public void canCheckForUpdates() throws IOException{
		
		List<String> up= Utils.checkUpdates(50000);
			
		assertEquals(2, up.size());
		assertTrue(Character.isDigit(up.get(0).charAt(0)));
		assertTrue(Character.isDigit(up.get(1).charAt(0)));
		
		assertEquals(0, Utils.versionCompare("1.0.0", "1.0.0"));
		assertEquals(1, Utils.versionCompare("1.0.0", "0.0.9")); // Running version ahead of repo

		assertEquals(-1, Utils.versionCompare("1.0.0", "1.0.1")); // Running version out of date
		assertEquals(-1, Utils.versionCompare("1.0.0", "1.0.0.1"));
		
		// This should timeput and throw a warning
		up= Utils.checkUpdates(1);
		
	}
	
	@Test
	public void testAwk() throws Exception{

		String fin= "/Users/berald01/Downloads/hg19.gencode_genes_v19.gtf"; // "test_data/hg19_genes_head.gtf"
		InputStream is= new FileInputStream(fin);
		
		String[] args= {"-F", "\t", "NR <= 100000"}; 
		System.err.println(Arrays.toString(args));
		
		PrintStream stdout = System.out;
		ByteArrayOutputStream baos = new ByteArrayOutputStream();
		try{
			PrintStream os= new PrintStream(baos);
			System.err.println("START");
			long t0= System.currentTimeMillis();
			new org.jawk.Main(args, is, os, System.err);
			long t1= System.currentTimeMillis();
			System.err.println(t1-t0);

		} finally{
			System.setOut(stdout);
			is.close();
		}
		System.err.println("DONE");
		long t0= System.currentTimeMillis();
		String content = new String(baos.toByteArray(), StandardCharsets.UTF_8);
		// assertEquals(10, content.split("\n").length);
		
		List<String> full= FileUtils.readLines(new File(fin));
		Set<String> keepme= new HashSet<String>(Arrays.asList(content.split("\n")));
		for(String line : full){
			if(keepme.contains(line)){
				// System.err.println(line);
			}
		}
		long t1= System.currentTimeMillis();
		System.err.println(t1-t0);
 
	}
	
	@Test
	public void testHistory() throws IOException{
		ConsoleReader console= new ConsoleReader();
		//History history= new History(new File(System.getProperty("user.home") + File.separator + ".asciigenome_history"));
		History history= new MemoryHistory();
		history.add("foobar");
		history.add("baz");
		console.setHistory(history);
		System.out.println(console.getHistory());
	}
	
	@Test
	public void canGlobFiles() throws IOException{

		ArrayList<String> cmdInput = Utils.tokenize("test_data/ear*{bam,tdf} README.*", " ");
		List<String> globbed= Utils.globFiles(cmdInput);
		assertEquals(3, globbed.size());
		
		cmdInput = Utils.tokenize("test_data", " ");
		globbed= Utils.globFiles(cmdInput);
		assertTrue(globbed.size() > 10);
		
		cmdInput = Utils.tokenize("test_data/*", " ");
		globbed= Utils.globFiles(cmdInput);
		assertTrue(globbed.size() > 10);
		
		cmdInput = Utils.tokenize("test_data//*", " ");
		globbed= Utils.globFiles(cmdInput);
		assertTrue(globbed.size() > 10);
		
		cmdInput = Utils.tokenize("test_data/../test_data", " ");
		globbed= Utils.globFiles(cmdInput);
		assertTrue(globbed.size() > 10);
		
		cmdInput = Utils.tokenize("test_data/*.gtf.*", " ");
		globbed= Utils.globFiles(cmdInput);
		System.out.println(globbed);
		
		// With URLs
		cmdInput = Utils.tokenize("ftp://ftp.ensembl.org/pub/release-84/gtf/homo_sapiens/Homo_sapiens.GRCh38.84.abinitio.gtf.gz", " ");
		globbed= Utils.globFiles(cmdInput);
		assertEquals(1, globbed.size());
	}
	
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
	
	@Test
	public void createFastaIndex() throws IOException, UnindexableFastaFileException{
		String fastaFile= "test_data/noindex.fa";
		if(new File("test_data/noindex.fa.fai").isFile()){
			new File("test_data/noindex.fa.fai").delete();
		}
		
		Utils.checkFasta(fastaFile);
		
		assertTrue( (new File("test_data/noindex.fa.fai")).isFile() );
		assertTrue( (new File("test_data/noindex.fa.fai")).length() > 10 );
	} 
	
	@Test
	public void canMergeIntervals() throws InvalidGenomicCoordsException{
		
		// Zero len list
		List<IntervalFeature> intv= new ArrayList<IntervalFeature>();
		assertEquals(0, Utils.mergeIntervalFeatures(intv, false).size());
		
		// Fully contained feature
		intv.clear();
		intv.add(new IntervalFeature("chr1 . . 100 1000 . . .".replaceAll(" ", "\t"), TrackFormat.GTF));
		intv.add(new IntervalFeature("chr1 . . 200 300 . . .".replaceAll(" ", "\t"), TrackFormat.GTF));
		assertEquals(1, Utils.mergeIntervalFeatures(intv, false).size());
		assertEquals(100, Utils.mergeIntervalFeatures(intv, false).get(0).getFrom());
		assertEquals(1000, Utils.mergeIntervalFeatures(intv, false).get(0).getTo());
		
		// Partial overlap contained feature
		intv.clear();
		intv.add(new IntervalFeature("chr1 . . 100 1000 . . .".replaceAll(" ", "\t"), TrackFormat.GTF));
		intv.add(new IntervalFeature("chr1 . . 200 300 . . .".replaceAll(" ", "\t"), TrackFormat.GTF));
		intv.add(new IntervalFeature("chr1 . . 500 5000 . . .".replaceAll(" ", "\t"), TrackFormat.GTF));
		assertEquals(1, Utils.mergeIntervalFeatures(intv, false).size());
		assertEquals(100, Utils.mergeIntervalFeatures(intv, false).get(0).getFrom());
		assertEquals(5000, Utils.mergeIntervalFeatures(intv, false).get(0).getTo());
		
		/* MEMO: Start of bed features must be augmented by 1 */
		// One feature
		intv.clear();
		intv.add(new IntervalFeature("chr1 0 10 x1".replaceAll(" ", "\t"), TrackFormat.BED));
		assertEquals(1, Utils.mergeIntervalFeatures(intv, false).get(0).getFrom());
		// Test the name is taken from the original feature since only one interval is merged (i.e. no merging at all)
		assertEquals(intv.get(0).getName(), Utils.mergeIntervalFeatures(intv, false).get(0).getName());
		
		// One feature overalapping
		intv.add(new IntervalFeature("chr1 5 10".replaceAll(" ", "\t"), TrackFormat.BED));
		IntervalFeature expected= new IntervalFeature("chr1 0 10".replaceAll(" ", "\t"), TrackFormat.BED);
		
		assertEquals(expected.getFrom(), Utils.mergeIntervalFeatures(intv, false).get(0).getFrom());
		assertTrue(expected.equals(Utils.mergeIntervalFeatures(intv, false).get(0)));

		intv.add(new IntervalFeature("chr1 20 100".replaceAll(" ", "\t"), TrackFormat.BED));
		assertEquals(2, Utils.mergeIntervalFeatures(intv, false).size());
		assertEquals(21, Utils.mergeIntervalFeatures(intv, false).get(1).getFrom());
		assertEquals(100, Utils.mergeIntervalFeatures(intv, false).get(1).getTo());
		
		intv.add(new IntervalFeature("chr1 30 110".replaceAll(" ", "\t"), TrackFormat.BED));
		intv.add(new IntervalFeature("chr1 50 110".replaceAll(" ", "\t"), TrackFormat.BED));
		assertEquals(2, Utils.mergeIntervalFeatures(intv, false).size());
		assertEquals(21, Utils.mergeIntervalFeatures(intv, false).get(1).getFrom());
		assertEquals(110, Utils.mergeIntervalFeatures(intv, false).get(1).getTo());
		
		// Touching features get merged into a single one
		intv.clear();
		intv.add(new IntervalFeature("chr1 0 10".replaceAll(" ", "\t"), TrackFormat.BED));
		intv.add(new IntervalFeature("chr1 10 20".replaceAll(" ", "\t"), TrackFormat.BED));
		assertEquals(1, Utils.mergeIntervalFeatures(intv, false).size());
		assertEquals(1, Utils.mergeIntervalFeatures(intv, false).get(0).getFrom());
		assertEquals(20, Utils.mergeIntervalFeatures(intv, false).get(0).getTo());

		// Touching GFF feature 
		intv.clear();
		intv.add(new IntervalFeature("chr1 . . 1 10 . . .".replaceAll(" ", "\t"), TrackFormat.GTF));
		intv.add(new IntervalFeature("chr1 . . 11 20 . . .".replaceAll(" ", "\t"), TrackFormat.GTF));
		assertEquals(1, Utils.mergeIntervalFeatures(intv, false).size());
		assertEquals(1, Utils.mergeIntervalFeatures(intv, false).get(0).getFrom());
		assertEquals(20, Utils.mergeIntervalFeatures(intv, false).get(0).getTo());

		// Nothing to merge 
		intv.clear();
		intv.add(new IntervalFeature("chr1 . . 1 10 . . .".replaceAll(" ", "\t"), TrackFormat.GTF));
		intv.add(new IntervalFeature("chr1 . . 20 30 . . .".replaceAll(" ", "\t"), TrackFormat.GTF));
		intv.add(new IntervalFeature("chr1 . . 40 50 . . .".replaceAll(" ", "\t"), TrackFormat.GTF));
		assertEquals(3, Utils.mergeIntervalFeatures(intv, false).size());
		
		intv.add(new IntervalFeature("chr1 . . 40 50 . . .".replaceAll(" ", "\t"), TrackFormat.GTF));
		assertEquals(3, Utils.mergeIntervalFeatures(intv, false).size());
		
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
		assertEquals(0, Utils.parseZoom("zo -3", 1));   // < 0 reset to 0
		assertEquals(0, Utils.parseZoom("zo 3.3", 1)); // Invalid INT reset to zero
		assertEquals(0, Utils.parseZoom("zo foo", 1)); // Invalid INT reset to zero
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
	public void canGetBamReadCount() throws IOException{
		assertEquals(15098, Utils.getAlignedReadCount("test_data/ds051.actb.bam"));
		// Painfully slow!
		// assertEquals(6337212, Utils.getAlignedReadCount("http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeCaltechRnaSeq/wgEncodeCaltechRnaSeqGm12878R2x75Il400SplicesRep2V2.bam"));
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

		// This ftp file has index but see https://github.com/samtools/htsjdk/issues/797
		assertTrue( ! Utils.hasTabixIndex("ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/analysis_results/input_call_sets/ALL.wex.union_illumina_wcmc_bcm_bc_bi.20110521.snps.exome.sites.vcf.gz"));
		
		// HTTP is ok.
		// If this file does not exist, put any valid tabix file and its index on Dropbox/Public and use
		// the dropbox link here.
		assertTrue(Utils.hasTabixIndex("http://genome.ucsc.edu/goldenPath/help/examples/vcfExample.vcf.gz"));

		// NB: Uncompressed files give a OutOfMemoryError: Java heap space
		assertTrue(! Utils.hasTabixIndex("ftp://ftp.solgenomics.net/tomato_genome/annotation/ITAG3.0_release/ITAG3.0_RepeatModeler_repeats_light.gff"));
	}
		
	@Test
	public void canGetFileTypeFromName(){
		
		assertEquals(TrackFormat.VCF,
		Utils.getFileTypeFromName("test/gz.vcf.bgz"));
		
		assertEquals(TrackFormat.BIGWIG,
		Utils.getFileTypeFromName("http://foo/bar/wgEncodeCaltechRnaSeqGm12878R2x75Il400SigRep2V2.bigWig"));
	} 

	@Test
	public void canInitRegion() throws IOException, InvalidGenomicCoordsException, ClassNotFoundException, InvalidCommandLineException, InvalidRecordException, SQLException{
		
		assertEquals("chrM", Utils.initRegionFromFile("test_data/ds051.short.bam"));
		assertEquals("chr9", Utils.initRegionFromFile("test_data/hg18_var_sample.wig.v2.1.30.tdf"));
		assertEquals("chr1:10536", Utils.initRegionFromFile("test_data/wgEncodeCaltechRnaSeqGm12878R2x75Il400SigRep2V2.bigWig"));
		assertEquals("chr1:67208779", Utils.initRegionFromFile("test_data/refSeq.hg19.short.bed"));
		assertEquals("chr1:8404074", Utils.initRegionFromFile("test_data/refSeq.hg19.short.sort.bed.gz"));
		assertEquals("chr1:11874", Utils.initRegionFromFile("test_data/hg19_genes_head.gtf.gz"));
		assertEquals("chr1:564666", Utils.initRegionFromFile("test_data/wgEncodeDukeDnase8988T.fdr01peaks.hg19.bb"));
		assertEquals("chr1", Utils.initRegionFromFile("http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeCaltechRnaSeq/wgEncodeCaltechRnaSeqGm12878R2x75Il400SplicesRep2V2.bam"));
		
		boolean pass= false;
		try{
			Utils.initRegionFromFile(
					"http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeCaltechRnaSeq/wgEncodeCaltechRnaSeqGm12878R2x75Il400SigRep2V2.bigWig");
		} catch(InvalidGenomicCoordsException e){
			pass= true;
		}
		assertTrue(pass);
		
		// assertTrue(Utils.initRegionFromFile("hg19:refGene").startsWith("chr1:"));
	}
	
	@Test
	public void canInitRegionFromURLBam() throws IOException, InvalidGenomicCoordsException, ClassNotFoundException, InvalidCommandLineException, InvalidRecordException, SQLException{
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
		Double[] rounded= Utils.roundToSignificantDigits(x, y, nSignif);
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
		
		// Empty dict -> Empty string
		assertEquals("", Utils.printSamSeqDict(emptyDict, 30));
		assertEquals("", Utils.printSamSeqDict(null, 30));
	}
	
//	@Test
//	public void canSplitCommandLineInTokens(){
//		
//		List<String> x= Lists.newArrayList(Splitter.on(Pattern.compile("&&(?=([^']*'[^']*')*[^']*$)"))
//			.trimResults()
//			.split("awk '$4 > 910000 && $4 < 1000000' && grep -i foo"));
//		
//		System.err.println(x);
//		
//		List<String> cmdInput= Utils.tokenize("awk '$4 > 910000 && $4 < 1000000'", "&&");
//		System.err.println(cmdInput);
//		assertEquals(1, cmdInput.size());
//		
//		// Note single quotes
//		cmdInput= Utils.tokenize("foo && bar && 'baz && baz'", "&&");
//		assertEquals(3, cmdInput.size());
//		assertEquals("baz && baz", cmdInput.get(2));
//	}
	
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
	public void canGetWritableFileOrNull() throws InvalidGenomicCoordsException, IOException{
		GenomicCoords gc= new GenomicCoords("chr7:1-200", samSeqDict, fastaFile);
		String x= Utils.parseCmdinputToGetSnapshotFile("save", gc);
		assertEquals("chr7_1-200.txt", x);
		x= Utils.parseCmdinputToGetSnapshotFile("save /tmp/foo.txt", gc);
		assertEquals("/tmp/foo.txt", x);
	}

//	@Test
//	public void testPng() throws IOException{
//		Utils.convertTextFileToGraphic(new File("test_data/chr7_5564857-5570489.txt"), new File("tmp.png"));
//		new File("tmp.png").deleteOnExit();
//	}
	
	@Test
	public void canConvertCoordsToString(){
		assertEquals("chr1:1-100", Utils.coordinatesToString("chr1", 1, 100));
		assertEquals("chr1:1", Utils.coordinatesToString("chr1", 1, null));
		assertEquals("chr1:1", Utils.coordinatesToString("chr1", 1, null));
		assertEquals("chr1:1", Utils.coordinatesToString("chr1", null, null));
		assertEquals("chr1:1", Utils.coordinatesToString("chr1", null, -1));
		assertEquals("chr1:10", Utils.coordinatesToString("chr1", 10, 9));
	}
	
	@Test
	public void canTestForUcscSource(){
		assertTrue(Utils.isUcscGenePredSource("dm6:refGene"));
		assertTrue(Utils.isUcscGenePredSource("test_data/refGene.hg19.chr7.txt.gz"));
		assertTrue(Utils.isUcscGenePredSource("http://hgdownload.soe.ucsc.edu/goldenPath/dm6/database/refGene.txt.gz"));
		assertTrue(!Utils.isUcscGenePredSource("test_data/hg19_genes.gtf.gz"));
	}
	
	@Test
	public void canExpandTildeToHomeDir(){
		
		//No change
		assertEquals("/foo/bar/baz", Utils.tildeToHomeDir("/foo/bar/baz"));
		
		// Exand to home dir
		assertTrue(Utils.tildeToHomeDir("~/foo/bar/baz").startsWith(File.separator));

		// No change
		assertEquals("/~foo/~bar/baz", Utils.tildeToHomeDir("/~foo/~bar/baz"));
		
		// No change
		assertEquals("~foo", Utils.tildeToHomeDir("~foo"));
		
		//Expanded to /Users/berald01/
		assertEquals(System.getProperty("user.home") + File.separator, Utils.tildeToHomeDir("~/"));
		
		// This will fail most likelly fail on systems other than *nix:
		assertTrue(Utils.tildeToHomeDir("~/foo/bar/baz").startsWith("/Users/") || Utils.tildeToHomeDir("~/foo/bar/baz").startsWith("/home/"));
	}
	
}