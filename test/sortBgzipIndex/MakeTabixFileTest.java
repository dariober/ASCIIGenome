package sortBgzipIndex;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;

import org.junit.Test;

import com.google.common.io.Files;

import exceptions.InvalidGenomicCoordsException;
import exceptions.InvalidRecordException;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.tribble.index.tabix.TabixFormat;
import htsjdk.tribble.readers.TabixReader;
import htsjdk.tribble.readers.TabixReader.Iterator;
import htsjdk.tribble.util.TabixUtils;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

public class MakeTabixFileTest {

	public void vcfTester(String inVcf){
		// Credit: https://www.biostars.org/p/262943/
		VCFFileReader r=new VCFFileReader(new File(inVcf));
		CloseableIterator<VariantContext> t=r.iterator();
		while(t.hasNext()){
		    t.next();
		}
		t.close();
		r.close();
	} 
	
	@Test
	public void canCompressAndIndexHeaderlessVCF() throws ClassNotFoundException, IOException, InvalidRecordException, SQLException{

		String infile= "test_data/noheader.vcf";
		File outfile= new File("test_data/noheader.vcf.gz");
		outfile.deleteOnExit();
		
		File expectedTbi= new File(outfile.getAbsolutePath() + TabixUtils.STANDARD_INDEX_EXTENSION); 
		expectedTbi.deleteOnExit();

		new MakeTabixIndex(infile, outfile, TabixFormat.VCF);
		
		assertTrue(outfile.exists());
		assertTrue(outfile.length() > 200);
		assertTrue(expectedTbi.exists());
		assertTrue(expectedTbi.length() > 100);

		TabixReader tbx = new TabixReader(outfile.getAbsolutePath());
		Iterator x = tbx.query("1", 1, 10000000);
		assertTrue(x.next().startsWith("1"));

	}	
	
	@Test
	public void testRealFileSizeVCF() throws ClassNotFoundException, IOException, InvalidRecordException, SQLException{
		
		// See test_data/README.md for this file. This is fairly large and we want to check it is processed
		// in a reasonable amount of time.
		String infile= "test_data/ALL.wex.union_illumina_wcmc_bcm_bc_bi.20110521.snps.exome.sites.vcf";
		
		File outfile= new File("deleteme.vcf.gz");
		outfile.deleteOnExit();
		File expectedTbi= new File(outfile.getAbsolutePath() + TabixUtils.STANDARD_INDEX_EXTENSION); 
		expectedTbi.deleteOnExit();
		
		long t0= System.currentTimeMillis();
		new MakeTabixIndex(infile, outfile, TabixFormat.VCF);
		long t1= System.currentTimeMillis();
		
		assertTrue(outfile.exists());
		assertTrue(outfile.length() > 1000);
		assertTrue((t1 - t0) < 20000); // Should be << than 20 sec, ~2 sec

		// Check you can read ok
		this.vcfTester(outfile.getAbsolutePath());
	}
	
	@Test
	public void canCompressAndIndexVCF_CEU() throws ClassNotFoundException, IOException, InvalidRecordException, SQLException{

		String infile= "test_data/CEU.exon.2010_06.genotypes.vcf";
		
		File outfile= new File("deleteme.vcf.gz");
		outfile.deleteOnExit();
		File expectedTbi= new File(outfile.getAbsolutePath() + TabixUtils.STANDARD_INDEX_EXTENSION); 
		expectedTbi.deleteOnExit();
		
		new MakeTabixIndex(infile, outfile, TabixFormat.VCF);

		assertTrue(outfile.exists());
		assertTrue(outfile.length() > 1000);
		
		// Check you can read ok
		this.vcfTester(outfile.getAbsolutePath());
	}
	
	// @Test
	public void handlingInvalidFile() throws ClassNotFoundException, IOException, InvalidRecordException, SQLException{
		
		String infile= "test_data/chr7.fa";
		File outfile= new File("deleteme.gtf.gz");
		outfile.deleteOnExit();
		File expectedTbi= new File(outfile.getAbsolutePath() + TabixUtils.STANDARD_INDEX_EXTENSION); 
		expectedTbi.deleteOnExit();
		
		new MakeTabixIndex(infile, outfile, TabixFormat.BED);
	}
	
	@Test 
	public void canHandleEmptyFile() throws ClassNotFoundException, IOException, InvalidRecordException, SQLException, InvalidGenomicCoordsException{
		String infile= "test_data/empty.bed";
		File outfile= new File("test_data/empty.tmp.bed.gz");
		outfile.deleteOnExit();
		File expectedTbi= new File(outfile.getAbsolutePath() + TabixUtils.STANDARD_INDEX_EXTENSION); 
		expectedTbi.deleteOnExit();
		
		new MakeTabixIndex(infile, outfile, TabixFormat.BED);
	
		assertTrue(outfile.exists());
		assertTrue(expectedTbi.exists());
		
		assertEquals(28, outfile.length()); // Checked with `bgzip empty.bed && ls -l empty.bed.gz`
		assertTrue(expectedTbi.length() > 0);

	}
	
	@Test
	public void canCompressAndIndexSortedFile() throws IOException, InvalidRecordException, ClassNotFoundException, SQLException {
		
		String infile= "test_data/overlapped.bed";
		File outfile= new File("test_data/tmp.bed.gz");
		outfile.deleteOnExit();
		
		File expectedTbi= new File(outfile.getAbsolutePath() + TabixUtils.STANDARD_INDEX_EXTENSION); 
		expectedTbi.deleteOnExit();
		
		new MakeTabixIndex(infile, outfile, TabixFormat.BED);
		
		assertTrue(outfile.exists());
		assertTrue(outfile.length() > 80);
		assertTrue(expectedTbi.exists());
		assertTrue(expectedTbi.length() > 80);
		
		TabixReader tbx = new TabixReader(outfile.getAbsolutePath());
		Iterator x = tbx.query("chr1", 1, 1000000);
		assertTrue(x.next().startsWith("chr1"));
		
	}

	@Test
	public void canCompressAndIndexSortedGzipFile() throws IOException, InvalidRecordException, ClassNotFoundException, SQLException {
		
		String infile= "test_data/hg19_genes.gtf.gz";
		File outfile= new File("test_data/tmp2.bed.gz");
		outfile.deleteOnExit();
		
		File expectedTbi= new File(outfile.getAbsolutePath() + TabixUtils.STANDARD_INDEX_EXTENSION); 
		expectedTbi.deleteOnExit();
		
		new MakeTabixIndex(infile, outfile, TabixFormat.GFF);

		assertTrue(outfile.exists());
		assertTrue(outfile.length() > 7000000);
		assertTrue(expectedTbi.exists());
		assertTrue(expectedTbi.length() > 500000);
		
		TabixReader tbx = new TabixReader(outfile.getAbsolutePath());
		Iterator x = tbx.query("chr1", 1, 1000000);
		assertTrue(x.next().startsWith("chr1"));
		
	}
	
	
	@Test
	public void canCompressAndIndexUnsortedFile() throws IOException, InvalidRecordException, ClassNotFoundException, SQLException {
		
		String infile= "test_data/refSeq.hg19.bed.gz";
		File outfile= new File("test_data/tmp3.bed.gz");
		outfile.deleteOnExit();
		
		File expectedTbi= new File(outfile.getAbsolutePath() + TabixUtils.STANDARD_INDEX_EXTENSION); 
		expectedTbi.deleteOnExit();
		
		new MakeTabixIndex(infile, outfile, TabixFormat.BED);
		
		assertTrue(outfile.exists());
		assertTrue(outfile.length() > 200000);
		assertTrue(expectedTbi.exists());
		assertTrue(expectedTbi.length() > 100000);
		
	}
	
	@Test
	public void canCompressAndIndexSortedURL() throws IOException, InvalidRecordException, ClassNotFoundException, SQLException {
		
		String infile= "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeDacMapabilityConsensusExcludable.bed.gz";
		File outfile= new File("test_data/tmp4.bed.gz");
		outfile.deleteOnExit();
		
		File expectedTbi= new File(outfile.getAbsolutePath() + TabixUtils.STANDARD_INDEX_EXTENSION); 
		expectedTbi.deleteOnExit();
		
		new MakeTabixIndex(infile, outfile, TabixFormat.BED);
		
		assertTrue(outfile.exists());
		assertTrue(outfile.length() > 1000);
		assertTrue(expectedTbi.exists());
		assertTrue(expectedTbi.length() > 1000);
		
	}
	
	@Test
	public void canCompressAndIndexUnsortedURL() throws IOException, InvalidRecordException, ClassNotFoundException, SQLException {
		
		String infile= "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeSydhTfbs/wgEncodeSydhTfbsGm12878P300bStdPk.narrowPeak.gz";
		File outfile= new File("test_data/tmp5.bed.gz");
		outfile.deleteOnExit();
		
		File expectedTbi= new File(outfile.getAbsolutePath() + TabixUtils.STANDARD_INDEX_EXTENSION); 
		expectedTbi.deleteOnExit();
		
		new MakeTabixIndex(infile, outfile, TabixFormat.BED);
		
		assertTrue(outfile.exists());
		assertTrue(outfile.length() > 1000);
		assertTrue(expectedTbi.exists());
		assertTrue(expectedTbi.length() > 1000);
		
	}
	
	@Test
	public void canCompressAndIndexVCF() throws ClassNotFoundException, IOException, InvalidRecordException, SQLException{

		String infile= "test_data/CHD.exon.2010_03.sites.unsorted.vcf";
		File outfile= new File("test_data/tmp6.bed.gz");
		outfile.deleteOnExit();
		
		File expectedTbi= new File(outfile.getAbsolutePath() + TabixUtils.STANDARD_INDEX_EXTENSION); 
		expectedTbi.deleteOnExit();

		new MakeTabixIndex(infile, outfile, TabixFormat.VCF);
		
		assertTrue(outfile.exists());
		assertTrue(outfile.length() > 1000);
		assertTrue(expectedTbi.exists());
		assertTrue(expectedTbi.length() > 1000);

		TabixReader tbx = new TabixReader(outfile.getAbsolutePath());
		Iterator x = tbx.query("1", 20000000, 30000000);
		assertTrue(x.next().startsWith("1"));

		// Check you can read ok
		this.vcfTester(outfile.getAbsolutePath());
	}
	
	@Test
	public void blockCompressFileInPlace() throws IOException, ClassNotFoundException, InvalidRecordException, SQLException{
		// Test we can compress and index a file and overwrite the original file. 
		
		File testFile= new File("deleteme.bed.gz");
		testFile.deleteOnExit();
		Files.copy(new File("test_data/refSeq.hg19.bed.gz"), testFile);
		
		File expectedTbi= new File(testFile.getAbsolutePath() + TabixUtils.STANDARD_INDEX_EXTENSION); 
		expectedTbi.deleteOnExit();
		
		new MakeTabixIndex(testFile.getAbsolutePath(), testFile, TabixFormat.BED);
		
		assertTrue(testFile.exists());
		assertTrue(testFile.length() > 200000);
		assertTrue(expectedTbi.exists());
		assertTrue(expectedTbi.length() > 100000);

	}
	
	@Test
	public void canProcessBedgraphWithTrackLine() throws ClassNotFoundException, IOException, InvalidRecordException, SQLException{

		String testFile= "test_data/test2.bedGraph";

		String outfile= testFile + ".tmp.bedGraph.gz";
		File expectedTbi= new File(outfile + TabixUtils.STANDARD_INDEX_EXTENSION); 
		new File(outfile).deleteOnExit();
		expectedTbi.deleteOnExit();
		
		new MakeTabixIndex(testFile, new File(outfile), TabixFormat.BED);
		assertTrue(expectedTbi.exists());
		assertTrue(expectedTbi.length() > 50);
	}
	
}
