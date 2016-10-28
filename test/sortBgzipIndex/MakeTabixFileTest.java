package sortBgzipIndex;

import static org.junit.Assert.*;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;

import org.junit.Test;

import com.google.common.io.Files;

import exceptions.InvalidGenomicCoordsException;
import exceptions.InvalidRecordException;
import htsjdk.tribble.index.tabix.TabixFormat;
import htsjdk.tribble.readers.TabixReader;
import htsjdk.tribble.readers.TabixReader.Iterator;

public class MakeTabixFileTest {

	@Test 
	public void canHandleEmptyFile() throws ClassNotFoundException, IOException, InvalidRecordException, SQLException, InvalidGenomicCoordsException{
		String infile= "test_data/empty.bed";
		File outfile= new File("test_data/empty.tmp.bed.gz");
		outfile.deleteOnExit();
		File expectedTbi= new File(outfile.getAbsolutePath() + ".tbi"); 
		expectedTbi.deleteOnExit();
		
		new MakeTabixIndex(infile, outfile, TabixFormat.BED);
		
	}
	
	@Test
	public void canCompressAndIndexSortedFile() throws IOException, InvalidRecordException, ClassNotFoundException, SQLException {
		
		String infile= "test_data/overlapped.bed";
		File outfile= new File("test_data/tmp.bed.gz");
		outfile.deleteOnExit();
		
		File expectedTbi= new File(outfile.getAbsolutePath() + ".tbi"); 
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
		
		File expectedTbi= new File(outfile.getAbsolutePath() + ".tbi"); 
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
		
		File expectedTbi= new File(outfile.getAbsolutePath() + ".tbi"); 
		expectedTbi.deleteOnExit();
		
		new MakeTabixIndex(infile, outfile, TabixFormat.BED);
		
		assertTrue(outfile.exists());
		assertTrue(outfile.length() > 200000);
		assertTrue(expectedTbi.exists());
		assertTrue(expectedTbi.length() > 100000);
		
	}
	
	//@Test
	public void canCompressAndIndexSortedURL() throws IOException, InvalidRecordException, ClassNotFoundException, SQLException {
		
		String infile= "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeDacMapabilityConsensusExcludable.bed.gz";
		File outfile= new File("test_data/tmp4.bed.gz");
		outfile.deleteOnExit();
		
		File expectedTbi= new File(outfile.getAbsolutePath() + ".tbi"); 
		expectedTbi.deleteOnExit();
		
		new MakeTabixIndex(infile, outfile, TabixFormat.BED);
		
		assertTrue(outfile.exists());
		assertTrue(outfile.length() > 1000);
		assertTrue(expectedTbi.exists());
		assertTrue(expectedTbi.length() > 1000);
		
	}
	
	// @Test
	public void canCompressAndIndexUnsortedURL() throws IOException, InvalidRecordException, ClassNotFoundException, SQLException {
		
		String infile= "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeSydhTfbs/wgEncodeSydhTfbsGm12878P300bStdPk.narrowPeak.gz";
		File outfile= new File("test_data/tmp5.bed.gz");
		outfile.deleteOnExit();
		
		File expectedTbi= new File(outfile.getAbsolutePath() + ".tbi"); 
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
		//outfile.deleteOnExit();
		
		File expectedTbi= new File(outfile.getAbsolutePath() + ".tbi"); 
		//expectedTbi.deleteOnExit();

		new MakeTabixIndex(infile, outfile, TabixFormat.VCF);
		
		assertTrue(outfile.exists());
		assertTrue(outfile.length() > 1000);
		assertTrue(expectedTbi.exists());
		assertTrue(expectedTbi.length() > 1000);

		TabixReader tbx = new TabixReader(outfile.getAbsolutePath());
		Iterator x = tbx.query("1", 20000000, 30000000);
		assertTrue(x.next().startsWith("1"));

	}
	
	@Test
	public void blockCompressFileInPlace() throws IOException, ClassNotFoundException, InvalidRecordException, SQLException{
		// Test we can compress and index a file and overwrite the original file. 
		
		File testFile= new File("deleteme.bed.gz");
		testFile.deleteOnExit();
		Files.copy(new File("test_data/refSeq.hg19.bed.gz"), testFile);
		
		File expectedTbi= new File(testFile.getAbsolutePath() + ".tbi"); 
		expectedTbi.deleteOnExit();
		
		new MakeTabixIndex(testFile.getAbsolutePath(), testFile, TabixFormat.BED);
		
		assertTrue(testFile.exists());
		assertTrue(testFile.length() > 200000);
		assertTrue(expectedTbi.exists());
		assertTrue(expectedTbi.length() > 100000);

	}
	
}
