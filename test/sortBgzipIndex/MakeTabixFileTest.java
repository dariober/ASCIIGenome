package sortBgzipIndex;

import static org.junit.Assert.*;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;

import org.junit.Test;

import exceptions.InvalidRecordException;
import htsjdk.tribble.index.tabix.TabixFormat;

public class MakeTabixFileTest {

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
	
	@Test
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
	
	@Test
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
	
}
