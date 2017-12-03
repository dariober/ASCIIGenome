package utils;

import static org.junit.Assert.*;

import java.io.File;
import java.io.IOException;
import java.net.MalformedURLException;

import org.junit.Test;

public class MakeBedgraphTest {

	@Test
	public void canConvertBedToBedgraph() throws MalformedURLException, IOException {
		File bdgOut= new File("test.bedGraph.gz");
		File bdgTbi= new File("test.bedGraph.gz.tbi");
		bdgOut.deleteOnExit();
		bdgTbi.deleteOnExit();
		new MakeBedgraph("test_data/refSeq.hg19.short.sort.bed", 5, bdgOut);
		assertTrue(bdgOut.exists());
		assertTrue(bdgTbi.exists());
		bdgOut.delete(); bdgTbi.delete(); 
		
	}

	@Test
	public void canConvertUrlBedToBedgraph() throws MalformedURLException, IOException {
		File bdgOut= new File("test.bedGraph.gz");
		File bdgTbi= new File("test.bedGraph.gz.tbi");
		bdgOut.deleteOnExit();
		bdgTbi.deleteOnExit();
		new MakeBedgraph("https://raw.githubusercontent.com/dariober/ASCIIGenome/master/test_data/refSeq.hg19.short.sort.bed", 5, bdgOut);
		assertTrue(bdgOut.exists());
		assertTrue(bdgTbi.exists());
		bdgOut.delete(); bdgTbi.delete();
	}
}
