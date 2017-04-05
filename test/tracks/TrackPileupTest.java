package tracks;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.List;
import java.util.Map;

import org.apache.commons.io.FileUtils;
import org.junit.Test;

import com.google.common.base.Splitter;

import exceptions.InvalidGenomicCoordsException;
import exceptions.InvalidRecordException;
import samTextViewer.GenomicCoords;

public class TrackPileupTest {

	@Test
	public void canCollectCoverage() throws ClassNotFoundException, IOException, InvalidGenomicCoordsException, InvalidRecordException, SQLException{
		GenomicCoords gc= new GenomicCoords("chr7:5566736-5566856", null, null);
		TrackPileup tr= new TrackPileup("test_data/ds051.short.bam", gc);
		
		assertEquals(79, tr.getDepth().size());
		assertEquals(1, (int)tr.getDepth().get(5566778)); // Depths checked against mpileup
		assertEquals(5, (int)tr.getDepth().get(5566782));
		assertEquals(18, (int)tr.getDepth().get(5566856));
		
		gc= new GenomicCoords("chr7:5522059-5612125", null, null);
		long t0= System.currentTimeMillis();
		tr= new TrackPileup("test_data/ear045.oxBS.actb.bam", gc);
		
		Map<Integer, Integer> depth = tr.getDepth();
		long t1= System.currentTimeMillis();
		assertTrue(t1-t0 < 10000); // Processing time (in ms) is acceptably small
		assertTrue(t1-t0 > 100); // But not suspiciously small
		System.err.println("Time to parse " + depth.size() + " positions: " + (t1-t0) + " ms");

		// See test_data/README.md for obtaining this test file (samtools mpileup ...)
		String expPileup= FileUtils.readFileToString(new File("test_data/ear045.oxBS.actb.pileup"));
		List<String> expList = Splitter.on("\n").omitEmptyStrings().splitToList(expPileup);
	
		// mpileup and TrackPileup hit the same positions
		assertEquals(expList.size(), depth.size());

		int i= 0;
		for(int obsPos : depth.keySet()){
			int obsDepth=depth.get(obsPos);
			int expPos= Integer.parseInt(Splitter.on("\t").splitToList(expList.get(i)).get(1));
			int expDepth= Integer.parseInt(Splitter.on("\t").splitToList(expList.get(i)).get(3));
			try{
				assertEquals(expPos, obsPos);
				assertTrue(Math.abs((expDepth - obsDepth)) < 10);
			} catch(AssertionError e){
				System.err.println(i);
				System.err.println(expList.get(i));
				System.err.println(obsPos);
				throw e; 
			}
			i++;
		}		
//		
//		System.err.println(expList.get(0));
		
	}

//	@Test
//	public void canAddSAMRecord() {
//	
//		// Prepare a sam header
//		String sqHeader = "@HD\tSO:coordinate\tVN:1.0\n"
//	                        + "@SQ\tSN:chr1\tAS:HG18\tLN:10000000\n";
//
//		// Prepare one read with a 500,000 bases skipped
//		String cigar= "10S" + "5M" + "5D" + "8M" + "5000N" + "3M";
//		String s1 = "read1\t0\tchr1\t100\t255\t" + cigar + "\t*\t0\t0\tnnnnnnnnnnACCTACGTTCAATATTACAGG\t*\n";
//		
//		
//		// Prepare sam input and samReader
//		String exampleSam = sqHeader + s1 + s1;
//		ByteArrayInputStream inputStream = new ByteArrayInputStream(exampleSam.getBytes());
//
//		SamReader samReader= SamReaderFactory.makeDefault().open(SamInputResource.of(inputStream));
//		SAMRecord rec= samReader.iterator().next();
//
//		GenomicCoords gc= new GenomicCoords("chr1:1-10000", null, null);
//
//		TrackPileup pileup= new TrackPileup(gc);
//		pileup.add(rec);
//		// First position is 100 as per sam record
//		assertEquals(100, (int)pileup.getDepth().keySet().iterator().next());
//		
//		// Number of position covered (sum of M and D ops)
//		assertEquals(5+5+8+3, pileup.getDepth().keySet().size());
//		
//		// All pos depth 1:
//		for(int x : pileup.getDepth().values() ){
//			assertEquals(1, x);
//		}
//	}

}
