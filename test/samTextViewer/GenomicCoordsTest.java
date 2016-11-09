package samTextViewer;

import static org.junit.Assert.*;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.List;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import tracks.TrackWiggles;

import org.junit.Test;

import exceptions.InvalidGenomicCoordsException;
import exceptions.InvalidRecordException;

public class GenomicCoordsTest {
	
	static SamReaderFactory srf=SamReaderFactory.make();
	static SamReader samReader= srf.open(new File("test_data/ds051.short.bam"));
	public static SAMSequenceDictionary samSeqDict= samReader.getFileHeader().getSequenceDictionary();
	
	public static String fastaFile= "test_data/chr7.fa";
	
	@Test
	public void canGetStringRegion() throws InvalidGenomicCoordsException, IOException{
		assertEquals("chr7:1-100", new GenomicCoords("chr7:1-100", null, null).toStringRegion());
		
		GenomicCoords gc= new GenomicCoords("chr7:1", null, null);
		assertTrue(gc.toStringRegion().matches("^chr7:1-\\d+$")); // String region is something like chr7:1-79
		
		gc= new GenomicCoords("chr7", null, null);
		assertTrue(gc.toStringRegion().matches("^chr7:1-\\d+$")); // chr7:1-79
	}
	
	@Test
	public void canPrintChromMap() throws InvalidGenomicCoordsException, IOException{
			
		GenomicCoords gc= new GenomicCoords("chr7:1-100", samSeqDict, null);
		
		String chromMap= gc.getChromIdeogram(10, true);		
		assertTrue(chromMap.startsWith("*---------"));

		gc= new GenomicCoords("chr7:1-100", samSeqDict, null);
		chromMap= gc.getChromIdeogram(10, true);
		assertEquals(79, chromMap.length()); // 79 Should be the size of eclipse's terminal
		
		gc= new GenomicCoords("chr7:1-1000000000", samSeqDict, null);
		chromMap= gc.getChromIdeogram(10, true);		
		assertTrue(chromMap.startsWith("**********"));
		
		gc= new GenomicCoords("chr7:200000000-300000000", samSeqDict, null);
		chromMap= gc.getChromIdeogram(10, true);		
		assertTrue(chromMap.startsWith("1--------"));
		
	}
	
	@Test
	public void printRefSeq() throws InvalidGenomicCoordsException, IOException{
		GenomicCoords gc= new GenomicCoords("chr7:5540580-5540590", null, "test_data/chr7.fa");
		assertEquals("ggccggctggg\n", gc.printableRefSeq(true));
	}
	
	@Test
	public void canTestForEqualCoords() throws InvalidGenomicCoordsException, IOException{
		GenomicCoords gc= new GenomicCoords("chr1:1-10", null, null);
		GenomicCoords other= new GenomicCoords("chr1:1-10", null, null);
		assertTrue(gc.equalCoords(other));
		
		GenomicCoords other2= new GenomicCoords("chr2:1-10", null, null);
		assertTrue(!gc.equalCoords(other2));
		other2= new GenomicCoords("chr1:2-10", null, null);
		assertTrue(!gc.equalCoords(other2));
		other2= new GenomicCoords("chr1:1-100", null, null);
		assertTrue(!gc.equalCoords(other2));
		
		GenomicCoords gc2= (GenomicCoords) gc.clone();
		System.out.println(gc);
		System.out.println(gc2);
		gc.zoomOut();
		System.out.println(gc);
		System.out.println(gc2);
		// assertTrue(gc2.equalCoords(gc));		
	}
	
	@Test
	public void canInitializeSamSeqDictFromGenomeFile() throws IOException, InvalidGenomicCoordsException{
	
		List<String> insam= new ArrayList<String>();
		insam.add("hg19");
		// From resource:
		GenomicCoords gc= new GenomicCoords("chr1", null, null);
		gc.setGenome(insam);
		assertEquals(93, gc.getSamSeqDict().size());
		// From bam header:
		insam.set(0, "test_data/ds051.short.bam");
		gc.setGenome(insam);
		assertEquals(25, gc.getSamSeqDict().size());
		
		// Check we get the full path to source file.
		assertTrue(gc.getSamSeqDictSource().length() > "test_data/ds051.short.bam".length());
		
	}
	
//	@Test
//	public void canGetSamSeqDict() throws IOException{
//		List<String> insam= new ArrayList<String>();
//		insam.add("test_data/ds051.short.bam.bai"); // This will not produce anything
//		insam.add("test_data/ds051.short.bam");
//		SAMSequenceDictionary ssd = GenomicCoords.getSamSeqDictFromAnyFile(insam, null, null);
//		assertEquals(25, ssd.size());
//		
//		// From indexed fasta
//		insam= new ArrayList<String>();
//		insam.add("test_data/ds051.short.bam.bai"); // This will not produce anything
//		ssd = GenomicCoords.getSamSeqDictFromAnyFile(null, fastaFile, null);
//		assertEquals(1, ssd.size());		
//	}
	
	@Test
	public void canPrintRefSeq() throws InvalidGenomicCoordsException, IOException{
		GenomicCoords gc= new GenomicCoords("chr7:5566770-5566790", samSeqDict, fastaFile);
		assertEquals("CACTTGGCCTCATTTTTAAGG\n", gc.printableRefSeq(true));
		// with format
		gc= new GenomicCoords("chr7:5566770-5566772", samSeqDict, fastaFile);
		assertTrue(gc.printableRefSeq(false).contains("["));
	}
	
	@Test
	public void canConstructGenomicCoords() throws InvalidGenomicCoordsException, IOException{
		
		GenomicCoords gc= new GenomicCoords("chr7:1-100", samSeqDict, fastaFile);
		assertEquals("chr7", gc.getChrom());
		assertEquals(1, (int)gc.getFrom());
		assertEquals(100, (int)gc.getTo());

		gc= new GenomicCoords("chr7:1-100", samSeqDict, fastaFile);
		assertEquals("chr7", gc.getChrom());
		assertEquals(1, (int)gc.getFrom());
		assertEquals(100, (int)gc.getTo());

		gc= new GenomicCoords("chr7:11", samSeqDict, fastaFile);		
		assertEquals(11 + 79 - 1, (int)gc.getTo());

		gc= new GenomicCoords("chr7", samSeqDict, fastaFile);		
		assertEquals(1, (int)gc.getFrom());
		assertEquals(79, (int)gc.getTo());

		gc= new GenomicCoords("chr7:1000000000-1000000000", samSeqDict, fastaFile); // Reset to size of chrom
		assertEquals(159138663, (int)gc.getFrom());

		gc= new GenomicCoords("chr7:100", samSeqDict, fastaFile);
		assertEquals(100, (int)gc.getFrom());
		assertEquals(100 + 79 - 1, (int)gc.getTo());

		gc= new GenomicCoords("chr7:1,000,000,000", samSeqDict, null);
		assertEquals(159138663-79+1, (int)gc.getFrom()); // Reset to chrom size
		assertEquals(159138663, (int)gc.getTo());

		gc= new GenomicCoords("chr7:1,000,000,000", null, null);
		assertEquals(1000000000, (int)gc.getFrom()); // Fine, no dict to check against.
		assertEquals(1000000000 + 79 -1, (int)gc.getTo());

	}
	
	@Test
	public void canGetRefSeq() throws InvalidGenomicCoordsException, IOException{
		GenomicCoords gc= new GenomicCoords("chr7:5566770-5566790", samSeqDict, fastaFile);
		assertEquals("CACTTGGCCTCATTTTTAAGG", new String(gc.getRefSeq()));
		gc= new GenomicCoords("chr7:1-80", samSeqDict, fastaFile);
		assertEquals(null, gc.getRefSeq());
	
		// Return seq even if len(seq) > windowSize
		//assertEquals("CACTTGGCCTCATTTTTAAGG", new String(gc.getRefSeq()));
	}
	
//	@Test(expected = InvalidGenomicCoordsException.class)
//	public void canThrowNullChrom() throws InvalidGenomicCoordsException, IOException {
//		new GenomicCoords(null, samSeqDict, fastaFile);
//	}
	
//	@Test(expected = InvalidGenomicCoordsException.class)
//	public void canThrowInvalidCoords() throws InvalidGenomicCoordsException, IOException {
//		new GenomicCoords("chr7", -1, 100, samSeqDict, fastaFile);
//	}
	
	@Test(expected = InvalidGenomicCoordsException.class)
	public void canThrowChromNotInDict() throws InvalidGenomicCoordsException, IOException {
		new GenomicCoords("nonsense:1-100", samSeqDict, fastaFile);
	}
	
	@Test
	public void canZoom() throws InvalidGenomicCoordsException, IOException{
		GenomicCoords gc= new GenomicCoords("chr1:101-105", samSeqDict, null); 
		gc.zoomOut();
		assertEquals(99, (int)gc.getFrom()); // exp 95,111 if zoom fact is x2
		assertEquals(107, (int)gc.getTo()); // 
		for(int i= 0; i < 40; i++){
			gc.zoomOut();
		}
		assertEquals(1, (int)gc.getFrom());
		assertEquals(samSeqDict.getSequence("chr1").getSequenceLength(), (int)gc.getTo()); // Doesn't extend beyond chrom
		
		gc= new GenomicCoords("chr1:101-1000", samSeqDict, null);
		gc.zoomIn();
		assertEquals(326, (int)gc.getFrom());
		assertEquals(776, (int)gc.getTo());
		
		// Zoom-in in small interval has no effect
		gc= new GenomicCoords("chr1:101-130", samSeqDict, null);
		gc.zoomIn();
		assertEquals(101, (int)gc.getFrom());
		assertEquals(130, (int)gc.getTo());
		
		gc= new GenomicCoords("chr1:1-50", samSeqDict, null);
		gc.zoomIn();
		assertEquals(1, (int)gc.getFrom());
		assertEquals(50, (int)gc.getTo());

		gc= new GenomicCoords("chrM:16561-16571", samSeqDict, null); // End of chrom
		gc.zoomIn();
		assertEquals(16561, (int)gc.getFrom());
		assertEquals(16571, (int)gc.getTo());
	}
	
	@Test
	public void canPrepareRuler() throws InvalidGenomicCoordsException, IOException{
		GenomicCoords gc= new GenomicCoords("chr1:101-110", samSeqDict, null);
		int userWindowSize= gc.getUserWindowSize();
		
		assertEquals(10, gc.getMapping(userWindowSize).size());
		assertEquals(101.0, gc.getMapping(userWindowSize).get(0), 0.01);
		assertEquals(102.0, gc.getMapping(userWindowSize).get(1), 0.01);
		assertEquals(1.0, gc.getBpPerScreenColumn(), 0.01);
	}
	
	@Test
	public void canPrintRuler() throws InvalidGenomicCoordsException, IOException{
		
		GenomicCoords gc= new GenomicCoords("chr1:101-200", samSeqDict, null);
		assertEquals(79, gc.printableRuler(10, true).length());

	}
	
	// @Test // Do not test until gcProfile is sorted
	public void canGetGCProfileInRegion() throws InvalidGenomicCoordsException, IOException, InvalidRecordException, ClassNotFoundException, SQLException{
				
		GenomicCoords gc= new GenomicCoords("chr7:1000000-1000500", samSeqDict, null);
		assertEquals(null, gc.getGCProfile()); // null fasta
		
		gc= new GenomicCoords("chr7:1-100", samSeqDict, fastaFile);
		System.out.println("START");
		System.out.println(gc.getGCProfile().getTrackTag());
		System.out.println(gc.getGCProfile().printToScreen());
		
		gc= new GenomicCoords("seq:1-120", null, "test_data/seq_cg.fa");
		TrackWiggles gcCnt= gc.getGCProfile();
		gcCnt.setyMaxLines(2);
		String exp= "                                  ::::::::::::::::\n" +
                    "::::::::::::::::.________________:::::::::::::::::";
		gcCnt.setNoFormat(true);
		System.out.println(gcCnt.printToScreen());
		assertEquals(exp, gcCnt.printToScreen());
		System.out.println(gcCnt.getTitle());
		System.out.println(gcCnt.printToScreen());
		
	}

	@Test
	public void canCenterAndExtendGenomicCoords() throws InvalidGenomicCoordsException, IOException{
		GenomicCoords gc= new GenomicCoords("chr1:10000-20000", null, null);
		int size= 100;
		double slop= 5.0;
		gc.centerAndExtendGenomicCoords(gc, size, slop);
		assertEquals(9550, (int)gc.getFrom());
		assertEquals(10550, (int)gc.getTo());
		
		gc= new GenomicCoords("chr1:10000-20000", null, null);
		size= 3;
		gc.centerAndExtendGenomicCoords(gc, size, 3.3); // Extended size smaller then windowSize
		assertTrue(gc.getUserWindowSize() <= gc.getGenomicWindowSize());
		
		gc= new GenomicCoords("chr1:1-300", null, null);
		size= 100;
		gc.centerAndExtendGenomicCoords(gc, size, 5.0); 
		assertEquals(1, (int)gc.getFrom());
		
	}

}
