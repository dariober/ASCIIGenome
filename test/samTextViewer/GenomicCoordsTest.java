package samTextViewer;

import static org.junit.Assert.*;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import tracks.TrackSet;
import tracks.TrackWiggles;

import org.junit.Test;

import exceptions.InvalidGenomicCoordsException;

public class GenomicCoordsTest {
	
	static SamReaderFactory srf=SamReaderFactory.make();
	static SamReader samReader= srf.open(new File("test_data/ds051.short.bam"));
	public static SAMSequenceDictionary samSeqDict= samReader.getFileHeader().getSequenceDictionary();
	
	public static String fastaFile= "test_data/chr7.fa";
	
	@Test
	public void canPrintChromMap() throws InvalidGenomicCoordsException, IOException{
			
		GenomicCoords gc= new GenomicCoords("chr7", 1, 1, samSeqDict, 10, null);
		
		String chromMap= gc.getChromIdeogram(10);		
		assertEquals("*---------", chromMap);

		gc= new GenomicCoords("chr7", 1, 1, samSeqDict, 117, null);
		chromMap= gc.getChromIdeogram(10);
		assertEquals(117, chromMap.length());
		
		gc= new GenomicCoords("chr7", 1, 1000000000, samSeqDict, 10, null);
		chromMap= gc.getChromIdeogram(10);		
		assertEquals("**********", chromMap);
		
		gc= new GenomicCoords("chr7", 200000000, 200000000, samSeqDict, 10, null);
		chromMap= gc.getChromIdeogram(10);		
		assertEquals("1--------*", chromMap);
		
		
		gc= new GenomicCoords("chr7", 20000000, 55000000, samSeqDict, 16, null);
		chromMap= gc.getChromIdeogram(10);
		assertEquals("1-****----110M--", chromMap);
	}
	
	@Test
	public void printRefSeq() throws InvalidGenomicCoordsException, IOException{
		GenomicCoords gc= new GenomicCoords("chr7", 5540580, 5540590, null, 100, "test_data/chr7.fa");
		//GenomicCoords gc= new GenomicCoords("chr7", 5540580, 5540590, samSeqDict, 100, null);
		System.out.println(gc.printableRefSeq(true));
	}
	
	@Test
	public void canTestForEqualCoords() throws InvalidGenomicCoordsException, IOException{
		GenomicCoords gc= new GenomicCoords("chr1", 1, 10, null, 100, null);
		GenomicCoords other= new GenomicCoords("chr1", 1, 10, null, 1000, null);
		assertTrue(gc.equalCoords(other));
		
		GenomicCoords other2= new GenomicCoords("chr2", 1, 10, null, 1000, null);
		assertTrue(!gc.equalCoords(other2));
		other2= new GenomicCoords("chr1", 2, 10, null, 1000, null);
		assertTrue(!gc.equalCoords(other2));
		other2= new GenomicCoords("chr1", 1, 100, null, 1000, null);
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
	public void canInitializeSamSeqDictFromGenomeFile() throws IOException{
	
		List<String> insam= new ArrayList<String>();
		// From resource:
		assertEquals(93, GenomicCoords.getSamSeqDictFromAnyFile(insam, null, "hg19").size());
		// From bam header:
		assertEquals(25, GenomicCoords.getSamSeqDictFromAnyFile(insam, null, "test_data/ds051.short.bam").size());
	}
	
	@Test
	public void canGetSamSeqDict() throws IOException{
		List<String> insam= new ArrayList<String>();
		insam.add("test_data/ds051.short.bam.bai"); // This will not produce anything
		insam.add("test_data/ds051.short.bam");
		SAMSequenceDictionary ssd = GenomicCoords.getSamSeqDictFromAnyFile(insam, null, null);
		assertEquals(25, ssd.size());
		
		// From indexed fasta
		insam= new ArrayList<String>();
		insam.add("test_data/ds051.short.bam.bai"); // This will not produce anything
		ssd = GenomicCoords.getSamSeqDictFromAnyFile(null, fastaFile, null);
		assertEquals(1, ssd.size());		
	}
	
	@Test
	public void canPrintRefSeq() throws InvalidGenomicCoordsException, IOException{
		GenomicCoords gc= new GenomicCoords("chr7", 5566770, 5566790, samSeqDict, 100, fastaFile);
		assertEquals("CACTTGGCCTCATTTTTAAGG\n", gc.printableRefSeq(true));
		// with format
		gc= new GenomicCoords("chr7", 5566770, 5566772, samSeqDict, 100, fastaFile);
		assertEquals("[107;31mC[0m[107;34mA[0m[107;31mC[0m\n", gc.printableRefSeq(false));
	}
	
	@Test
	public void canConstructGenomicCoords() throws InvalidGenomicCoordsException, IOException{
		
		GenomicCoords gc= new GenomicCoords("chr7", 1, 100, samSeqDict, 1000, fastaFile);
		assertEquals("chr7", gc.getChrom());
		assertEquals(1, (int)gc.getFrom());
		assertEquals(100, (int)gc.getTo());

		gc= new GenomicCoords("chr7:1-100", samSeqDict, 1000, fastaFile);
		assertEquals("chr7", gc.getChrom());
		assertEquals(1, (int)gc.getFrom());
		assertEquals(100, (int)gc.getTo());

		gc= new GenomicCoords("chr7", 11, null, samSeqDict, 1000, fastaFile);		
		assertEquals(1010, (int)gc.getTo());

		gc= new GenomicCoords("chr7", null, null, samSeqDict, 1000, fastaFile);		
		assertEquals(1, (int)gc.getFrom());
		assertEquals(1000, (int)gc.getTo());

		gc= new GenomicCoords("chr7", 1000000000, 1000000000, samSeqDict, 1000, fastaFile); // Reset to size of chrom
		assertEquals(159138663, (int)gc.getFrom());

		gc= new GenomicCoords("chr7:100", samSeqDict, 1000, fastaFile);
		assertEquals(100, (int)gc.getFrom());
		assertEquals(100+1000-1, (int)gc.getTo());

		gc= new GenomicCoords("chr7:1,000,000,000", samSeqDict, 1000, null);
		assertEquals(159138663-1000+1, (int)gc.getFrom()); // Reset to chrom size
		assertEquals(159138663, (int)gc.getTo());

		gc= new GenomicCoords("chr7:1,000,000,000", null, 1000, null);
		assertEquals(1000000000, (int)gc.getFrom()); // Fine, no dict to check against.
		assertEquals(1000000000+1000-1, (int)gc.getTo());

		
	}
	
	@Test
	public void canGetRefSeq() throws InvalidGenomicCoordsException, IOException{
		GenomicCoords gc= new GenomicCoords("chr7", 5566770, 5566790, samSeqDict, 1000, fastaFile);
		assertEquals("CACTTGGCCTCATTTTTAAGG", new String(gc.getRefSeq()));
		gc= new GenomicCoords("chr7", 5566770, 5566790, samSeqDict, 20, fastaFile);
		// System.out.println(gc.getBpPerScreenColumn());
		assertEquals(null, gc.getRefSeq());
	}
	
	@Test(expected = InvalidGenomicCoordsException.class)
	public void canThrowNullChrom() throws InvalidGenomicCoordsException, IOException {
		new GenomicCoords(null, 1, 100, samSeqDict, 100, fastaFile);
	}
	
	@Test(expected = InvalidGenomicCoordsException.class)
	public void canThrowInvalidCoords() throws InvalidGenomicCoordsException, IOException {
		new GenomicCoords("chr7", -1, 100, samSeqDict, 100, fastaFile);
	}
	
	@Test(expected = InvalidGenomicCoordsException.class)
	public void canThrowChromNotInDict() throws InvalidGenomicCoordsException, IOException {
		new GenomicCoords("nonsense", 1, 100, samSeqDict, 100, fastaFile);
	}
	
	@Test
	public void canZoom() throws InvalidGenomicCoordsException, IOException{
		GenomicCoords gc= new GenomicCoords("chr1:101-105", samSeqDict, 100, null); 
		gc.zoomOut();
		assertEquals(99, (int)gc.getFrom()); // exp 95,111 if zoom fact is x2
		assertEquals(107, (int)gc.getTo()); // 
		for(int i= 0; i < 40; i++){
			gc.zoomOut();
		}
		assertEquals(1, (int)gc.getFrom());
		assertEquals(samSeqDict.getSequence("chr1").getSequenceLength(), (int)gc.getTo()); // Doesn't extend beyond chrom
		
		gc= new GenomicCoords("chr1:101-1000", samSeqDict, 100, null);
		gc.zoomIn();
		assertEquals(326, (int)gc.getFrom());
		assertEquals(776, (int)gc.getTo());
		
		// Zoom-in in small interval has no effect
		gc= new GenomicCoords("chr1:101-200", samSeqDict, 200, null);
		gc.zoomIn();
		assertEquals(101, (int)gc.getFrom());
		assertEquals(200, (int)gc.getTo());
		
		gc= new GenomicCoords("chr1:1-200", samSeqDict, 200, null);
		gc.zoomIn();
		assertEquals(1, (int)gc.getFrom());
		assertEquals(200, (int)gc.getTo());

		gc= new GenomicCoords("chrM:16561-16571", samSeqDict, 200, null); // End of chrom
		gc.zoomIn();
		assertEquals(16561, (int)gc.getFrom());
		assertEquals(16571, (int)gc.getTo());
	}
	
	@Test
	public void canPrepareRuler() throws InvalidGenomicCoordsException, IOException{
		GenomicCoords gc= new GenomicCoords("chr1:101-110", samSeqDict, 7, null);
		assertEquals(7, gc.getMapping().size());
		assertEquals(101.0, gc.getMapping().get(0), 0.01);
		assertEquals(102.5, gc.getMapping().get(1), 0.01);
		assertEquals(102.5, gc.getMapping().get(1), 0.01);
		assertEquals(1.42857142, gc.getBpPerScreenColumn(), 0.01);
	}
	
	@Test
	public void canPrintRuler() throws InvalidGenomicCoordsException, IOException{
		
		GenomicCoords gc= new GenomicCoords("chr1:101-200", samSeqDict, 50, null);
		assertEquals(50, gc.printableRuler(10).length());
	
	}
	
	@Test
	public void canGetGCProfileInRegion() throws InvalidGenomicCoordsException, IOException{
				
		GenomicCoords gc= new GenomicCoords("chr7", 1000000, 1000500, samSeqDict, 50, null);
		assertEquals(null, gc.getGCProfile()); // null fasta
		
		gc= new GenomicCoords("chr7", 1, 100, samSeqDict, 50, fastaFile);
		System.out.println("START");
		System.out.println(gc.getGCProfile().getFileTag());
		System.out.println(gc.getGCProfile().printToScreen());
		
		gc= new GenomicCoords("seq", 1, 120, null, 50, "test_data/seq_cg.fa");
		TrackWiggles gcCnt= gc.getGCProfile();
		gcCnt.setyMaxLines(2);
		String exp= "                                     .::::::::::::\n" +
                    "::::::::::::::::.....____________.::::::::::::::::";
		assertEquals(exp, gcCnt.printToScreen());
		System.out.println(gcCnt.getTitle());
		System.out.println(gcCnt.printToScreen());
	}

	@Test
	public void canCenterAndExtendGenomicCoords() throws InvalidGenomicCoordsException, IOException{
		GenomicCoords gc= new GenomicCoords("chr1", 10000, 20000, null, 100, null);
		int size= 100;
		double slop= 5.0;
		gc.centerAndExtendGenomicCoords(gc, size, slop);
		assertEquals(9550, (int)gc.getFrom());
		assertEquals(10550, (int)gc.getTo());
		
		gc= new GenomicCoords("chr1", 10000, 20000, null, 100, null);
		size= 3;
		gc.centerAndExtendGenomicCoords(gc, size, 3.3); // Extended size smaller then windowSize
		assertTrue(gc.getUserWindowSize() <= gc.getGenomicWindowSize());
		
		gc= new GenomicCoords("chr1", 1, 300, null, 100, null);
		size= 100;
		gc.centerAndExtendGenomicCoords(gc, size, 5.0); 
		assertEquals(1, (int)gc.getFrom());

	}

}
