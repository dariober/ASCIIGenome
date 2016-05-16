package tracks;

import static org.junit.Assert.*;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.junit.Test;

import exceptions.InvalidGenomicCoordsException;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.filter.SamRecordFilter;
import samTextViewer.GenomicCoords;

public class TrackReadsTest {

	public static List<SamRecordFilter> filters= new ArrayList<SamRecordFilter>();

	static SamReaderFactory srf=SamReaderFactory.make();
	static SamReader samReader= srf.open(new File("test_data/ds051.actb.bam"));
	public static SAMSequenceDictionary samSeqDict= samReader.getFileHeader().getSequenceDictionary();

	
	public static String fastaFile= "test_data/chr7.fa";
	
	@Test
	public void testOneLineStack() throws InvalidGenomicCoordsException, IOException{
		String bam= "test_data/adjacent.bam";
		int yMaxLines= 10;
		boolean bs= false;
		boolean noFormat= true;
		int windowSize= 100;
		int from= 1;
		int to= 100;
		GenomicCoords gc= new GenomicCoords("chr7", from, to, samSeqDict, windowSize, null);
			
		TrackReads tr= new TrackReads(bam, gc, filters, 1000);
		tr.setBs(bs);
		tr.setNoFormat(noFormat);
		String exp= 
		"AAAAAAAAAA           GGGGGGGGGG TTTTTTTTTT\n"+
		"          CCCCCCCCCC";
		tr.setyMaxLines(yMaxLines);
		assertEquals(exp, tr.printToScreen());
		
		System.out.println(tr.printToScreen());
		
		windowSize= 99;
		gc= new GenomicCoords("chr7", from, to, samSeqDict, windowSize, null);
		tr= new TrackReads(bam, gc, filters, 1000);
		tr.setBs(bs);
		tr.setNoFormat(noFormat);
		tr.setyMaxLines(yMaxLines);
		assertEquals(">>>>>>>>>>>>>>>>>>>> >>>>>>>>>> >>>>>>>>>>", tr.printToScreen());
	}
	String exp= "AAAAAAAAAA           GGGGGGGGGG TTTTTTTTTT\n"+
	"          CCCCCCCCCC";
	@Test
	public void test() throws InvalidGenomicCoordsException, IOException {
		String bam= "test_data/ds051.actb.bam";
		int maxReadsStack= 1000;
		int yMaxLines= 50;
		boolean bs= false;
		boolean noFormat= true;
		int windowSize= 200;
		int from= 5566770;
		int to= from + 100;
		GenomicCoords gc= new GenomicCoords("chr7", from, to, samSeqDict, windowSize, fastaFile);
			
		TrackReads tr= new TrackReads(bam, gc, filters, maxReadsStack);
		tr.setBs(bs);
		tr.setNoFormat(noFormat);
		tr.setyMaxLines(yMaxLines);
		System.out.println(tr.printToScreen());
	}

	@Test
	public void testNoReadsInRegion() throws InvalidGenomicCoordsException, IOException {
		String bam= "test_data/ds051.actb.bam";
		int maxReadsStack= 1000;
		int yMaxLines= 50;
		boolean bs= false;
		boolean noFormat= true;
		int windowSize= 200;
		int from= 1;
		int to= from + 100;
		GenomicCoords gc= new GenomicCoords("chr7", from, to, samSeqDict, windowSize, fastaFile);
			
		TrackReads tr= new TrackReads(bam, gc, filters, maxReadsStack);
		tr.setyMaxLines(yMaxLines);
		tr.setBs(bs);
		tr.setNoFormat(noFormat);
		System.out.println(tr.printToScreen());
	}
	
	@Test
	public void canResetToZeroLargeWindow() throws IOException, InvalidGenomicCoordsException {
		// If the genomic window is too large do not process the bam file and return zero height track.
		GenomicCoords gc= new GenomicCoords("chr7", 1, 100000000, samSeqDict, 100, fastaFile);
		TrackReads tr= new TrackReads("test_data/ds051.actb.bam", gc, filters, 1000);
		assertEquals("", tr.printToScreen());
	}
}
