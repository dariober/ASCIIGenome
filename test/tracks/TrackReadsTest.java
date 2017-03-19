package tracks;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;

import org.junit.Test;

import exceptions.InvalidColourException;
import exceptions.InvalidGenomicCoordsException;
import exceptions.InvalidRecordException;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import samTextViewer.GenomicCoords;

public class TrackReadsTest {

	//public static List<SamRecordFilter> filters= new ArrayList<SamRecordFilter>();

	static SamReaderFactory srf=SamReaderFactory.make();
	static SamReader samReader= srf.open(new File("test_data/ds051.actb.bam"));
	public static SAMSequenceDictionary samSeqDict= samReader.getFileHeader().getSequenceDictionary();

	
	public static String fastaFile= "test_data/chr7.fa";
	
	@Test
	public void canFilterReadsWithAwk() throws InvalidGenomicCoordsException, IOException, ClassNotFoundException, InvalidRecordException, SQLException, InvalidColourException{
		GenomicCoords gc= new GenomicCoords("chr7:5566000-5567000", samSeqDict, null);
		TrackReads tr= new TrackReads("test_data/ds051.short.bam", gc);
		tr.setNoFormat(true);
		tr.setyMaxLines(1000);

		assertEquals(22, tr.printToScreen().split("\n").length); // N. reads stacked in this interval before filtering
		
		tr.setAwk("'$1 ~ \"NCNNNCCC\"'");

		assertEquals(6, tr.printToScreen().split("\n").length);

		System.err.println(Track.awkFunc);
		int i= 0;
		StringBuilder foo= new StringBuilder();
		
		while(i < 100000){
			foo.append("baz");
			i++;
		}
		foo.length();	
	}
	
	@Test
	public void canGetTitle() throws InvalidGenomicCoordsException, InvalidColourException, IOException, ClassNotFoundException, InvalidRecordException, SQLException{
		String bam= "test_data/adjacent.bam";
		GenomicCoords gc= new GenomicCoords("chr7:1-100", samSeqDict, null);
		TrackReads tr= new TrackReads(bam, gc);
		
		tr.setNoFormat(true);
		tr.setTrackTag("aln.bam#1");
		
		assertEquals("aln.bam#1", tr.getTitle().trim());
		
	}
	
	@Test
	public void testOneLineStack() throws InvalidGenomicCoordsException, IOException, ClassNotFoundException, InvalidRecordException, SQLException, InvalidColourException{
		String bam= "test_data/adjacent.bam";
		int yMaxLines= 10;
		boolean bs= false;
		boolean noFormat= true;
		GenomicCoords gc= new GenomicCoords("chr7:1-50", samSeqDict, null);
			
		TrackReads tr= new TrackReads(bam, gc);
		tr.setBisulf(bs);
		tr.setNoFormat(noFormat);

		// NB: The success of this test depends on the screen width of eclipse
		String exp= 
		"AAAAAAAAAA           GGGGGGGGGG TTTTTTTTTT\n"+
		"          CCCCCCCCCC";
		tr.setyMaxLines(yMaxLines);
		assertEquals(exp, tr.printToScreen());
		
		System.out.println(tr.printToScreen());
		
		gc= new GenomicCoords("chr7:1-100", samSeqDict, null);
		tr= new TrackReads(bam, gc);
		tr.setBisulf(bs);
		tr.setNoFormat(noFormat);
		tr.setyMaxLines(yMaxLines);
		assertTrue(tr.printToScreen().split("\n").length == 1);
	}

	@Test
	public void test() throws InvalidGenomicCoordsException, IOException, ClassNotFoundException, InvalidRecordException, SQLException, InvalidColourException {
		String bam= "test_data/ds051.actb.bam";
		int yMaxLines= 50;
		boolean bs= false;
		boolean noFormat= true;
		GenomicCoords gc= new GenomicCoords("chr7:5566770-5566870", samSeqDict, fastaFile);
			
		TrackReads tr= new TrackReads(bam, gc);
		tr.setBisulf(bs);
		tr.setNoFormat(noFormat);
		tr.setyMaxLines(yMaxLines);
		System.out.println(tr.printToScreen());
	}

	@Test
	public void testNoReadsInRegion() throws InvalidGenomicCoordsException, IOException, ClassNotFoundException, InvalidRecordException, SQLException, InvalidColourException {
		String bam= "test_data/ds051.actb.bam";
		int yMaxLines= 50;
		boolean bs= false;
		boolean noFormat= true;
		GenomicCoords gc= new GenomicCoords("chr7:1-101", samSeqDict, fastaFile);
			
		TrackReads tr= new TrackReads(bam, gc);
		tr.setyMaxLines(yMaxLines);
		tr.setBisulf(bs);
		tr.setNoFormat(noFormat);
		assertEquals("", tr.printToScreen());
	}
	
	@Test
	public void canResetToZeroLargeWindow() throws IOException, InvalidGenomicCoordsException, ClassNotFoundException, InvalidRecordException, SQLException, InvalidColourException {
		// If the genomic window is too large do not process the bam file and return zero height track.
		GenomicCoords gc= new GenomicCoords("chr7:1-100000000", samSeqDict, fastaFile);
		TrackReads tr= new TrackReads("test_data/ds051.actb.bam", gc);
		assertEquals("", tr.printToScreen());
	}
}
