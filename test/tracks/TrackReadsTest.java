package tracks;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.List;

import org.junit.BeforeClass;
import org.junit.Test;

import com.google.common.base.Stopwatch;

import coloring.Config;
import coloring.ConfigKey;
import coloring.Xterm256;
import exceptions.InvalidColourException;
import exceptions.InvalidCommandLineException;
import exceptions.InvalidConfigException;
import exceptions.InvalidGenomicCoordsException;
import exceptions.InvalidRecordException;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import net.sourceforge.argparse4j.impl.Arguments;
import samTextViewer.GenomicCoords;

public class TrackReadsTest {

    @BeforeClass
    public static void init() throws IOException, InvalidConfigException {
        new Config(null);
    }

	static SamReaderFactory srf=SamReaderFactory.make();
	static SamReader samReader= srf.open(new File("test_data/ds051.actb.bam"));
	public static SAMSequenceDictionary samSeqDict= samReader.getFileHeader().getSequenceDictionary();
	public static String fastaFile= "test_data/chr7.fa";

	@Test
	public void canShadeLowBaseQuality() throws InvalidGenomicCoordsException, InvalidColourException, ClassNotFoundException, IOException, InvalidRecordException, SQLException, InvalidCommandLineException, InvalidConfigException{
		
		String shade = Config.get(ConfigKey.shade_low_mapq);
		final Xterm256 xterm256= new Xterm256();
		String xshade = Integer.toString(xterm256.colorNameToXterm256(shade));
		
		GenomicCoords gc= new GenomicCoords("chr7:999-1041", 80, null, null);
		TrackReads tr= new TrackReads("test_data/missingReadSeq.bam", gc);
		System.err.println(tr.printToScreen());
		assertTrue(tr.printToScreen().trim().startsWith("[48;5;" + xshade));

		// Read with soft clipped bases
		gc= new GenomicCoords("chr7:9999-10050", 80, null, null);
		tr= new TrackReads("test_data/missingReadSeq.bam", gc);
		assertTrue(tr.printToScreen().contains(xshade));

		// Read with deletions and skipped bases
		gc= new GenomicCoords("chr7:19999-20050", 80, null, null);
		tr= new TrackReads("test_data/missingReadSeq.bam", gc);
		tr.printToScreen();
	}

	@Test
	public void canChangeReadColourOnRegex() throws InvalidGenomicCoordsException, IOException, ClassNotFoundException, InvalidRecordException, SQLException, InvalidColourException{
		GenomicCoords gc= new GenomicCoords("chr7:5566778-5566943", 80, null, null);
		TrackReads tr= new TrackReads("test_data/ds051.short.bam", gc);
		List<Argument> list= new ArrayList<Argument>();
		Argument re= new Argument("NCNNNCCC", "red1", false);
		list.add(re);
		tr.changeFeatureColor(list);
		assertTrue(tr.printToScreen().contains("196;"));
	}
	
	@Test
	public void canReadReadsWithMissingSequence() throws ClassNotFoundException, IOException, InvalidGenomicCoordsException, InvalidRecordException, SQLException, InvalidColourException{

		// Window size larger then 1 bp per column
		GenomicCoords gc= new GenomicCoords("chr7:1-1000", 80, null, null);
		TrackReads tr= new TrackReads("test_data/missingReadSeq.bam", gc);
		tr.setNoFormat(true);
		assertTrue(tr.printToScreen().trim().startsWith(">"));
		assertTrue(tr.printToScreen().trim().endsWith(">"));
		assertTrue(tr.printToScreen().trim().length() > 1);
		
		// Window size < 1 bp per column
		gc= new GenomicCoords("chr7:100-120", 80, null, null);
		tr= new TrackReads("test_data/missingReadSeq.bam", gc);
		tr.setNoFormat(true);
		assertTrue(tr.printToScreen().trim().startsWith("N"));
		assertTrue(tr.printToScreen().trim().endsWith("N"));
		assertTrue(tr.printToScreen().trim().length() > 1);

	}
	
//	@Test 
	public void testSpeed() throws InvalidGenomicCoordsException, IOException, ClassNotFoundException, InvalidRecordException, SQLException, InvalidColourException, InvalidConfigException{
		GenomicCoords gc= new GenomicCoords("chr7:5567700-5567872", 80, null, null);
		SamReader sr= srf.open(new File("test_data/MT.bam"));
//		SAMRecordIterator reads = sr.query("chr7", gc.getFrom(), gc.getTo(), false);
//		int n= 0;
//		long t0= System.currentTimeMillis();
//		while(reads.hasNext()){
//			SAMRecord rec= reads.next();
//			Cigar cigar= rec.getCigar();
//			List<AlignmentBlock> blocks = rec.getAlignmentBlocks();
//			rec.getReferencePositionAtReadPosition(0);
//			n++;
//			TextRead txr= new TextRead(rec, gc);
//			txr.getPrintableTextRead(false, true, false, gc.getBpPerScreenColumn());
//		}
//		long t1= System.currentTimeMillis();
//		System.out.println((t1-t0)/1000.0);
//		System.out.println(n);
		
		long t2= System.currentTimeMillis();
		TrackReads tr= new TrackReads("test_data/MT.bam", gc);
		long t3= System.currentTimeMillis();
		System.out.println(tr.getRecordsAsStrings().size());
		System.out.println((t3-t2)/1000.0);
		
	}
	
	@Test
	public void canShowReadsAsPairs() throws Exception{
		GenomicCoords gc= new GenomicCoords("chr7:1-80", 80, null, null);
		TrackReads tr= new TrackReads("test_data/pairs.sam", gc);
		tr.setNoFormat(true);
		// Still unpaired
		assertTrue(tr.printToScreen().startsWith("NANAN  gcgcgcgcgc  ntntn "));
		assertTrue(tr.printToScreen().contains("ATATATATAT  "));
		
		tr.setReadsAsPairs(true); // Switch on pairing
		
		// Properly paired
		assertTrue(tr.printToScreen().contains(" GGGGG~~~~~~~~~~~~~~~ggggg "));
		
		// Overlapping pair
		assertTrue(tr.printToScreen().contains("ATATATAgcgcgcgcgc "));
		System.err.println(tr.printToScreen());
		
		// Pair where one read is fully contained in the other
		assertTrue(tr.printToScreen().contains(" CCaaaaaaCC "));
		
		// Properly paired but mate is not in the window
		assertTrue(tr.printToScreen().contains("TTTTT  "));
		
		// Properly paired flag is set but mateAlignmentStart is not set. Treat as unpaired:
		assertTrue(tr.printToScreen().contains("NANAN  "));
		assertTrue(tr.printToScreen().contains(" ntntn "));
			
	}
	
	@Test
	public void canReturnReadsAsRawStrings() throws ClassNotFoundException, IOException, InvalidGenomicCoordsException, InvalidRecordException, SQLException{
		GenomicCoords gc= new GenomicCoords("chr7:5566000-5567000", 80, null, null);
		TrackReads tr= new TrackReads("test_data/ds051.short.bam", gc);
		
		assertEquals(22, tr.getRecordsAsStrings().size());

		// No read in interval
		gc= new GenomicCoords("chr1:1-1000", 80, null, null);
		tr= new TrackReads("test_data/ds051.short.bam", gc);
		assertEquals(0, tr.getRecordsAsStrings().size());
		
	}
	
	@Test
	public void canPrintReads() throws InvalidGenomicCoordsException, IOException, ClassNotFoundException, InvalidRecordException, SQLException, InvalidColourException, InvalidConfigException, InvalidCommandLineException{
		
		new Config(null);
		
		GenomicCoords gc= new GenomicCoords("chr7:5566000-5567000", 80, null, null);
		TrackReads tr= new TrackReads("test_data/ds051.short.bam", gc);
		tr.setNoFormat(true);
		
		tr.setPrintMode(PrintRawLine.OFF);
		String printable = tr.printLines();
		assertEquals("", printable); // Empty string if printing is off.
		
		tr.setPrintMode(PrintRawLine.CLIP);
		
		tr.setPrintRawLineCount(5);
		printable = tr.printLines();
		assertTrue(printable.length() > 100);
		assertEquals(5+1, printable.split("\n").length); // Expect 6 lines: 5 for reads and 1 for info header.
	}
	
	@Test
	public void canShowReadsInWindow() throws Exception{
		GenomicCoords gc= new GenomicCoords("chr7:5566000-5567000",80, samSeqDict, null);
		TrackReads tr= new TrackReads("test_data/ds051.short.bam", gc);
		tr.setNoFormat(true);
		tr.setyMaxLines(1000);
		assertTrue(tr.getTitle().contains("22")); // N. reads stacked in this interval before filtering
	}
	
	@Test
	public void canFilterReadsWithAwk() throws InvalidGenomicCoordsException, IOException, ClassNotFoundException, InvalidRecordException, SQLException, InvalidColourException{
		GenomicCoords gc= new GenomicCoords("chr7:5566000-5567000",80, samSeqDict, null);
		TrackReads tr= new TrackReads("test_data/ds051.short.bam", gc);
		tr.setNoFormat(true);
		tr.setyMaxLines(1000);
		assertEquals(22, tr.printToScreen().split("\n").length); // N. reads stacked in this interval before filtering		
		tr.setAwk("'$1 ~ \"NCNNNCCC\"'");
		assertEquals(6, tr.printToScreen().split("\n").length);
		assertTrue(tr.getTitle().contains("awk"));
		System.err.println(tr.getTitle());
	}

	@Test
	public void canFilterReadsWithGrep() throws InvalidGenomicCoordsException, IOException, ClassNotFoundException, InvalidRecordException, SQLException, InvalidColourException{
		GenomicCoords gc= new GenomicCoords("chr7:5566000-5567000",80, samSeqDict, null);
		TrackReads tr= new TrackReads("test_data/ds051.short.bam", gc);
		tr.setNoFormat(true);
		tr.setyMaxLines(1000);
		assertEquals(22, tr.printToScreen().split("\n").length); // N. reads stacked in this interval before filtering		
		tr.setShowHideRegex("NCNNNCCC", "\\t5566779\\t");
		assertEquals(4, tr.printToScreen().split("\n").length);
	}

	@Test
	public void canFilterReadsWithGrepAndAwk() throws InvalidGenomicCoordsException, IOException, ClassNotFoundException, InvalidRecordException, SQLException, InvalidColourException{
		GenomicCoords gc= new GenomicCoords("chr7:5566000-5567000",80, samSeqDict, null);
		TrackReads tr= new TrackReads("test_data/ds051.short.bam", gc);
		tr.setNoFormat(true);
		tr.setyMaxLines(1000);
		assertEquals(22, tr.printToScreen().split("\n").length); // N. reads stacked in this interval before filtering		
		tr.setShowHideRegex("NCNNNCCC", Track.HIDE_REGEX);
		tr.setAwk("'$4 != 5566779'");
		assertEquals(4, tr.printToScreen().split("\n").length);
	}

	@Test
	public void canShowReadCount() throws Exception{
		GenomicCoords gc= new GenomicCoords("chr7:5565600-5567600", 80, null, null);
		TrackReads tr= new TrackReads("test_data/ear045.oxBS.actb.bam", gc);
		
		// Same as: samtools view -F 4 -c ear045.oxBS.actb.bam chr7:5565600-5567600
		// AND excluding also reads fully soft-clipped
		assertTrue(tr.getTitle().contains("2436"));
	}
	
	@Test
	public void canGetTitle() throws InvalidGenomicCoordsException, InvalidColourException, IOException, ClassNotFoundException, InvalidRecordException, SQLException{
		String bam= "test_data/adjacent.bam";
		GenomicCoords gc= new GenomicCoords("chr7:1-100",80, samSeqDict, null);
		TrackReads tr= new TrackReads(bam, gc);
		
		tr.setNoFormat(true);
		tr.setTrackTag("aln.bam#1");
		
		assertTrue(tr.getTitle().trim().startsWith("aln.bam#1"));
		
	}
	
	@Test
	public void testOneLineStack() throws InvalidGenomicCoordsException, IOException, ClassNotFoundException, InvalidRecordException, SQLException, InvalidColourException{
		String bam= "test_data/adjacent.bam";
		int yMaxLines= 10;
		boolean bs= false;
		boolean noFormat= true;
		GenomicCoords gc= new GenomicCoords("chr7:1-50",80, samSeqDict, null);
			
		TrackReads tr= new TrackReads(bam, gc);
		tr.setBisulf(bs);
		tr.setNoFormat(noFormat);

		// NB: The success of this test depends on the screen width of eclipse
//		String exp= 
//		"AAAAAAAAAA           GGGGGGGGGG TTTTTTTTTT \n"+
//		"          CCCCCCCCCC                       ";
		tr.setyMaxLines(yMaxLines);
		assertEquals(2, tr.printToScreen().split("\n").length); // Two lines
		
		System.out.println(tr.printToScreen());
		
		gc= new GenomicCoords("chr7:1-100",80, samSeqDict, null);
		tr= new TrackReads(bam, gc);
		tr.setBisulf(bs);
		tr.setNoFormat(noFormat);
		tr.setyMaxLines(yMaxLines);
		assertTrue(tr.printToScreen().split("\n").length == 1);
	}

	@Test
	public void test() throws InvalidGenomicCoordsException, IOException, ClassNotFoundException, InvalidRecordException, SQLException, InvalidColourException, InvalidConfigException {
		
		new Config(null);
		
		String bam= "test_data/ds051.actb.bam";
		int yMaxLines= 50;
		boolean bs= false;
		GenomicCoords gc= new GenomicCoords("chr7:5566770-5566870",80, samSeqDict, fastaFile);
			
		TrackReads tr= new TrackReads(bam, gc);
		tr.setBisulf(bs);
		tr.setNoFormat(true);
		tr.setyMaxLines(yMaxLines);
		System.out.print(tr.printToScreen());
	}

	@Test
	public void testNoReadsInRegion() throws InvalidGenomicCoordsException, IOException, ClassNotFoundException, InvalidRecordException, SQLException, InvalidColourException {
		String bam= "test_data/ds051.actb.bam";
		int yMaxLines= 50;
		boolean bs= false;
		boolean noFormat= true;
		GenomicCoords gc= new GenomicCoords("chr7:1-101",80, samSeqDict, fastaFile);
			
		TrackReads tr= new TrackReads(bam, gc);
		tr.setyMaxLines(yMaxLines);
		tr.setBisulf(bs);
		tr.setNoFormat(noFormat);
		assertEquals("", tr.printToScreen());
	}
	
	@Test
	public void canResetToZeroLargeWindow() throws IOException, InvalidGenomicCoordsException, ClassNotFoundException, InvalidRecordException, SQLException, InvalidColourException {
		// If the genomic window is too large do not process the bam file and return zero height track.
		GenomicCoords gc= new GenomicCoords("chr7:1-100000000",80, samSeqDict, fastaFile);
		TrackReads tr= new TrackReads("test_data/ds051.actb.bam", gc);
		assertEquals("", tr.printToScreen());
	}
	
	@Test
	public void canConstructFromUnsortedInput() throws IOException, InvalidGenomicCoordsException, ClassNotFoundException, InvalidRecordException, SQLException, InvalidColourException {
		// If the genomic window is too large do not process the bam file and return zero height track.
		GenomicCoords gc= new GenomicCoords("chr7:1-100000000",80, samSeqDict, fastaFile);
		TrackReads tr= new TrackReads("test_data/ds051.noindex.sam", gc);
		assertEquals("", tr.printToScreen());
	}
	
}
