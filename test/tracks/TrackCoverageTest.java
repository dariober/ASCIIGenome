package tracks;

import static org.junit.Assert.*;

import java.io.File;
import java.io.IOException;
import java.net.URL;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import org.junit.Test;

import exceptions.InvalidGenomicCoordsException;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import samTextViewer.GenomicCoords;
import samTextViewer.SamLocusIterator;
import tracks.TrackCoverage;

public class TrackCoverageTest {
	
	public static List<SamRecordFilter> filters= new ArrayList<SamRecordFilter>();
	
	static SamReaderFactory srf=SamReaderFactory.make();
	static SamReader sr= srf.open(new File("test_data/ds051.short.bam"));
	public static SAMSequenceDictionary samSeqDict= sr.getFileHeader().getSequenceDictionary();
	public static String fastaFile= "test_data/chr7.fa";
	
	//@Test
	public void testSpeedSamLocIter() throws IOException{
		
		SamReader samReader= srf.open(new File("/Volumes/My_Passport_for_Mac/tmp/rhh147-148_untreat_14102014_atac_hacat.bam"));
		SAMFileHeader fh= samReader.getFileHeader();

		// FAST
		// IntervalList il= new IntervalList(fh);
		// int x= 5000000;
		// while(x < 10000000){
		//	il.add(new Interval("chr7", x, x));
		//	x += 10;
		// }
		//IntervalList il= new IntervalList(fh);
		// SamLocusIterator samLocIter= new SamLocusIterator(samReader, il); // SLOW
		//System.out.println(il.size());

		IntervalList il= new IntervalList(fh);
		il.add(new Interval("chr7", 5000000, 60000000));
		SamLocusIterator samLocIter= new SamLocusIterator(samReader, il); // FAST
		System.out.println("SamLocIter DONE");
		
		long t0= System.currentTimeMillis();
		Iterator<samTextViewer.SamLocusIterator.LocusInfo> iter= samLocIter.iterator(); // FAST
		long t1= System.currentTimeMillis();
		System.out.println("Iterator done in: " + (t1- t0));
		
		// FAST
		int i= 0;
		long nbp= 0;
		while(iter.hasNext()){
			i++;
			samTextViewer.SamLocusIterator.LocusInfo locusInfo= iter.next();
			nbp += locusInfo.getRecordAndPositions().size();
		}
		long t2= System.currentTimeMillis();
		System.out.println(i + " loci Done in: " + (t2- t1) + " ms; counted " + nbp);
		samLocIter.close();
		samReader.close();
	}
	
	@Test
	public void testRPMnorm() throws InvalidGenomicCoordsException, IOException{
		
		int yMaxLines= 11;
		int windowSize= 101;

		GenomicCoords gc= new GenomicCoords("chr7", 5566770, 5566870, samSeqDict, windowSize, fastaFile);
		TrackCoverage tc= new TrackCoverage("test_data/ds051.short.bam", gc, filters, false);
		tc.setyMaxLines(yMaxLines);
		tc.setRpm(true);
		// assertEquals(1000000, tc.getMaxDepth(), 0.1);
	}
	
	@Test
	public void canPrintCoverageTrack() throws IOException, InvalidGenomicCoordsException {
		int yMaxLines= 11;
		int windowSize= 101;

		GenomicCoords gc= new GenomicCoords("chr7", 5566770, 5566870, samSeqDict, windowSize, fastaFile);
		TrackCoverage tc= new TrackCoverage("test_data/ds051.short.bam", gc, filters, false);
		tc.setyMaxLines(yMaxLines);
		System.out.println(gc.toString());
		System.out.println(tc.printToScreen());
		// System.out.println(tc.getScorePerDot());
	}

	@Test
	public void canReadBAMFromURL() throws IOException, InvalidGenomicCoordsException {
		
		String urlStr= "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeHaibTfbs/wgEncodeHaibTfbsA549Atf3V0422111Etoh02AlnRep1.bam";

		long t0= System.currentTimeMillis();
		
		SamReaderFactory srf=SamReaderFactory.make();
		srf.validationStringency(ValidationStringency.SILENT);
		SamReader samReader= SamReaderFactory.makeDefault().open(SamInputResource.of(new URL(urlStr)).index(new URL(urlStr + ".bai")));
		// SamReader samReader= srf.open(SamInputResource.of(new File("test_data/ear045.oxBS.actb.bam")).index(new File(urlStr + ".bai")));
		long t1= System.currentTimeMillis();
		System.out.println(t1-t0);
		SAMRecordIterator iter= samReader.iterator();
		int n= 0;
		while(iter.hasNext()){
			SAMRecord rec= iter.next();
			if(n % 10000 == 0){
				System.out.println(n);
			}
			n++;
		}
		samReader.close();
		
		// int windowSize= 101;
		// GenomicCoords gc= new GenomicCoords("chr7", 1, 10000, samSeqDict, windowSize, fastaFile);
		// TrackCoverage tc= new TrackCoverage(urlStr, gc, filters, false);
	}
	
	@Test
	public void canPrintCoverageTrackWithZeroHeight() throws IOException, InvalidGenomicCoordsException {
		int yMaxLines= 0;
		int windowSize= 10;

		GenomicCoords gc= new GenomicCoords("chr7", 5566770, 5566870, samSeqDict, windowSize, fastaFile);
		TrackCoverage tc= new TrackCoverage("test_data/ds051.short.bam", gc, filters, false);
		tc.setyMaxLines(yMaxLines);
		assertEquals("", tc.printToScreen());
	}
	
	@Test
	public void canResetToZeroLargeWindow() throws IOException, InvalidGenomicCoordsException {
		// If the genomic window is too large do not process the bam file and return zero height track.
		GenomicCoords gc= new GenomicCoords("chr7", 1, 100000000, samSeqDict, 100, fastaFile);
		TrackCoverage tc= new TrackCoverage("test_data/ds051.actb.bam", gc, filters, false);
		assertEquals(0, tc.getScreenLocusInfoList().size());
		assertTrue(tc.printToScreen().startsWith("Track"));
	}
	
	@Test 
	public void testWithZeroReadsRegion() throws InvalidGenomicCoordsException, IOException{
		
		GenomicCoords  gc= new GenomicCoords("chr7", 1, 1000, samSeqDict, 20, null);
		TrackCoverage tc= new TrackCoverage("test_data/ds051.short.bam", gc, filters, false);
		tc.setyMaxLines(2);
		assertEquals(
				   "                    \n"                   
				+  "____________________", tc.printToScreen());
	}

	@Test
	public void canPrintCoverage() throws IOException, InvalidGenomicCoordsException {
		int yMaxLines= 11;
		int windowSize= 101;

		GenomicCoords gc= new GenomicCoords("chr7:5568363-5568390", samSeqDict, windowSize, fastaFile);
		TrackCoverage tc= new TrackCoverage("test_data/ds051.actb.bam", gc, filters, false);
		tc.setyMaxLines(yMaxLines);
		System.out.println(gc.toString());
		System.out.println(tc.printToScreen());
	}
}
