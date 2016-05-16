package tracks;

import static org.junit.Assert.*;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import org.junit.Test;

import exceptions.InvalidGenomicCoordsException;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import samTextViewer.GenomicCoords;

public class TrackMethylationTest {

	public static List<SamRecordFilter> filters= new ArrayList<SamRecordFilter>();

	static SamReaderFactory srf=SamReaderFactory.make();
	static SamReader samReader= srf.open(new File("test_data/ds051.short.bam"));
	public static SAMSequenceDictionary samSeqDict= samReader.getFileHeader().getSequenceDictionary();

	public static String fastaFile= "test_data/chr7.fa";
	
	@Test
	public void canPrintMethylationProfile() throws IOException, InvalidGenomicCoordsException {

		int windowSize= 10;
		
		IndexedFastaSequenceFile faSeqFile = new IndexedFastaSequenceFile(new File("test_data/chr7.fa"));
		GenomicCoords gc= new GenomicCoords("chr7", 5566770, 5566870, samSeqDict, windowSize, fastaFile);
		
		int yMaxLines= 5;
		TrackCoverage tc= new TrackCoverage("test_data/ds051.short.bam", gc, filters, true);
		TrackMethylation tm= new TrackMethylation(tc.getFilename(), tc.getScreenLocusInfoList());
		tm.setyMaxLines(yMaxLines);
		tm.setNoFormat(true);
		String tmProfile= tm.printToScreen();
		System.out.println(tmProfile);
		System.out.println(tm.getScorePerDot());
	
		yMaxLines= 25;
		windowSize= 101;
		gc= new GenomicCoords("chr7", 5566770, 5566870, samSeqDict, windowSize, fastaFile);
		tc= new TrackCoverage("test_data/ds051.short.bam", gc, filters, true);
		tm= new TrackMethylation(tc.getFilename(), tc.getScreenLocusInfoList());
		tm.setyMaxLines(yMaxLines);
		tm.setNoFormat(true);
		tmProfile= tm.printToScreen();
		System.out.println(tmProfile);
		System.out.println(tm.getScorePerDot());
	
		long t0= System.currentTimeMillis();
		int i= 0;
		Random rand= new Random();
		while(i < 10000000){
			rand.nextFloat();
			i++;
		}
		long t1= System.currentTimeMillis();
		System.out.println(t1-t0);
		
	}

}
