package tracks;

import static org.junit.Assert.*;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.Random;

import org.junit.Test;

import exceptions.InvalidGenomicCoordsException;
import exceptions.InvalidRecordException;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import samTextViewer.GenomicCoords;

public class TrackMethylationTest {

	//public static List<SamRecordFilter> filters= new ArrayList<SamRecordFilter>();

	static SamReaderFactory srf=SamReaderFactory.make();
	static SamReader samReader= srf.open(new File("test_data/ds051.short.bam"));
	public static SAMSequenceDictionary samSeqDict= samReader.getFileHeader().getSequenceDictionary();

	public static String fastaFile= "test_data/chr7.fa";
	
	@Test
	public void canPrintMethylationProfile() throws IOException, InvalidGenomicCoordsException, ClassNotFoundException, InvalidRecordException, SQLException {

		// IndexedFastaSequenceFile faSeqFile = new IndexedFastaSequenceFile(new File("test_data/chr7.fa"));
		// faSeqFile.close();
		GenomicCoords gc= new GenomicCoords("chr7:5566770-5566870", samSeqDict, fastaFile);
		
		int yMaxLines= 5;
		TrackCoverage tc= new TrackCoverage("test_data/ds051.short.bam", gc, true);
		TrackMethylation tm= new TrackMethylation(tc.getFilename(), tc.getScreenLocusInfoList());
		tm.setyMaxLines(yMaxLines);
		tm.setNoFormat(true);
		String tmProfile= tm.printToScreen();
		System.out.println(tmProfile);
		// System.out.println(tm.getScorePerDot());
	
		yMaxLines= 25;
		gc= new GenomicCoords("chr7:5566770-5566870", samSeqDict, fastaFile);
		tc= new TrackCoverage("test_data/ds051.short.bam", gc, true);
		tm= new TrackMethylation(tc.getFilename(), tc.getScreenLocusInfoList());
		tm.setyMaxLines(yMaxLines);
		tm.setNoFormat(true);
		tmProfile= tm.printToScreen();
		System.out.println(tmProfile);
		// System.out.println(tm.getScorePerDot());
	
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
