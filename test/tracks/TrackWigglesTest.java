package tracks;

import static org.junit.Assert.*;

import java.io.IOException;
import java.security.NoSuchAlgorithmException;
import java.util.List;

import org.broad.igv.bbfile.BBFileReader;
import org.broad.igv.bbfile.BigWigIterator;
import org.broad.igv.tdf.TDFUtils;
import org.junit.Test;

import exceptions.InvalidGenomicCoordsException;
import samTextViewer.GenomicCoords;

public class TrackWigglesTest {

	@Test
	public void canReadBigWigFromRemote() throws IOException{
		// String urlStr= "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeHaibTfbs/wgEncodeHaibTfbsA549Atf3V0422111Etoh02RawRep1.bigWig";
		String urlStr= "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeHaibTfbs/wgEncodeHaibTfbsA549Cebpbsc150V0422111RawRep1.bigWig";
		BBFileReader reader=new BBFileReader(urlStr);
		System.out.println(reader.getChromosomeNames());
		BigWigIterator iter = reader.getBigWigIterator("chr1", 1000000, "chr1", 2000000, true);
		while(iter.hasNext()){
			System.out.println(iter.next().getStartBase());
		}
		System.out.println("NEW");
		iter = reader.getBigWigIterator("chr10", 1000000, "chr10", 2000000, true);
			while(iter.hasNext()){
				System.out.println(iter.next().getStartBase());
			}
		reader.close();
	}
	
	@Test
	public void canGetDataColumnIndexForBedGraph() throws IOException, NoSuchAlgorithmException, InvalidGenomicCoordsException{
		
		String url= "test_data/test.bedGraph";
		int windowSize= 160;
		GenomicCoords gc= new GenomicCoords("chr1:1-30", null, windowSize, null);
		TrackWiggles tw= new TrackWiggles(url, gc, 5);
		assertEquals(0, tw.getScreenScores().get(0), 0.0001);
	}
	
	
	@Test
	public void canParseNonBGZFFile() throws IOException, InvalidGenomicCoordsException{
		
		String url= "test_data/test2.bedGraph";
		int windowSize= 160;
		GenomicCoords gc= new GenomicCoords("chr1:1-30", null, windowSize, null);
		TrackWiggles tw= new TrackWiggles(url, gc, 4);
				
	}
	
	@Test
	public void testYLimits() throws InvalidGenomicCoordsException, IOException{

		String url= "test_data/test.bedGraph.gz";
		int windowSize= 160;
		GenomicCoords gc= new GenomicCoords("chr1:1-30", null, windowSize, null);
		TrackWiggles tw= new TrackWiggles(url, gc, 4);
		tw.setYLimitMax(10.0);
		tw.setYLimitMin(-10.0);
		tw.setyMaxLines(10);
		String prof= tw.printToScreen();
		System.out.println(prof);
		
		tw.setTitle("bla");
		System.out.println(tw.getTitle());	
	}
	
	@Test
	public void testCloseToBorder() throws InvalidGenomicCoordsException, IOException{
		String url= "test_data/test.bedGraph.gz";
		int yMaxLines= 10;
		int windowSize= 160;
		GenomicCoords gc= new GenomicCoords("chr1:1-800", null, windowSize, null);
		TrackWiggles tw= new TrackWiggles(url, gc, 4);
		tw.setYLimitMax(Double.NaN);
		tw.setYLimitMin(Double.NaN);
		tw.setyMaxLines(yMaxLines);
		String prof= tw.printToScreen();
		System.out.println(prof);
	} 
	
	
	@Test 
	public void canPrintBedGraph() throws InvalidGenomicCoordsException, IOException{
		
		String url= "test_data/test.bedGraph.gz";
		int yMaxLines= 5;
		int windowSize= 22;
		GenomicCoords gc= new GenomicCoords("chr1", 1, 22, null, windowSize, null);
		TrackWiggles tw= new TrackWiggles(url, gc, 4);
		tw.setYLimitMax(Double.NaN);
		tw.setYLimitMin(Double.NaN);
		tw.setyMaxLines(yMaxLines);
		String prof= tw.printToScreen();
		System.out.println(prof);
		
		tw= new TrackWiggles("test_data/positive.bedGraph.gz", gc, 4);
		tw.setYLimitMax(Double.NaN);
		tw.setYLimitMin(Double.NaN);
		tw.setyMaxLines(5);
		prof= tw.printToScreen();
		System.out.println(prof);
		
		tw= new TrackWiggles("test_data/negative.bedGraph.gz", gc, 4);
		tw.setYLimitMax(Double.NaN);
		tw.setYLimitMin(Double.NaN);
		tw.setyMaxLines(5);
		// prof= tw.printToScreen();

		System.out.println(prof);
		
		gc= new GenomicCoords("chr1", 1, 52, null, 52, null);
		tw= new TrackWiggles("test_data/posNeg.bedGraph.gz", gc, 4);
		tw.setYLimitMax(Double.NaN);
		tw.setYLimitMin(Double.NaN);
		tw.setyMaxLines(14);
		System.out.println(tw.printToScreen());
	}
	
	// @Test
	public void canPrintWiggleTrack() throws InvalidGenomicCoordsException, IOException {
		
		// * Check big* are 0 or 1 based
		
		// String url= "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeCaltechRnaSeq/wgEncodeCaltechRnaSeqGm12892R2x75Th1014Il200SigRep3V4.bigWig";
		String url= "/Users/berald01/Downloads/wgEncodeCaltechRnaSeqGm12878R2x75Il400SigRep2V2.bigWig";
		
		int yMaxLines= 11;
		int windowSize= 101;
		GenomicCoords gc= new GenomicCoords("chrM", 1, 1000, null, windowSize, null);
		
		TrackWiggles tw= new TrackWiggles(url, gc, 4);
		// System.out.println(tw.printToScreen(yMaxLines));
		
		//System.out.println(tw.getMaxDepth());
		//System.out.println(tw.getScorePerDot());
	
	}
	
	// @Test
	public void canPrintFromTdf() throws IOException, InvalidGenomicCoordsException{

		GenomicCoords gc= new GenomicCoords("chr8", 1, 100, null, 100, null);
		String tdfFile= "test_data/hg18_var_sample.wig.v2.1.30.tdf";
		List<ScreenWiggleLocusInfo> screenLocInfo = 
		TDFUtils.tdfRangeToScreen(tdfFile, gc.getChrom(), gc.getFrom(), gc.getTo(), gc.getMapping());
		// assertEquals(0.925, screenLocInfo.get(1).getMeanScore(), 0.1);

	
		gc= new GenomicCoords("chrM:1-16000", null, 100, null);
		tdfFile= "/Volumes/My_Passport_for_Mac/tmp/rhh_hacat_0508-1406_FAIRE.tdf";
		screenLocInfo = TDFUtils.tdfRangeToScreen(tdfFile, gc.getChrom(), gc.getFrom(), gc.getTo(), gc.getMapping());
		int i= 1;
		for(ScreenWiggleLocusInfo x : screenLocInfo){
			//System.out.println(i + " " + x);
			i++;
		}

		TrackWiggles tw= new TrackWiggles(tdfFile, gc, 4);
		tw.setYLimitMax(Double.NaN);
		tw.setYLimitMin(Double.NaN);
		tw.setyMaxLines(40);
		System.out.println(tw.printToScreen());
		
	}

	
}
