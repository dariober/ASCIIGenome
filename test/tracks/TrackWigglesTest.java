package tracks;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.IOException;
import java.security.NoSuchAlgorithmException;
import java.sql.SQLException;
import java.util.List;

import org.broad.igv.bbfile.BBFileReader;
import org.broad.igv.bbfile.BigWigIterator;
import org.broad.igv.tdf.TDFUtils;
import org.junit.Test;

import exceptions.InvalidColourException;
import exceptions.InvalidGenomicCoordsException;
import exceptions.InvalidRecordException;
import samTextViewer.GenomicCoords;

public class TrackWigglesTest {

	@Test
	public void canPrintChromosomeNames() throws InvalidGenomicCoordsException, IOException, ClassNotFoundException, InvalidRecordException, SQLException{

		GenomicCoords gc= new GenomicCoords("chr7:5540000-5570000", null, null);
		TrackWiggles tw= new TrackWiggles("test_data/hg18_var_sample.wig.v2.1.30.tdf", gc, 4);
		assertTrue(tw.getChromosomeNames().size() > 10);
		
		tw= new TrackWiggles("test_data/test.bedGraph", gc, 4);
		assertTrue(tw.getChromosomeNames().size() > 0);
		
		tw= new TrackWiggles("http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeHaibTfbs/wgEncodeHaibTfbsA549Cebpbsc150V0422111RawRep1.bigWig", gc, 4);
		assertTrue(tw.getChromosomeNames().size() > 10);
	}
	
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
	public void canGetDataColumnIndexForBedGraph() throws IOException, NoSuchAlgorithmException, InvalidGenomicCoordsException, InvalidRecordException, ClassNotFoundException, SQLException{
		
		String url= "test_data/test.bedGraph";
		GenomicCoords gc= new GenomicCoords("chr1:1-30", null, null);
		TrackWiggles tw= new TrackWiggles(url, gc, 5);
		assertEquals(0, tw.getScreenScores().get(0), 0.0001);
	}
	
	
	@Test
	public void canParseNonBGZFFile() throws IOException, InvalidGenomicCoordsException, InvalidRecordException, ClassNotFoundException, SQLException{
		
		String url= "test_data/test2.bedGraph";
		GenomicCoords gc= new GenomicCoords("chr1:1-30", null, null);
		TrackWiggles tw= new TrackWiggles(url, gc, 4);
				
	}
	
	@Test
	public void testYLimits() throws InvalidGenomicCoordsException, IOException, InvalidColourException, InvalidRecordException, ClassNotFoundException, SQLException{

		String url= "test_data/test.bedGraph.gz";
		GenomicCoords gc= new GenomicCoords("chr1:1-30", null, null);
		TrackWiggles tw= new TrackWiggles(url, gc, 4);
		tw.setYLimitMax(10.0);
		tw.setYLimitMin(-10.0);
		tw.setyMaxLines(10);
		String prof= tw.printToScreen();
		System.out.println(prof);
		
	}
	
	@Test
	public void testCloseToBorder() throws InvalidGenomicCoordsException, InvalidColourException, IOException, InvalidRecordException, ClassNotFoundException, SQLException{
		String url= "test_data/test.bedGraph.gz";
		int yMaxLines= 10;
		GenomicCoords gc= new GenomicCoords("chr1:1-800", null, null);
		TrackWiggles tw= new TrackWiggles(url, gc, 4);
		tw.setYLimitMax(Double.NaN);
		tw.setYLimitMin(Double.NaN);
		tw.setyMaxLines(yMaxLines);
		String prof= tw.printToScreen();
		System.out.println(prof);
	} 
	
	
	@Test 
	public void canPrintBedGraph() throws InvalidGenomicCoordsException, IOException, InvalidRecordException, ClassNotFoundException, SQLException, InvalidColourException{
		
		String url= "test_data/test.bedGraph.gz";
		int yMaxLines= 5;
		GenomicCoords gc= new GenomicCoords("chr1:1-22", null, null);
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
		
		gc= new GenomicCoords("chr1:1-52", null, null);
		tw= new TrackWiggles("test_data/posNeg.bedGraph.gz", gc, 4);
		tw.setYLimitMax(Double.NaN);
		tw.setYLimitMin(Double.NaN);
		tw.setyMaxLines(14);
		System.out.println(tw.printToScreen());
	}
	
	// @Test
	public void canPrintWiggleTrack() throws InvalidGenomicCoordsException, IOException, InvalidRecordException, ClassNotFoundException, SQLException {
		
		// * Check big* are 0 or 1 based
		
		// String url= "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeCaltechRnaSeq/wgEncodeCaltechRnaSeqGm12892R2x75Th1014Il200SigRep3V4.bigWig";
		String url= "/Users/berald01/Downloads/wgEncodeCaltechRnaSeqGm12878R2x75Il400SigRep2V2.bigWig";
		
		GenomicCoords gc= new GenomicCoords("chrM:1-1000", null, null);
		
		TrackWiggles tw= new TrackWiggles(url, gc, 4);
		// System.out.println(tw.printToScreen(yMaxLines));
		
		//System.out.println(tw.getMaxDepth());
		//System.out.println(tw.getScorePerDot());
	
	}
	
	// @Test
	public void canPrintFromTdf() throws IOException, InvalidGenomicCoordsException, InvalidRecordException, ClassNotFoundException, SQLException, InvalidColourException{

		GenomicCoords gc= new GenomicCoords("chr8:1-100", null, null);
		String tdfFile= "test_data/hg18_var_sample.wig.v2.1.30.tdf";
		List<ScreenWiggleLocusInfo> screenLocInfo = 
		TDFUtils.tdfRangeToScreen(tdfFile, gc.getChrom(), gc.getFrom(), gc.getTo(), gc.getMapping());
		// assertEquals(0.925, screenLocInfo.get(1).getMeanScore(), 0.1);

	
		gc= new GenomicCoords("chrM:1-16000", null, null);
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

	@Test
	/** Snippet to extract totalCount from TDF, useful for normalizing signal. 
	 * */
	public void canNomrmalizeTDFtoRPM() throws InvalidGenomicCoordsException, IOException, InvalidRecordException, ClassNotFoundException, SQLException{

		System.out.println("START");
		GenomicCoords gc= new GenomicCoords("chr7:5540000-5570000", null, null);
		TrackWiggles tw= new TrackWiggles("test_data/ear045.oxBS.actb.tdf", gc, 4);
		Double raw= tw.getScreenScores().get(0);
		tw.setRpm(true);
		Double rpm= tw.getScreenScores().get(0);
		assertTrue(rpm > raw);
	}
		
}
