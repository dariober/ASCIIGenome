package tracks;

import static org.junit.Assert.*;

import java.io.IOException;

import org.junit.Test;

import exceptions.InvalidGenomicCoordsException;
import samTextViewer.GenomicCoords;

public class TrackIntervalFeatureVCFTest {

	@Test
	public void canConstructTrack() throws InvalidGenomicCoordsException, IOException {
		String intervalFileName= "test_data/CHD.exon.2010_03.sites.vcf.gz";
		GenomicCoords gc= new GenomicCoords("1:1-10000000", null, 70, null);
		TrackIntervalFeatureVCF tif= new TrackIntervalFeatureVCF(intervalFileName, gc);
		tif.setNoFormat(true);
		
		IntervalFeatureSet x = tif.getIntervalFeatureSet();
		
		System.out.println(x.getFeaturesInInterval(gc.getChrom(), gc.getFrom(), gc.getTo()));
		
		//String exp= "||||||||||          ||||||||||                                        ";
		//assertEquals(exp, tif.printToScreen());

		//gc= new GenomicCoords("chr1:1-70", null, 35, null);
		//tif= new TrackIntervalFeature(intervalFileName, gc);
		//tif.setNoFormat(true);
		//exp= "|||||     |||||                    ";
		//assertEquals(exp, tif.printToScreen());	
	}

}
