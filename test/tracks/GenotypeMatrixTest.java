package tracks;

import static org.junit.Assert.*;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.Arrays;
import java.util.List;

import org.junit.Before;
import org.junit.Test;

import coloring.Config;
import coloring.ConfigKey;
import coloring.Xterm256;
import exceptions.InvalidColourException;
import exceptions.InvalidConfigException;
import exceptions.InvalidGenomicCoordsException;
import exceptions.InvalidRecordException;
import samTextViewer.GenomicCoords;

public class GenotypeMatrixTest {

	@Before
	public void config() throws IOException, InvalidConfigException{
		new Config(null);
		new Xterm256();
	}

	@Test
	public void bigData()throws InvalidGenomicCoordsException, IOException, ClassNotFoundException, InvalidRecordException, SQLException, InvalidColourException {

		GenomicCoords gc= new GenomicCoords("1:572807-755079", 200, null, null);
		TrackIntervalFeature vcf= new TrackIntervalFeature("/Users/db291g/Downloads/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz", gc);
		vcf.setNoFormat(false);
		long t0= System.currentTimeMillis();
		System.err.println(vcf.printToScreen().length());
		long t1= System.currentTimeMillis();
		System.err.println(t1-t0);
	}	
	
	@Test
	public void canInitMatrix() throws InvalidGenomicCoordsException, IOException, ClassNotFoundException, InvalidRecordException, SQLException, InvalidColourException {
		
		GenomicCoords gc= new GenomicCoords("1:572807-755079", 80, null, null);
		TrackIntervalFeature vcf= new TrackIntervalFeature("test_data/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz", gc);
		List<IntervalFeature> linf = vcf.getIntervalFeatureList();

		GenotypeMatrix gm= new GenotypeMatrix();
		
		// No data at all
		assertTrue(gm.printToScreen(true).isEmpty());
		
		gm.makeMatrix(linf, 80);
		String x= gm.printToScreen(true);
		
		// Check sample name
		assertTrue(x.startsWith("HG00096"));
		
		String[] rows = x.split("\n");
		assertEquals(3, rows.length);
		assertEquals(80, rows[0].length());
		
		// Check genotype coding
		assertTrue(rows[0].contains("?")); // Missing genotype
		assertTrue(rows[0].contains(".")); // HOM ref
		assertTrue(rows[1].contains("E")); // HET
		assertTrue(rows[2].contains("0")); // HOM alt
	}

	@Test
	public void canSelectSamplesByRegex() throws Exception{

		GenomicCoords gc= new GenomicCoords("1:572807-755079", 80, null, null);
		TrackIntervalFeature vcf= new TrackIntervalFeature("test_data/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz", gc);
		List<IntervalFeature> linf = vcf.getIntervalFeatureList();

		GenotypeMatrix gm= new GenotypeMatrix();
		gm.makeMatrix(linf, 80);
		gm.setSelectSampleRegex("96|99");

		String x= gm.printToScreen(true);
		String[] rows = x.split("\n");
		assertEquals(2, rows.length);
		
		// Exclude all samples
		gm.setSelectSampleRegex("^$");
		assertTrue(gm.printToScreen(true).isEmpty());
	}
	
	@Test
	public void behaviourWithHaploid(){
		// TODO
	}

	@Test
	public void canHandleMissingGenotypes(){
		// TODO
	}
	
	@Test
	public void canHandleIntervalWithNoFeatures() throws InvalidGenomicCoordsException, IOException, ClassNotFoundException, InvalidRecordException, SQLException, InvalidColourException{
		GenomicCoords gc= new GenomicCoords("2:1-1000", 80, null, null);
		TrackIntervalFeature vcf= new TrackIntervalFeature("test_data/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz", gc);
		List<IntervalFeature> linf = vcf.getIntervalFeatureList();

		GenotypeMatrix gm= new GenotypeMatrix();
		gm.makeMatrix(linf, 80);
		String x= gm.printToScreen(true);
		assertTrue(x.isEmpty());
	}

	@Test
	public void canHandleVCFWithNoSamples() throws Exception{
		GenomicCoords gc= new GenomicCoords("1:1-10000000", 80, null, null);
		TrackIntervalFeature vcf= new TrackIntervalFeature("test_data/CHD.exon.2010_03.sites.vcf.gz", gc);
		List<IntervalFeature> linf = vcf.getIntervalFeatureList();
		assertTrue(linf.size() > 0); // Make sure we do have some features otherwise the test is meaningless.
		
		GenotypeMatrix gm= new GenotypeMatrix();
		gm.makeMatrix(linf, 80);
		String x= gm.printToScreen(true);
		assertTrue(x.isEmpty());
	}

}
