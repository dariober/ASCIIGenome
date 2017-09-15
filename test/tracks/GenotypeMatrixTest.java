package tracks;

import static org.junit.Assert.*;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.List;

import org.junit.Before;
import org.junit.Test;

import com.google.common.base.Splitter;

import coloring.Config;
import coloring.Xterm256;
import exceptions.InvalidColourException;
import exceptions.InvalidConfigException;
import exceptions.InvalidGenomicCoordsException;
import exceptions.InvalidRecordException;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import samTextViewer.GenomicCoords;

public class GenotypeMatrixTest {

	@Before
	public void config() throws IOException, InvalidConfigException{
		new Config(null);
		new Xterm256();
	}

	@Test
	public void bigData()throws InvalidGenomicCoordsException, IOException, ClassNotFoundException, InvalidRecordException, SQLException, InvalidColourException {

		VCFFileReader reader = new VCFFileReader(new File("test_data/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz"));
		VCFHeader vcfHeader = reader.getFileHeader();
		reader.close();
		GenomicCoords gc= new GenomicCoords("1:5934301-6000000", 80, null, null);
		TrackIntervalFeature vcf= new TrackIntervalFeature("test_data/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz", gc);

		List<IntervalFeature> linf = vcf.getIntervalFeatureList();

		GenotypeMatrix gm= new GenotypeMatrix();
		gm.setJsScriptFilter("{GT} == 10 || 1>2");
		gm.printToScreen(true, linf, 80, vcfHeader);
	}	

	@Test
	public void canFilterWithJavascript() throws InvalidGenomicCoordsException, IOException, ClassNotFoundException, InvalidRecordException, SQLException, InvalidColourException, IOException {
		
		VCFFileReader reader = new VCFFileReader(new File("test_data/info_formats.vcf.gz"));
		VCFHeader vcfHeader = reader.getFileHeader();
		reader.close();
		
		GenomicCoords gc= new GenomicCoords("1:17822074-17822184", 80, null, null);
		TrackIntervalFeature vcf= new TrackIntervalFeature("test_data/info_formats.vcf.gz", gc);

		List<IntervalFeature> linf = vcf.getIntervalFeatureList();

		GenotypeMatrix gm= new GenotypeMatrix();
		
		// No filter
		assertEquals(2, gm.printToScreen(true, linf, 80, vcfHeader).split("\n").length);

		// Filter by array in FORMAT
		gm.setJsScriptFilter("{FF}[1] > 0");
		String y= gm.printToScreen(true, linf, 80, vcfHeader);
		assertTrue(y.contains("sample1"));
		assertTrue( ! y.contains("sample2") );
		
		// One sample satisfies DP at at least one record.
		gm.setJsScriptFilter("{DP} >= 100");
		String x= gm.printToScreen(true, linf, 80, vcfHeader);
		assertEquals(1, x.split("\n").length);

		// Filter by array in INFO
		gm.setJsScriptFilter("{XA}[0] > 0");
		x= gm.printToScreen(true, linf, 80, vcfHeader);
		assertTrue(x.contains("sample1"));

		// Exclude all samples.
		gm.setJsScriptFilter("{DP} >= 10000");
		x= gm.printToScreen(true, linf, 80, vcfHeader);
		assertTrue(x.isEmpty());

		// Exclude all sample II: Non sense comparison
		gm.setJsScriptFilter("{POS} == 'G' && {DP} > 5");
		x= gm.printToScreen(true, linf, 80, vcfHeader);
		assertTrue(x.isEmpty());

		// Remove filter
		gm.setJsScriptFilter("1 > 2"); // First remove all samples
		x= gm.printToScreen(true, linf, 80, vcfHeader); 
		assertTrue(x.isEmpty());
		
		gm.setJsScriptFilter(null); // Then remove filter and check samples are back
		x= gm.printToScreen(true, linf, 80, vcfHeader);
		assertTrue( ! x.isEmpty());
		
		gm.setJsScriptFilter(""); // Empty also works to remove filter
		x= gm.printToScreen(true, linf, 80, vcfHeader);
		assertTrue( ! x.isEmpty());
		
		// Handle missing values
		gm.setJsScriptFilter("{ID} == '.'"); // Use literal '.', not the null value.
		x= gm.printToScreen(true, linf, 80, vcfHeader);
		assertTrue( x.trim().length() > 50);

		// Filter by FLAG field type
		gm.setJsScriptFilter("{XB} && {DP} == 100");
		x= gm.printToScreen(true, linf, 80, vcfHeader);
		assertTrue(x.contains("sample2"));
		assertTrue( ! x.contains("sample1"));

		// Filter by ALT: Note array slicing
		// NB: Although all samples satisfy ALT, only sample2 also satisfies DP
		gm.setJsScriptFilter("{ALT}[0] == 'G' && {DP} > 5");
		x= gm.printToScreen(true, linf, 80, vcfHeader);
		assertEquals(1, x.split("\n").length);
		assertTrue(x.contains("sample2"));

		// Same using POS
		gm.setJsScriptFilter("{POS} == 17822094 && {DP} > 5");
		x= gm.printToScreen(true, linf, 80, vcfHeader);
		assertEquals(1, x.split("\n").length);
		assertTrue(x.contains("sample2"));
	}

	@Test
	public void genotypeAsNumericInJS() throws ClassNotFoundException, IOException, InvalidGenomicCoordsException, InvalidRecordException, SQLException, InvalidColourException, IOException{
		VCFFileReader reader = new VCFFileReader(new File("test_data/info_formats.vcf.gz"));
		VCFHeader vcfHeader = reader.getFileHeader();
		reader.close();
		GenomicCoords gc= new GenomicCoords("1:17822074-17822184", 80, null, null);
		TrackIntervalFeature vcf= new TrackIntervalFeature("test_data/info_formats.vcf.gz", gc);
		List<IntervalFeature> linf = vcf.getIntervalFeatureList();

		GenotypeMatrix gm= new GenotypeMatrix();
		
		gm.setJsScriptFilter("{GT} == '1/0'");
		String x= gm.printToScreen(true, linf, 80, vcfHeader);
		assertTrue(x.contains("sample1"));
		assertTrue( ! x.contains("sample2") );
		
		gm.setJsScriptFilter("{GT} == './.'");
		x= gm.printToScreen(true, linf, 80, vcfHeader);
		assertTrue( x.contains("sample1"));
		assertTrue( ! x.contains("sample2") );
	}
	
	@Test
	public void genotypeKetwordInJS() throws ClassNotFoundException, IOException, InvalidGenomicCoordsException, InvalidRecordException, SQLException, InvalidColourException, IOException{

		VCFFileReader reader = new VCFFileReader(new File("test_data/info_formats.vcf.gz"));
		VCFHeader vcfHeader = reader.getFileHeader();
		reader.close();
		GenomicCoords gc= new GenomicCoords("1:17822074-17822184", 80, null, null);
		TrackIntervalFeature vcf= new TrackIntervalFeature("test_data/info_formats.vcf.gz", gc);
		List<IntervalFeature> linf = vcf.getIntervalFeatureList();

		GenotypeMatrix gm= new GenotypeMatrix();
		
		gm.setJsScriptFilter("{HOM}");
		String x= gm.printToScreen(true, linf, 80, vcfHeader);
		assertTrue(x.contains("sample2"));
		assertTrue( ! x.contains("sample1") );
		
		gm.setJsScriptFilter("{NO_CALL}");
		x= gm.printToScreen(true, linf, 80, vcfHeader);
		assertTrue(x.contains("sample1"));
		assertTrue( ! x.contains("sample2") );
		
		gm.setJsScriptFilter("{HET_NON_REF}");
		x= gm.printToScreen(true, linf, 80, vcfHeader);
		assertTrue(x.contains("sample2"));
		assertTrue( ! x.contains("sample1") );
		
		gm.setJsScriptFilter("{CALLED} && {POS} == 17822094");
		x= gm.printToScreen(true, linf, 80, vcfHeader);
		assertTrue(x.contains("sample2"));
		assertTrue( ! x.contains("sample1") );
		
	}
	
	@Test
	public void invalidJSscripts() throws ClassNotFoundException, IOException, InvalidGenomicCoordsException, InvalidRecordException, SQLException, InvalidColourException, IOException{

		VCFFileReader reader = new VCFFileReader(new File("test_data/info_formats.vcf.gz"));
		VCFHeader vcfHeader = reader.getFileHeader();
		reader.close();
		GenomicCoords gc= new GenomicCoords("1:17822074-17822184", 80, null, null);
		TrackIntervalFeature vcf= new TrackIntervalFeature("test_data/info_formats.vcf.gz", gc);
		List<IntervalFeature> linf = vcf.getIntervalFeatureList();

		GenotypeMatrix gm= new GenotypeMatrix();
		
		// JS script must return boolean
		gm.setJsScriptFilter("{DP} > 5 && 10 + 3");
		gm.printToScreen(true, linf, 80, vcfHeader);
		assertNull(gm.getJsScriptFilter()); // After faulty script reset to null. 

		// Invalid JS syntax. E.g. when {TAG} does not exist.
		gm.setJsScriptFilter("{FOOBAR} > 10");
		gm.printToScreen(true, linf, 80, vcfHeader);
		assertNull(gm.getJsScriptFilter());
		
		// After exception the invalid filter has been removed
		gm.printToScreen(true, linf, 80, vcfHeader);
	}
	
	@Test
	public void canFilterGenotypeWithJS() throws IOException, ClassNotFoundException, IOException, InvalidGenomicCoordsException, InvalidRecordException, SQLException, InvalidColourException{

		VCFFileReader reader = new VCFFileReader(new File("test_data/info_formats.vcf.gz"));
		VCFHeader vcfHeader = reader.getFileHeader();
		reader.close();
		GenomicCoords gc= new GenomicCoords("1:17822074-17822184", 80, null, null);
		TrackIntervalFeature vcf= new TrackIntervalFeature("test_data/info_formats.vcf.gz", gc);
		List<IntervalFeature> linf = vcf.getIntervalFeatureList();

		GenotypeMatrix gm= new GenotypeMatrix();
		
		// Missing alleles
		gm.setJsScriptFilter("{GT} == './.'");
		String x= gm.printToScreen(true, linf, 80, vcfHeader);
		assertTrue(x.contains("sample1"));
		assertTrue( ! x.contains("sample2"));
	}
	
	@Test
	public void canInitMatrix() throws IOException, InvalidGenomicCoordsException, IOException, ClassNotFoundException, InvalidRecordException, SQLException, InvalidColourException {
		
		VCFFileReader reader = new VCFFileReader(new File("test_data/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz"));
		VCFHeader vcfHeader = reader.getFileHeader();
		reader.close();
		
		GenomicCoords gc= new GenomicCoords("1:572807-755079", 80, null, null);
		TrackIntervalFeature vcf= new TrackIntervalFeature("test_data/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz", gc);

		List<IntervalFeature> linf = vcf.getIntervalFeatureList();

		GenotypeMatrix gm= new GenotypeMatrix();
		
		// No feature at all
		assertTrue(gm.printToScreen(true, new ArrayList<IntervalFeature>(), 80, vcfHeader).isEmpty());
		
		String x= gm.printToScreen(true, linf, 80, vcfHeader);
		
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
	public void overalappingSymbols() throws IOException, InvalidGenomicCoordsException, IOException, ClassNotFoundException, InvalidRecordException, SQLException, InvalidColourException {
		
		VCFFileReader reader = new VCFFileReader(new File("test_data/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz"));
		VCFHeader vcfHeader = reader.getFileHeader();
		reader.close();
		
		// Same genotype: no ambiguity.
		GenomicCoords gc= new GenomicCoords("1:199882-200100", 80, null, null);
		TrackIntervalFeature vcf= new TrackIntervalFeature("test_data/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz", gc);
		List<IntervalFeature> linf = vcf.getIntervalFeatureList();
		
		GenotypeMatrix gm= new GenotypeMatrix();
		String x= gm.printToScreen(false, linf, 80, vcfHeader);
		assertTrue(x.contains("O")); // Two genotype at the same position: Both the same
		assertTrue(x.contains("*")); // Different
	}
	
	@Test
	public void canSelectSamplesByRegex() throws Exception{

		GenomicCoords gc= new GenomicCoords("1:572807-755079", 80, null, null);
		TrackIntervalFeature vcf= new TrackIntervalFeature("test_data/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz", gc);
		List<IntervalFeature> linf = vcf.getIntervalFeatureList();

		GenotypeMatrix gm= new GenotypeMatrix();
		gm.setSelectSampleRegex("96|99");
		// gm.makeMatrix(linf, 80, null);

		String x= gm.printToScreen(true, linf, 80, null);
		String[] rows = x.split("\n");
		assertEquals(2, rows.length);
		
		// Exclude all samples
		gm.setSelectSampleRegex("^$");
		// gm.makeMatrix(linf, 80, null);
		assertTrue(gm.printToScreen(true, linf, 80, null).isEmpty());
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
	public void canHandleIntervalWithNoFeatures() throws IOException, InvalidGenomicCoordsException, IOException, ClassNotFoundException, InvalidRecordException, SQLException, InvalidColourException{
		GenomicCoords gc= new GenomicCoords("2:1-1000", 80, null, null);
		TrackIntervalFeature vcf= new TrackIntervalFeature("test_data/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz", gc);
		List<IntervalFeature> linf = vcf.getIntervalFeatureList();

		GenotypeMatrix gm= new GenotypeMatrix();
		// gm.makeMatrix(linf, 80, null);
		String x= gm.printToScreen(true, linf, 80, null);
		assertTrue(x.isEmpty());
	}

	@Test
	public void canHandleVCFWithNoSamples() throws Exception{
		GenomicCoords gc= new GenomicCoords("1:1-10000000", 80, null, null);
		TrackIntervalFeature vcf= new TrackIntervalFeature("test_data/CHD.exon.2010_03.sites.vcf.gz", gc);
		List<IntervalFeature> linf = vcf.getIntervalFeatureList();
		assertTrue(linf.size() > 0); // Make sure we do have some features otherwise the test is meaningless.
		
		GenotypeMatrix gm= new GenotypeMatrix();
		String x= gm.printToScreen(true, linf, 80, null);
		assertTrue(x.isEmpty());
	}

	@Test
	public void sampleOrderIsTheSameAsInVcf() throws Exception{
		GenomicCoords gc= new GenomicCoords("1:199982-200052", 70, null, null);
		TrackIntervalFeature vcf= new TrackIntervalFeature("test_data/sample_order.vcf", gc);
		List<IntervalFeature> linf = vcf.getIntervalFeatureList();

		GenotypeMatrix gm= new GenotypeMatrix();
		String x= gm.printToScreen(true, linf, 80, vcf.getVcfHeader());

		assertTrue(x.split("\n")[0].startsWith("HG03"));
		assertTrue(x.split("\n")[1].startsWith("HG01"));
		assertTrue(x.split("\n")[2].startsWith("HG02"));
	}
	
}
