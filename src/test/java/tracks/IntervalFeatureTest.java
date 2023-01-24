package tracks;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.apache.commons.lang3.math.NumberUtils;
import org.junit.Before;
import org.junit.Test;

import coloring.Config;
import exceptions.InvalidColourException;
import exceptions.InvalidConfigException;
import exceptions.InvalidGenomicCoordsException;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import samTextViewer.Utils;

public class IntervalFeatureTest {

	private String ideogramToString(List<FeatureChar> fchars, boolean noFormat) throws InvalidColourException{
		StringBuilder sb= new StringBuilder();
		for(FeatureChar x : fchars){
			sb.append(x.format(noFormat));
		}
		return sb.toString();
	}
	
	@Before
	public void setConfig() throws IOException, InvalidConfigException{
		new Config(null);
	}
	
	@Test
	public void beahviourOfIsNumber() {

		// Valid numbers
		assertTrue(NumberUtils.isCreatable("1.1"));
		assertTrue(NumberUtils.isCreatable("001"));
		assertTrue(NumberUtils.isCreatable("0x0004")); // Also valid
		assertTrue(NumberUtils.isCreatable("-1.1e9"));
		
		// Invalid numbers
		assertFalse(NumberUtils.isCreatable("001.1"));
		assertFalse(NumberUtils.isCreatable("1.1 ")); // Note trailing space
		assertFalse(NumberUtils.isCreatable(""));
	}

	@Test 
	public void canGetMidPointOfFeature() throws InvalidGenomicCoordsException{

		String line= "chr1 0 100".replaceAll(" ", "\t"); // Genomic coords are irrelavant to this test
		IntervalFeature f= new IntervalFeature(line, TrackFormat.BED, null, -1);
		f.setScreenFrom(0);
		f.setScreenTo(0);
		assertEquals(0, f.getScreenMid());
		
		f.setScreenFrom(10);
		f.setScreenTo(10);
		assertEquals(10, f.getScreenMid());
		
		f.setScreenFrom(10); // Even
		f.setScreenTo(13);
		assertEquals(11, f.getScreenMid());
		
		f.setScreenFrom(0); // Even
		f.setScreenTo(3);
		assertEquals(1, f.getScreenMid());
		
		f.setScreenFrom(0); // Odd
		f.setScreenTo(4);
		assertEquals(2, f.getScreenMid());
		
		f.setScreenFrom(10); // Odd
		f.setScreenTo(14);
		assertEquals(12, f.getScreenMid());
	}
	
//	@Test
//	public void longestRun() throws InvalidGenomicCoordsException{
//		
//		IntervalFeature plus= new IntervalFeature("chr1 0 10 x . +".replaceAll(" ", "\t"), TrackFormat.BED, null);
//		System.out.println(Arrays.toString(plus.maxRun("zzBBBccc".toCharArray())));
//		System.out.println(Arrays.toString(plus.maxRun("xAAAAAAzz".toCharArray())));
//		System.out.println(Arrays.toString(plus.maxRun("xxAAAAAAzz".toCharArray())));
//	}
	
	@Test
	public void canMakeIdeogram() throws InvalidGenomicCoordsException, InvalidColourException{
		
		String line= "chr1 0 10".replaceAll(" ", "\t");
		IntervalFeature f= new IntervalFeature(line, TrackFormat.BED, null, -1);
		f.setScreenFrom(0);
		f.setScreenTo(9);

		String ideogram= this.ideogramToString(f.getIdeogram(true, true), true);
		assertEquals("||||||||||", ideogram);

		f.setStrand('+');
		ideogram= this.ideogramToString(f.getIdeogram(true, true), true);
		assertEquals(">>>>>>>>>>", ideogram);
		
		// Can add name
		f.setName("foo");
		ideogram= this.ideogramToString(f.getIdeogram(true, true), true);
		assertTrue(ideogram.contains("foo"));
		assertTrue(f.getIdeogram(true, true).size() > 4);

		// With GTF feature
		line= "chr1 na exon 1 100 . + . ID=mrna0001;foo=myname".replaceAll(" ", "\t");
		f= new IntervalFeature(line, TrackFormat.GTF, null, -1);
		f.setScreenFrom(0);
		f.setScreenTo(9);
		f.setGtfAttributeForName("foo");
		ideogram= this.ideogramToString(f.getIdeogram(true, true), true);
		assertTrue(ideogram.contains("E_myname_E"));

		// No key value found
		f.setGtfAttributeForName("NOT");
		ideogram= this.ideogramToString(f.getIdeogram(true, true), true);
		assertEquals("EEEEEEEEEE", ideogram);

		// Use default
		f.setGtfAttributeForName(null);
		ideogram= this.ideogramToString(f.getIdeogram(true, true), true);
		assertTrue(ideogram.contains("mrna"));
		
	}

	@Test
	public void canSetIdeogram() throws InvalidGenomicCoordsException, InvalidColourException{

		String line= "chr1 0 10".replaceAll(" ", "\t");
		IntervalFeature f= new IntervalFeature(line, TrackFormat.BED, null, -1);
		f.setScreenFrom(0);
		f.setScreenTo(9);

		// Use an ideogram created from outside
		List<FeatureChar> xIdeogram= new ArrayList<FeatureChar>(); 
		for(int i= 0; i< "xxxxxxxxxx".length(); i++){
			FeatureChar x= new FeatureChar();
			x.setText('x');
			xIdeogram.add(x);
		}
		f.setName("foo");
		f.setIdeogram(xIdeogram, true);

		String ideogram = this.ideogramToString(f.getIdeogram(false, true), true);
		assertTrue(ideogram.contains("x_foo_x"));
		
		// Use default, internal ideogram
		f.setIdeogram(null, false);
		ideogram= this.ideogramToString(f.getIdeogram(true, true), true);
		assertTrue(ideogram.contains("|_foo_|"));
	
	}
	
	@Test
	public void canSetFeatureNameFromBedField() throws InvalidGenomicCoordsException, InvalidColourException{

		String line= "chr1 0 10 foo bar baz".replaceAll(" ", "\t");
		IntervalFeature f= new IntervalFeature(line, TrackFormat.BED, null, -1);
		f.setScreenFrom(0);
		f.setScreenTo(9);
		
		// Default: 4th field for name
		String ideogram= this.ideogramToString(f.getIdeogram(true, true), true);
		assertTrue(ideogram.contains("foo"));
		
		f.setBedFieldName(4);
		ideogram= this.ideogramToString(f.getIdeogram(true, true), true);
		assertTrue(ideogram.contains("bar"));
		
		f.setBedFieldName(5);
		ideogram= this.ideogramToString(f.getIdeogram(true, true), true);
		assertTrue(ideogram.contains("baz"));
		
		// Invalid index: No change
		f.setBedFieldName(99);
		ideogram= this.ideogramToString(f.getIdeogram(true, true), true);
		assertTrue(ideogram.contains("baz"));
		
		// Column index not available
		line= "chr1 0 10".replaceAll(" ", "\t");
		f= new IntervalFeature(line, TrackFormat.BED, null, -1);
		f.setScreenFrom(0);
		f.setScreenTo(9);
		f.setBedFieldName(4);
		ideogram= this.ideogramToString(f.getIdeogram(true, true), true);
		assertTrue(ideogram.replaceAll("\\|", "").isEmpty());
		
		// Name missing
		line= "chr1 0 10 .".replaceAll(" ", "\t");
		f= new IntervalFeature(line, TrackFormat.BED, null, -1);
		f.setScreenFrom(0);
		f.setScreenTo(9);
		f.setBedFieldName(3);
		ideogram= this.ideogramToString(f.getIdeogram(true, true), true);
		assertTrue(ideogram.replaceAll("\\|", "").isEmpty());
		
		// Name not wanted
		line= "chr1 0 10 foo".replaceAll(" ", "\t");
		f= new IntervalFeature(line, TrackFormat.BED, null, -1);
		f.setScreenFrom(0);
		f.setScreenTo(9);
		f.setBedFieldName(3);
		ideogram= this.ideogramToString(f.getIdeogram(true, true), true);
		assertTrue(ideogram.contains("foo"));
		f.setBedFieldName(-1);
		ideogram= this.ideogramToString(f.getIdeogram(true, true), true);
		assertTrue( ! ideogram.contains("foo"));
	}
	
	@Test
	public void canAssignNameToFeature() throws InvalidGenomicCoordsException{
		
		String line; 
		IntervalFeature f;

		// Name not wanted to display:
		line= "chr1 0 10 myname".replaceAll(" ", "\t");
		f= new IntervalFeature(line, TrackFormat.BED, null, -1);
		f.setGtfAttributeForName("-na");
		assertEquals(".", f.getName());
		
		// Do not use a name
		line= "chr1 0 10 myname".replaceAll(" ", "\t");
		f= new IntervalFeature(line, TrackFormat.BED, null, -1);
		f.setGtfAttributeForName(null);
		assertEquals("myname", f.getName());
		
		//Custom name from BED: Has no effect
		line= "chr1 0 10 myname".replaceAll(" ", "\t");
		f= new IntervalFeature(line, TrackFormat.BED, null, -1);
		f.setGtfAttributeForName("ID");
		assertEquals("myname", f.getName());
		
		//Custom name from GTF
		line= "chr1 na exon 1 10 . + . ID=mrna0001;foo=myname".replaceAll(" ", "\t");
		f= new IntervalFeature(line, TrackFormat.GTF, null, -1);
		f.setGtfAttributeForName("foo");
		assertEquals("myname", f.getName());
		
		//Custom name from GTF, with attribute not found
		line= "chr1 na exon 1 10 . + . ID=mrna0001;Name=myname".replaceAll(" ", "\t");
		f= new IntervalFeature(line, TrackFormat.GTF, null, -1);
		f.setGtfAttributeForName("foo");
		assertEquals(".", f.getName());
		
		// BED with and without name
		line= "chr1 0 10".replaceAll(" ", "\t");
		f= new IntervalFeature(line, TrackFormat.BED, null, -1);
		assertEquals(".", f.getName());
		
		line= "chr1 0 10 myname".replaceAll(" ", "\t");
		f= new IntervalFeature(line, TrackFormat.BED, null, -1);
		assertEquals("myname", f.getName());
		
		// GTF/GFF without attributes or without any valid filed to get name from
		line= "chr1 na exon 1 10 . + .".replaceAll(" ", "\t");
		f= new IntervalFeature(line, TrackFormat.GTF, null, -1);
		assertEquals(".", f.getName());
		
		
		line= "chr1 na exon 1 10 . + . foo=mrna0001;bar=myname".replaceAll(" ", "\t");
		f= new IntervalFeature(line, TrackFormat.GTF, null, -1);
		assertEquals(".", f.getName());
		
		// GFF with Name
		line= "chr1 na exon 1 10 . + . ID=mrna0001;Name=myname".replaceAll(" ", "\t");
		f= new IntervalFeature(line, TrackFormat.GTF, null, -1);
		assertEquals("myname", f.getName());

		// GFF use ID
		line= "chr1 na exon 1 10 . + . ID=mrna0001;foo=myname".replaceAll(" ", "\t");;
		f= new IntervalFeature(line, TrackFormat.GTF, null, -1);
		assertEquals("mrna0001", f.getName());
		
		//GTF
		line= "chr1\tna\texon\t1\t10\t.\t+\t.\tgene_id \"mygene\"; transcript_id \"mytranx\";";
		f= new IntervalFeature(line, TrackFormat.GTF, null, -1);
		assertEquals("mytranx", f.getName());

		line= "chr1\tna\texon\t1\t10\t.\t+\t.\tgene_id \"mygene\";";
		f= new IntervalFeature(line, TrackFormat.GTF, null, -1);
		assertEquals("mygene", f.getName());
		
	}
	
	@Test
	public void canTestForEqualCoords() throws InvalidGenomicCoordsException{
		IntervalFeature plus= new IntervalFeature("chr1 0 10 x . +".replaceAll(" ", "\t"), TrackFormat.BED, null, -1);
		IntervalFeature minus= new IntervalFeature("chr1 0 10 y . -".replaceAll(" ", "\t"), TrackFormat.BED, null, -1);
		
		assertTrue(plus.equals(minus)); // Strand not matters
		assertTrue(! plus.equalStranded(minus)); // Strand matters
		
		// Strand NA
		IntervalFeature na1= new IntervalFeature("chr1 0 10".replaceAll(" ", "\t"), TrackFormat.BED, null, -1);
		IntervalFeature na2= new IntervalFeature("chr1 0 10".replaceAll(" ", "\t"), TrackFormat.BED, null, -1);
		assertTrue(na1.equals(na2));
		assertTrue(na1.equalStranded(na2));
		
		assertTrue(plus.equals(na2));
		assertTrue(!plus.equalStranded(na2));
	}
	
	@Test
	public void canCreateIntervalFromGtfString() throws InvalidGenomicCoordsException{
		String gtfLine= "chr1\tunknown\texon\t11874\t12227\t.\t+\t.\tgene_id \"DDX11L1\"; transcript_id \"NR_046018_1\"; gene_name \"DDX11L1\"; tss_id \"TSS14523\";";
		IntervalFeature f= new IntervalFeature(gtfLine, TrackFormat.GTF, null, -1);
		assertEquals(11874, f.getFrom());
		assertEquals(12227, f.getTo());
		assertEquals("exon", f.getFeature());

	}
	
	@Test
	public void canGetAttribute() throws InvalidGenomicCoordsException{
		
		/* Enable this test by uncommenting IF you actually need to use getAttribute() */
		
		String gtfLine= "chr1\tunknown\texon\t11874\t12227\t.\t+\t.\tgene_id \"DDX11L1\"; transcript_id \"NR_046018_1\"; gene_name \"DDX11L1\"; tss_id \"TSS14523\";";	
		
		IntervalFeature f= new IntervalFeature(gtfLine, TrackFormat.GTF, null, -1);
		String x= f.getGFFValueFromKey("transcript_id");
		assertEquals("NR_046018_1", x);
	}

	
	@Test
	public void canCreateIntervalFromString() throws InvalidGenomicCoordsException{
		String bedLine= "chr1\t0\t1";
		IntervalFeature f= new IntervalFeature(bedLine, TrackFormat.BED, null, -1);
		assertEquals("chr1", f.getChrom());
		assertEquals(1, f.getFrom()); // Note start augmented by 1.
		
		bedLine= "chr1\t0\t1\tgene\t0.1\t+";
		f= new IntervalFeature(bedLine, TrackFormat.BED, null, -1);
		assertEquals("gene", f.getName());
		assertEquals('+', f.getStrand());

		bedLine= " chr1\t0\t1";
		f= new IntervalFeature(bedLine, TrackFormat.BED, null, -1); 
		assertEquals("chr1", f.getChrom()); // NB: spaces in chrom stripped.
	}
	
	@Test
	public void canSortByChromPos() throws InvalidGenomicCoordsException{
		List<IntervalFeature> flist= new ArrayList<IntervalFeature>();
		flist.add(new IntervalFeature("chrM\t10\t100\tg5", TrackFormat.BED, null, -1));
		flist.add(new IntervalFeature("chrM\t1\t100\tg4", TrackFormat.BED, null, -1));
		flist.add(new IntervalFeature("chrM\t1\t90\tg3", TrackFormat.BED, null, -1));
		flist.add(new IntervalFeature("chr1\t10\t90\tg2", TrackFormat.BED, null, -1));
		flist.add(new IntervalFeature("chr1\t1\t90\tg1", TrackFormat.BED, null, -1));
		Collections.sort(flist);
		assertEquals("g1", flist.get(0).getName());
		assertEquals("g5", flist.get(4).getName());
	}
	
	@Test
	public void canMapIntervalToRuler() throws InvalidGenomicCoordsException{
		List<Double> rulerMap= new ArrayList<Double>();
		for(int i= 10; i < 20; i += 2){
			rulerMap.add((double)(i + 0.3));
		} // [10.3, 12.3, 14.3, 16.3, 18.3]
		IntervalFeature f= new IntervalFeature("chrM\t10\t15", TrackFormat.BED, null, -1);
		f.mapToScreen(rulerMap);
		assertEquals(0, f.getScreenFrom());
		assertEquals(2, f.getScreenTo());

		// Feature is fully contained in just one text char on screen
		f= new IntervalFeature("chrM\t10\t11", TrackFormat.BED, null, -1);
		f.mapToScreen(rulerMap);
		assertEquals(0, f.getScreenFrom());
		assertEquals(0, f.getScreenTo());
		
		// Feature is not part of ruler:
		f= new IntervalFeature("chrM\t100\t500", TrackFormat.BED, null, -1);
		f.mapToScreen(rulerMap);
		assertEquals(-1, f.getScreenFrom());
		assertEquals(-1, f.getScreenTo());

		// NB: Feature is not part of ruler because start coord 18 is augmented by 1 to become 1-based.
		f= new IntervalFeature("chrM\t18\t30", TrackFormat.BED, null, -1);
		f.mapToScreen(rulerMap);
		assertEquals(-1, f.getScreenFrom());
		assertEquals(-1, f.getScreenTo());
		
		// Partial overlap:
		f= new IntervalFeature("chrM\t1\t16", TrackFormat.BED, null, -1);
		f.mapToScreen(rulerMap);
		assertEquals(0, f.getScreenFrom());
		assertEquals(3, f.getScreenTo());

		f= new IntervalFeature("chrM\t17\t30", TrackFormat.BED, null, -1);
		f.mapToScreen(rulerMap);
		assertEquals(4, f.getScreenFrom());
		assertEquals(4, f.getScreenTo());
	}
	
	@Test
	public void handleSpaceInVCFInfoAndInvalidKey(){
		// This reading would fail on the original htsjdk-1.141 
		VCFFileReader reader = new VCFFileReader(new File("test_data/malformed.vcf.gz"));
		reader.query("chr1", 1, 16000000);
		reader.close();
	}

	@Test
	public void handleMalformedHeader(){
		// This reading would fail on the original htsjdk-1.141 
		VCFFileReader reader = new VCFFileReader(new File("test_data/malformed_header.vcf.gz"));
		reader.close();
	}
	
	@Test
	public void handleUnsupportedVersion(){
		// This reading would fail on the original htsjdk-1.141 
		VCFFileReader reader = new VCFFileReader(new File("test_data/malformed_header2.vcf.gz"));
		reader.close();
	}

	@Test
	public void handleMissingVersionLine(){
		VCFFileReader reader = new VCFFileReader(new File("test_data/malformed_header4.vcf.gz"));
		reader.close();
	}
	
	@Test
	public void handleMissingInfoLinesInHeader(){
		// This is fine also with original htsjdk-1.141 
		VCFFileReader reader = new VCFFileReader(new File("test_data/malformed_header3.vcf.gz"));
		reader.close();
	}
	
	@Test
	public void canFormatVCFLine() throws InvalidGenomicCoordsException, InvalidColourException, IOException, InvalidConfigException{
		
		List<Double> rulerMap= new ArrayList<Double>();
		for(int i= 1; i < 100; i++){
			rulerMap.add((double)i);
		}
		
		// Prepare header
		VCFFileReader reader = new VCFFileReader(new File("test_data/CHD.exon.2010_03.sites.vcf.gz"));
		VCFHeader vcfHeader= reader.getFileHeader();
		reader.close();
		VCFCodec vcfCodec= new VCFCodec();
		vcfCodec.setVCFHeader(vcfHeader, Utils.getVCFHeaderVersion(vcfHeader));
		
		String vcfLine= "1 10 . C G 23 PASS AA=.,foo;AC=.;AN=.DP=.".replaceAll(" ", "\t");
		IntervalFeature ift= new IntervalFeature(vcfLine, TrackFormat.VCF, vcfCodec, -1);
		
		ift.mapToScreen(rulerMap);
		assertEquals(1, ift.getIdeogram(true, true).size());
		assertEquals('G', ift.getIdeogram(true, true).get(0).getText());
		assertTrue(ift.getIdeogram(true, true).get(0).format(false).contains("[")); // Just check there is a formatting char
		
		// Deletion
		vcfLine= "1 10 . CTTG C 23 PASS AA=.;AC=.;AN=.DP=.".replaceAll(" ", "\t");
		ift= new IntervalFeature(vcfLine, TrackFormat.VCF, vcfCodec, -1);
		ift.mapToScreen(rulerMap);
		assertEquals('D', ift.getIdeogram(true, true).get(0).getText());
		assertEquals(3, ift.getIdeogram(true, true).size());
		
		// Insertion
		vcfLine= "1 10 . C CTTG 23 PASS AA=.;AC=.;AN=.DP=.".replaceAll(" ", "\t");
		ift= new IntervalFeature(vcfLine, TrackFormat.VCF, vcfCodec, -1);
		ift.mapToScreen(rulerMap);
		assertEquals("I", ift.getIdeogram(true, true).get(0).format(true));
		assertEquals(3, ift.getIdeogram(true, true).size());
		
		// Multiple alleles
		vcfLine= "1 10 . C CTTG,A 23 PASS AA=.;AC=.;AN=.DP=.".replaceAll(" ", "\t");
		ift= new IntervalFeature(vcfLine, TrackFormat.VCF, vcfCodec, -1);
		ift.mapToScreen(rulerMap);
		assertEquals("|", ift.getIdeogram(true, true).get(0).format(true));
	}
	
	@Test
	public void canFormatVCFLineStructVar() throws InvalidGenomicCoordsException, InvalidColourException, IOException, InvalidConfigException{
		
		List<Double> rulerMap= new ArrayList<Double>();
		for(int i= 1; i < 100; i++){
			rulerMap.add((double)i);
		}
		
		// Prepare header
		VCFFileReader reader = new VCFFileReader(new File("test_data/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz"));
		VCFHeader vcfHeader= reader.getFileHeader();
		reader.close();
		VCFCodec vcfCodec= new VCFCodec();
		vcfCodec.setVCFHeader(vcfHeader, Utils.getVCFHeaderVersion(vcfHeader));

		String vcfLine= "1 668630 DUP_delly_DUP20532 G <CN2> . PASS AC=64;AF=0.0127795;AFR_AF=0.0015;AMR_AF=0;AN=5008;CIEND=-150,150;CIPOS=-150,150;CS=DUP_delly;EAS_AF=0.0595;END=850204;EUR_AF=0.001;IMPRECISE;NS=2504;SAS_AF=0.001;SITEPOST=1;SVTYPE=DUP GT 0|0 0|0 0|0".replaceAll(" ", "\t");
		IntervalFeature ift= new IntervalFeature(vcfLine, TrackFormat.VCF, vcfCodec, -1);
		ift.mapToScreen(rulerMap);
		assertEquals(850204, ift.getTo());
		assertEquals("|", ift.getIdeogram(true, true).get(0).format(true));
	}
	
	@Test
	public void canParseScorecolumn() throws InvalidGenomicCoordsException {
	    
	    String line= "chr1 0 1 9 8 FOO 10".replaceAll(" ", "\t");

	    // Default column indexes for scores
	    IntervalFeature ift= new IntervalFeature(line, TrackFormat.BED, null, -1);
        assertEquals(8, ift.getScore(), 0.001);
        ift= new IntervalFeature(line, TrackFormat.BEDGRAPH, null, -1);
        assertEquals(9, ift.getScore(), 0.001);
        
        ift= new IntervalFeature(line, TrackFormat.BEDGRAPH, null, 7);
        assertEquals(10, ift.getScore(), 0.001);
        ift= new IntervalFeature(line, TrackFormat.BED, null, 7);
        assertEquals(10, ift.getScore(), 0.001);
	
        ift= new IntervalFeature(line, TrackFormat.BED, null, 6);
        assertEquals(Float.NaN, ift.getScore(), 0.001);
        ift= new IntervalFeature(line, TrackFormat.BEDGRAPH, null, 6);
        assertEquals(Float.NaN, ift.getScore(), 0.001);
        
        ift= new IntervalFeature(line, TrackFormat.BED, null, 20);
        assertEquals(Float.NaN, ift.getScore(), 0.001);
        ift= new IntervalFeature(line, TrackFormat.BEDGRAPH, null, 20);
        assertEquals(Float.NaN, ift.getScore(), 0.001);
        
        line= "chr1 0 1".replaceAll(" ", "\t");
	}
	
//	@Test
//	public void canPrintNameOnFeature() throws InvalidGenomicCoordsException{
//
//		String bedLine= "chr1 0 50 ACTB".replaceAll(" ", "\t");
//		IntervalFeature ift= new IntervalFeature(bedLine, TrackFormat.BED, null);
//		ift.setScreenFrom(0);
//		ift.setScreenTo(4);
//		System.out.println(Arrays.toString(ift.makeIdeogramFormatted(true)));
//
//		
//		String gtfLine= "36 GeneDB exon 4809 7349 . - . gene_id".replaceAll(" ", "\t");
//		IntervalFeature ift= new IntervalFeature(gtfLine, TrackFormat.GTF, null);
//		assertEquals("e", ift.getFormattedTextForScreen(true)[0]);		
//		assertTrue(ift.getFormattedTextForScreen(false)[0].startsWith("["));
//	}
	
}
