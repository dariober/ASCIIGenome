package tracks;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.apache.commons.lang3.math.NumberUtils;
import org.biojava.nbio.genome.parsers.gff.Feature;
import org.biojava.nbio.genome.parsers.gff.Location;
import org.junit.Test;

import com.google.common.base.Joiner;

import coloring.Config;
import exceptions.InvalidColourException;
import exceptions.InvalidConfigException;
import exceptions.InvalidGenomicCoordsException;

public class IntervalFeatureTest {

	@Test
	public void beahviourOfIsNumber() {

		// Valid numbers
		assertTrue(NumberUtils.isNumber("1.1"));
		assertTrue(NumberUtils.isNumber("001"));
		assertTrue(NumberUtils.isNumber("0x0004")); // Also valid
		assertTrue(NumberUtils.isNumber("-1.1e9"));
		// Invalid numbers
		assertFalse(NumberUtils.isNumber("1.1 ")); // Note trailing space
		assertFalse(NumberUtils.isNumber(""));
		assertFalse(NumberUtils.isNumber("001.1"));
	}
	
//	@Test
//	public void longestRun() throws InvalidGenomicCoordsException{
//		
//		IntervalFeature plus= new IntervalFeature("chr1 0 10 x . +".replaceAll(" ", "\t"), TrackFormat.BED);
//		System.out.println(Arrays.toString(plus.maxRun("zzBBBccc".toCharArray())));
//		System.out.println(Arrays.toString(plus.maxRun("xAAAAAAzz".toCharArray())));
//		System.out.println(Arrays.toString(plus.maxRun("xxAAAAAAzz".toCharArray())));
//	}
	
	@Test
	public void canMakeIdeogram() throws InvalidGenomicCoordsException, InvalidColourException{
		
		String line= "chr1 0 10".replaceAll(" ", "\t");
		IntervalFeature f= new IntervalFeature(line, TrackFormat.BED);
		f.setScreenFrom(0);
		f.setScreenTo(9);

		assertEquals("||||||||||", Joiner.on("").join(f.makeIdeogramFormatted(true)));

		f.setStrand('+');
		assertEquals(">>>>>>>>>>", Joiner.on("").join(f.makeIdeogramFormatted(true)));
		
		// Can add name
		f.setName("foo");
		assertTrue(Joiner.on("").join(f.makeIdeogramFormatted(true)).contains("foo"));
		assertTrue(f.makeIdeogramFormatted(false).length > 4);

		// Set ideogram from outside and put name on it
		f.setIdeogram("xxxxxxxxxx".toCharArray());
		assertTrue(Joiner.on("").join(f.makeIdeogramFormatted(true)).contains("x_foo_x"));
		
		f.setIdeogram(null);
		assertTrue(Joiner.on("").join(f.makeIdeogramFormatted(true)).contains(">_foo_>"));

		// With GTF feature
		line= "chr1 na exon 1 100 . + . ID=mrna0001;foo=myname".replaceAll(" ", "\t");
		f= new IntervalFeature(line, TrackFormat.GTF);
		f.setScreenFrom(0);
		f.setScreenTo(9);
		f.setGtfAttributeForName("foo");
		assertTrue(Joiner.on("").join(f.makeIdeogramFormatted(true)).contains("E_myname_E"));

		// No key value found
		f.setGtfAttributeForName("NOT");
		assertEquals("EEEEEEEEEE", Joiner.on("").join(f.makeIdeogramFormatted(true)));

		// Use default
		f.setGtfAttributeForName(null);
		assertTrue(Joiner.on("").join(f.makeIdeogramFormatted(true)).contains("mrna"));
		
	}
	
	// @Test
	public void canAssignNameToFeature() throws InvalidGenomicCoordsException{
		
		Location location= new Location(10, 100);
		Feature gff= new Feature("chr1", "foo", "bar", location, (double) 0, 0, "ID=mirna");
		System.out.println(gff.getAttribute(null));		
		
		String line; 
		IntervalFeature f;

		// Name not wanted to display:
		line= "chr1 0 10 myname".replaceAll(" ", "\t");
		f= new IntervalFeature(line, TrackFormat.BED);
		f.setGtfAttributeForName("-na");
		assertEquals(".", f.getName());
		
		// Do not use a name
		line= "chr1 0 10 myname".replaceAll(" ", "\t");
		f= new IntervalFeature(line, TrackFormat.BED);
		f.setGtfAttributeForName(null);
		assertEquals("myname", f.getName());
		
		//Custom name from BED: Has no effect
		line= "chr1 0 10 myname".replaceAll(" ", "\t");
		f= new IntervalFeature(line, TrackFormat.BED);
		f.setGtfAttributeForName("ID");
		assertEquals("myname", f.getName());

		
		//Custom name from GTF
		line= "chr1 na exon 1 10 . + . ID=mrna0001;foo=myname".replaceAll(" ", "\t");
		f= new IntervalFeature(line, TrackFormat.GTF);
		f.setGtfAttributeForName("foo");
		assertEquals("myname", f.getName());
		
		//Custom name from GTF, with attribute not found
		line= "chr1 na exon 1 10 . + . ID=mrna0001;Name=myname".replaceAll(" ", "\t");
		f= new IntervalFeature(line, TrackFormat.GTF);
		f.setGtfAttributeForName("foo");
		assertEquals(".", f.getName());
		
		// BED with and without name
		line= "chr1 0 10".replaceAll(" ", "\t");
		f= new IntervalFeature(line, TrackFormat.BED);
		assertEquals(".", f.getName());
		
		line= "chr1 0 10 myname".replaceAll(" ", "\t");
		f= new IntervalFeature(line, TrackFormat.BED);
		assertEquals("myname", f.getName());
		
		// GTF/GFF without attributes or without any valid filed to get name from
		line= "chr1 na exon 1 10 . + .".replaceAll(" ", "\t");
		f= new IntervalFeature(line, TrackFormat.GTF);
		assertEquals(".", f.getName());
		
		
		line= "chr1 na exon 1 10 . + . foo=mrna0001;bar=myname".replaceAll(" ", "\t");
		f= new IntervalFeature(line, TrackFormat.GTF);
		assertEquals(".", f.getName());
		
		// GFF with Name
		line= "chr1 na exon 1 10 . + . ID=mrna0001;Name=myname".replaceAll(" ", "\t");
		f= new IntervalFeature(line, TrackFormat.GTF);
		assertEquals("myname", f.getName());

		// GFF use ID
		line= "chr1 na exon 1 10 . + . ID=mrna0001;foo=myname".replaceAll(" ", "\t");;
		f= new IntervalFeature(line, TrackFormat.GTF);
		assertEquals("mrna0001", f.getName());
		
		//GTF
		line= "chr1\tna\texon\t1\t10\t.\t+\t.\tgene_id \"mygene\"; transcript_id \"mytranx\";";
		f= new IntervalFeature(line, TrackFormat.GTF);
		assertEquals("mytranx", f.getName());

		line= "chr1\tna\texon\t1\t10\t.\t+\t.\tgene_id \"mygene\";";
		f= new IntervalFeature(line, TrackFormat.GTF);
		assertEquals("mygene", f.getName());
		
	}
	
	@Test
	public void canTestForEqualCoords() throws InvalidGenomicCoordsException{
		IntervalFeature plus= new IntervalFeature("chr1 0 10 x . +".replaceAll(" ", "\t"), TrackFormat.BED);
		IntervalFeature minus= new IntervalFeature("chr1 0 10 y . -".replaceAll(" ", "\t"), TrackFormat.BED);
		
		assertTrue(plus.equals(minus)); // Strand not matters
		assertTrue(! plus.equalStranded(minus)); // Strand matters
		
		// Strand NA
		IntervalFeature na1= new IntervalFeature("chr1 0 10".replaceAll(" ", "\t"), TrackFormat.BED);
		IntervalFeature na2= new IntervalFeature("chr1 0 10".replaceAll(" ", "\t"), TrackFormat.BED);
		assertTrue(na1.equals(na2));
		assertTrue(na1.equalStranded(na2));
		
		assertTrue(plus.equals(na2));
		assertTrue(!plus.equalStranded(na2));
	}
	
	@Test
	public void canCreateIntervalFromGtfString() throws InvalidGenomicCoordsException{
		String gtfLine= "chr1\tunknown\texon\t11874\t12227\t.\t+\t.\tgene_id \"DDX11L1\"; transcript_id \"NR_046018_1\"; gene_name \"DDX11L1\"; tss_id \"TSS14523\";";
		IntervalFeature f= new IntervalFeature(gtfLine, TrackFormat.GTF);
		assertEquals(11874, f.getFrom());
		assertEquals(12227, f.getTo());
		assertEquals("exon", f.getFeature());

	}
	
	@Test
	public void canGetAttribute() throws InvalidGenomicCoordsException{
		
		/* Enable this test by uncommenting IF you actually need to use getAttribute() */
		
		String gtfLine= "chr1\tunknown\texon\t11874\t12227\t.\t+\t.\tgene_id \"DDX11L1\"; transcript_id \"NR_046018_1\"; gene_name \"DDX11L1\"; tss_id \"TSS14523\";";	
		
		IntervalFeature f= new IntervalFeature(gtfLine, TrackFormat.GTF);
		String x= f.getGFFValueFromKey("transcript_id");
		assertEquals("NR_046018_1", x);
	}

	
	@Test
	public void canCreateIntervalFromString() throws InvalidGenomicCoordsException{
		String bedLine= "chr1\t0\t1";
		IntervalFeature f= new IntervalFeature(bedLine, TrackFormat.BED);
		assertEquals("chr1", f.getChrom());
		assertEquals(1, f.getFrom()); // Note start augmented by 1.
		
		bedLine= "chr1\t0\t1\tgene\t0.1\t+";
		f= new IntervalFeature(bedLine, TrackFormat.BED);
		assertEquals("gene", f.getName());
		assertEquals('+', f.getStrand());

		bedLine= " chr1\t0\t1";
		f= new IntervalFeature(bedLine, TrackFormat.BED); 
		assertEquals("chr1", f.getChrom()); // NB: spaces in chrom stripped.
	}
	
	@Test
	public void canSortByChromPos() throws InvalidGenomicCoordsException{
		List<IntervalFeature> flist= new ArrayList<IntervalFeature>();
		flist.add(new IntervalFeature("chrM\t10\t100\tg5", TrackFormat.BED));
		flist.add(new IntervalFeature("chrM\t1\t100\tg4", TrackFormat.BED));
		flist.add(new IntervalFeature("chrM\t1\t90\tg3", TrackFormat.BED));
		flist.add(new IntervalFeature("chr1\t10\t90\tg2", TrackFormat.BED));
		flist.add(new IntervalFeature("chr1\t1\t90\tg1", TrackFormat.BED));
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
		IntervalFeature f= new IntervalFeature("chrM\t10\t15", TrackFormat.BED);
		f.mapToScreen(rulerMap);
		assertEquals(0, f.getScreenFrom());
		assertEquals(2, f.getScreenTo());

		// Feature is fully contained in just one text char on screen
		f= new IntervalFeature("chrM\t10\t11", TrackFormat.BED);
		f.mapToScreen(rulerMap);
		assertEquals(0, f.getScreenFrom());
		assertEquals(0, f.getScreenTo());
		
		// Feature is not part of ruler:
		f= new IntervalFeature("chrM\t100\t500", TrackFormat.BED);
		f.mapToScreen(rulerMap);
		assertEquals(-1, f.getScreenFrom());
		assertEquals(-1, f.getScreenTo());

		// NB: Feature is not part of ruler because start coord 18 is augmented by 1 to become 1-based.
		f= new IntervalFeature("chrM\t18\t30", TrackFormat.BED);
		f.mapToScreen(rulerMap);
		assertEquals(-1, f.getScreenFrom());
		assertEquals(-1, f.getScreenTo());
		
		// Partial overlap:
		f= new IntervalFeature("chrM\t1\t16", TrackFormat.BED);
		f.mapToScreen(rulerMap);
		assertEquals(0, f.getScreenFrom());
		assertEquals(3, f.getScreenTo());

		f= new IntervalFeature("chrM\t17\t30", TrackFormat.BED);
		f.mapToScreen(rulerMap);
		assertEquals(4, f.getScreenFrom());
		assertEquals(4, f.getScreenTo());
		//System.out.println(f);
		//System.out.println(rulerMap);
	}
	
	@Test
	public void canFormatVCFLine() throws InvalidGenomicCoordsException, InvalidColourException, IOException, InvalidConfigException{
		
		new Config(null);
		
		List<Double> rulerMap= new ArrayList<Double>();
		for(int i= 1; i < 30; i++){
			rulerMap.add((double)i);
		}
		
		String vcfLine= "1 10 . C G 23 PASS AC=2;AN=4;DP=4718;NS=65 GT:VR:DP:FT".replaceAll(" ", "\t");
		IntervalFeature ift= new IntervalFeature(vcfLine, TrackFormat.VCF);
		ift.mapToScreen(rulerMap);
		assertEquals(1, ift.makeIdeogramFormatted(true).length);
		assertEquals("G", ift.makeIdeogramFormatted(true)[0]);
		assertTrue(ift.makeIdeogramFormatted(false)[0].contains("[")); // Just check there is a formatting char
		
		// Deletion
		vcfLine= "1 10 . CTTG C 23 PASS AC=2;AN=4;DP=4718;NS=65 GT:VR:DP:FT".replaceAll(" ", "\t");
		ift= new IntervalFeature(vcfLine, TrackFormat.VCF);
		ift.mapToScreen(rulerMap);
		assertEquals("D", ift.makeIdeogramFormatted(true)[0]);
		
		// Insertion
		vcfLine= "1 10 . C CTTG 23 PASS AC=2;AN=4;DP=4718;NS=65 GT:VR:DP:FT".replaceAll(" ", "\t");
		ift= new IntervalFeature(vcfLine, TrackFormat.VCF);
		ift.mapToScreen(rulerMap);
		assertEquals("I", ift.makeIdeogramFormatted(true)[0]);
	}
	
//	@Test
//	public void canPrintNameOnFeature() throws InvalidGenomicCoordsException{
//
//		String bedLine= "chr1 0 50 ACTB".replaceAll(" ", "\t");
//		IntervalFeature ift= new IntervalFeature(bedLine, TrackFormat.BED);
//		ift.setScreenFrom(0);
//		ift.setScreenTo(4);
//		System.out.println(Arrays.toString(ift.makeIdeogramFormatted(true)));
//
//		
//		String gtfLine= "36 GeneDB exon 4809 7349 . - . gene_id".replaceAll(" ", "\t");
//		IntervalFeature ift= new IntervalFeature(gtfLine, TrackFormat.GTF);
//		assertEquals("e", ift.getFormattedTextForScreen(true)[0]);		
//		assertTrue(ift.getFormattedTextForScreen(false)[0].startsWith("["));
//	}
	
}
