package tracks;

import static org.junit.Assert.*;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import org.apache.commons.lang3.math.NumberUtils;
import org.junit.Test;

import exceptions.InvalidGenomicCoordsException;
import tracks.IntervalFeature;

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
	
	@Test
	public void canAssignNameToFeature() throws InvalidGenomicCoordsException{
		
		String line; 
		IntervalFeature f;

		//Custom name from BED: Has no effect
		line= "chr1 0 10 myname".replaceAll(" ", "\t");
		f= new IntervalFeature(line, TrackFormat.BED);
		f.setGtfAttributeForName("ID");
		assertEquals("myname", f.getName());
		
		//Custom name from GTF
		line= "chr1 na exon 1 10 . + . ID=mrna0001;foo=myname".replaceAll(" ", "\t");
		f= new IntervalFeature(line, TrackFormat.GFF);
		f.setGtfAttributeForName("foo");
		assertEquals("myname", f.getName());
		
		//Custom name from GTF, with attribute not found
		line= "chr1 na exon 1 10 . + . ID=mrna0001;Name=myname".replaceAll(" ", "\t");
		f= new IntervalFeature(line, TrackFormat.GFF);
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
		f= new IntervalFeature(line, TrackFormat.GFF);
		assertEquals(".", f.getName());
		
		
		line= "chr1 na exon 1 10 . + . foo=mrna0001;bar=myname".replaceAll(" ", "\t");
		f= new IntervalFeature(line, TrackFormat.GFF);
		assertEquals(".", f.getName());
		
		// GFF with Name
		line= "chr1 na exon 1 10 . + . ID=mrna0001;Name=myname".replaceAll(" ", "\t");
		f= new IntervalFeature(line, TrackFormat.GFF);
		assertEquals("myname", f.getName());

		// GFF use ID
		line= "chr1 na exon 1 10 . + . ID=mrna0001;foo=myname".replaceAll(" ", "\t");;
		f= new IntervalFeature(line, TrackFormat.GFF);
		assertEquals("mrna0001", f.getName());
		
		//GTF
		line= "chr1\tna\texon\t1\t10\t.\t+\t.\tgene_id \"mygene\"; transcript_id \"mytranx\";";
		f= new IntervalFeature(line, TrackFormat.GFF);
		assertEquals("mytranx", f.getName());

		line= "chr1\tna\texon\t1\t10\t.\t+\t.\tgene_id \"mygene\";";
		f= new IntervalFeature(line, TrackFormat.GFF);
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
		IntervalFeature f= new IntervalFeature(gtfLine, TrackFormat.GFF);
		assertEquals(11874, f.getFrom());
		assertEquals(12227, f.getTo());
		assertEquals("exon", f.getFeature());
		System.out.println(f.toString());
	}
	
	@Test
	public void canGetAttribute() throws InvalidGenomicCoordsException{
		
		/* Enable this test by uncommenting IF you actually need to use getAttribute() */
		
		String gtfLine= "chr1\tunknown\texon\t11874\t12227\t.\t+\t.\tgene_id \"DDX11L1\"; transcript_id \"NR_046018_1\"; gene_name \"DDX11L1\"; tss_id \"TSS14523\";";	
		
		IntervalFeature f= new IntervalFeature(gtfLine, TrackFormat.GFF);
		//String x= f.getAttribute("transcript_id");
		//assertEquals("NR_046018_1", x);
		//System.out.println(f.toString());
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

		System.out.println(f.toString());
		
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
	public void canFormatVCFLine() throws InvalidGenomicCoordsException{
		String vcfLine= "1 113054374 . C G 23 PASS AC=2;AN=4;DP=4718;NS=65 GT:VR:DP:FT".replaceAll(" ", "\t");
		IntervalFeature ift= new IntervalFeature(vcfLine, TrackFormat.VCF);
		assertEquals("G", ift.assignTextToFeature(true));		
		assertTrue(ift.assignTextToFeature(false).trim().startsWith("["));
		
		// Deletion
		vcfLine= "1 113054374 . CTTG C 23 PASS AC=2;AN=4;DP=4718;NS=65 GT:VR:DP:FT".replaceAll(" ", "\t");
		ift= new IntervalFeature(vcfLine, TrackFormat.VCF);
		assertEquals("D", ift.assignTextToFeature(true));
		
		// Insertion
		vcfLine= "1 113054374 . C CTTG 23 PASS AC=2;AN=4;DP=4718;NS=65 GT:VR:DP:FT".replaceAll(" ", "\t");
		ift= new IntervalFeature(vcfLine, TrackFormat.VCF);
		assertEquals("I", ift.assignTextToFeature(true));
	}
	
	@Test
	public void canFormatGTFLine() throws InvalidGenomicCoordsException{
		String gtfLine= "36 GeneDB exon 4809 7349 . - . gene_id".replaceAll(" ", "\t");
		IntervalFeature ift= new IntervalFeature(gtfLine, TrackFormat.GFF);
		assertEquals("e", ift.assignTextToFeature(true));		
		assertTrue(ift.assignTextToFeature(false).trim().startsWith("["));
	}
	
}
