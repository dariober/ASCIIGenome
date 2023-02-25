package tracks;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.lang3.StringUtils;
import org.junit.Before;
import org.junit.Test;

import com.google.common.base.Splitter;

import coloring.Config;
import coloring.ConfigKey;
import exceptions.InvalidColourException;
import exceptions.InvalidConfigException;
import exceptions.InvalidGenomicCoordsException;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import samTextViewer.GenomicCoords;

public class TextReadTest {

	static SamReaderFactory srf=SamReaderFactory.make();
	static SamReader samReader= srf.open(new File("test_data/ds051.short.bam"));
	public static SAMSequenceDictionary samSeqDict= samReader.getFileHeader().getSequenceDictionary();
	
	public static String fastaFile= "test_data/chr7.fa";
	
	@Before
	public void setConfig() throws IOException, InvalidConfigException{
		new Config(null);
	}

    @Test
    public void canResetColorForStructuralVariant() throws InvalidGenomicCoordsException, IOException, InvalidColourException {
        GenomicCoords gc= new GenomicCoords("chr7:1-80", 80, null, null);
        SAMRecord rec= new SAMRecord(null);
        rec.setAlignmentStart(1);
        rec.setCigarString("1M");
        rec.setMappingQuality(30);
        rec.setReadBases("A".getBytes());
        rec.setAttribute("SA", "foo");
        TextRead tr= new TextRead(rec, gc, false);
        String txt = tr.getPrintableTextRead(false, false, false);
        
        // Default shading
        String shade = Config.get(ConfigKey.shade_structural_variant);
        String bg = Config.get(ConfigKey.background);
        assertTrue(txt.contains(";" + shade + ";"));
        
        // Omit shading
        Config.set(ConfigKey.shade_structural_variant, "false");
        txt = tr.getPrintableTextRead(false, false, false);
        assertTrue(!txt.contains(";" + shade + ";"));
        assertTrue(txt.contains(";" + bg + ";"));
        
        // Another shading
        Config.set(ConfigKey.shade_structural_variant, "123");
        txt = tr.getPrintableTextRead(false, false, false);
        assertTrue(txt.contains(";" + "123" + ";"));
        System.err.println(txt);
    }
	
	@Test
	public void getNumAlignedBasesFromCigar() {
		// A reminder of how cigar works.
		SAMRecord rec= new SAMRecord(null);
		rec.setCigarString("2M2D4M");
		int noRefBases= rec.getCigar().getReferenceLength();
		assertEquals(8, noRefBases);
		
		rec.setCigarString("2M2N4M"); // N same as D
		noRefBases= rec.getCigar().getReferenceLength();
		assertEquals(8, noRefBases);
		
		rec.setCigarString("2M2I4M");
		noRefBases= rec.getCigar().getReferenceLength();
		assertEquals(6, noRefBases);
		
		rec.setCigarString("2S10M4H");
		noRefBases= rec.getCigar().getReferenceLength();
		assertEquals(10, noRefBases);		
		
		// Cigar oprators consuming:
		rec.setCigarString("2H"); // Hard clipping doesn't consume bases on read or ref
		CigarElement el= rec.getCigar().getCigarElement(0);
		assertFalse(el.getOperator().consumesReadBases());		
		assertFalse(el.getOperator().consumesReferenceBases());

		rec.setCigarString("2S"); // Clips the read not the ref.
		el= rec.getCigar().getCigarElement(0);
		assertTrue(el.getOperator().consumesReadBases());		
		assertFalse(el.getOperator().consumesReferenceBases());

		rec.setCigarString("2D"); // This makes a gap in the read when aligned to ref
		el= rec.getCigar().getCigarElement(0);
		assertFalse(el.getOperator().consumesReadBases());		
		assertTrue(el.getOperator().consumesReferenceBases());

		rec.setCigarString("2I"); // Opposite of D
		el= rec.getCigar().getCigarElement(0);
		assertTrue(el.getOperator().consumesReadBases());		
		assertFalse(el.getOperator().consumesReferenceBases());
	}

	@Test
	public void canShowInsertion() throws InvalidGenomicCoordsException, IOException, InvalidColourException {
		
		GenomicCoords gc= new GenomicCoords("chr7:1-80", 80, null, null);
		SAMRecord rec= new SAMRecord(null);
		rec.setAlignmentStart(1);
		rec.setCigarString("1M3I5M1I3M");
		rec.setReadBases("AnnnACCGGnACT".getBytes());
		
		TextRead tr= new TextRead(rec, gc, false);
		// 7: code for "invert" colours. Checking for "[7;" is brittle since it assumes 7 is the first code.
		assertEquals(2, StringUtils.countMatches(tr.getPrintableTextRead(false, false, false), "[7;")); 

		// Ensure it is the first base being marked. I.e. the "A".
		assertTrue(tr.getPrintableTextRead(false, false, false).replaceAll("A.*", "").contains("[7;"));

		// Hide insertion if resolution is not single base
		GenomicCoords gc2= new GenomicCoords("chr7:1-80", 79, null, null);
		tr= new TextRead(rec, gc2, false);
		assertTrue( ! tr.getPrintableTextRead(false, false, false).contains("[7;"));
		
		rec.setCigarString("4M");
		rec.setReadBases("ACTG".getBytes());
		tr= new TextRead(rec, gc, false);
		assertTrue(! tr.getPrintableTextRead(false, false, false).contains("7"));
		
		rec.setCigarString("1M");
		rec.setReadBases("A".getBytes());
		tr= new TextRead(rec, gc, false);
	}	
	
	@Test
	public void canShadeBaseQuality() throws InvalidGenomicCoordsException, IOException, InvalidColourException {
		GenomicCoords gc= new GenomicCoords("chr7:1-80", 80, null, null);
		SAMRecord rec= new SAMRecord(null);
		rec.setAlignmentStart(1);
		rec.setCigarString("2M3D8M");
		rec.setMappingQuality(30);
		rec.setReadBases("AAAAATTTTT".getBytes());
		rec.setBaseQualities("!!!!!IIIII".getBytes());
		System.err.println(rec.getSAMString());
		TextRead tr= new TextRead(rec, gc, false);
		System.err.println(Splitter.on("m").omitEmptyStrings().splitToList(tr.getPrintableTextRead(false, false, false)));
	}

	@Test
	public void canShadeBaseQualityWithSoftClip() throws InvalidGenomicCoordsException, IOException, InvalidColourException {
		// Ensure this is the colour we expect with low base quality
		assertEquals("249", Config.get(ConfigKey.shade_low_mapq)); 
		
		GenomicCoords gc= new GenomicCoords("chr7:100-170", 80, null, null);
		SAMRecord rec= new SAMRecord(null);
		rec.setAlignmentStart(100);
		rec.setCigarString("2S4M");
		rec.setMappingQuality(30);
		rec.setReadBases("NNNNNN".getBytes());
		rec.setBaseQualityString("##ABC#");
		TextRead tr= new TextRead(rec, gc, false);
		List<String> printable= Splitter.on("N").omitEmptyStrings().splitToList(tr.getPrintableTextRead(false, false, false));
		assertTrue( ! printable.get(0).contains("249")); // Base qual= A
		assertTrue( ! printable.get(1).contains("249")); // Base qual= B
		assertTrue( ! printable.get(2).contains("249")); // Base qual= C
		assertTrue( printable.get(3).contains("249")); // Base qual= #
	}
	
//	@Test
//	public void canShowSoftClips() throws InvalidGenomicCoordsException, IOException, InvalidColourException {
//		
//		GenomicCoords gc= new GenomicCoords("chr7:1-80", 80, null, null);
//		SAMRecord rec= new SAMRecord(null);
//		rec.setAlignmentStart(1);
//		rec.setCigarString("1S10M5S3H");
//
//		TextRead textRead= new TextRead(rec, gc, true);
//
//		System.err.println(textRead.getSoftUnclippedAlignmentStart(textRead.getSamRecord()));
//		
//		assertEquals(15, textRead.getTextEnd());
//		assertEquals(1, textRead.getTextStart());
//		// We start with N even if the first operation is S. This is because the 
//		// alignment starts at 1 so the first soft clip base(s) is lost on the left. 
//		assertTrue(textRead.getPrintableTextRead(false, true, false).startsWith("N"));
//		assertTrue(textRead.getPrintableTextRead(false, true, false).endsWith("S"));
//		
//		rec.setCigarString("10S");
//		textRead= new TextRead(rec, gc, true);
//		assertEquals(10, textRead.getTextEnd());
//		assertEquals(1, textRead.getTextStart());
//		
//		rec.setCigarString("1S10M");
//		rec.setAlignmentStart(2);
//		textRead= new TextRead(rec, gc, true);
//		assertEquals(1, textRead.getTextStart());
//		
//		rec.setAlignmentStart(3);
//		textRead= new TextRead(rec, gc, true);
//		assertEquals(2, textRead.getTextStart());
//		
//		rec.setCigarString("1S2S10M");
//		rec.setAlignmentStart(4);
//		textRead= new TextRead(rec, gc, true);
//		assertEquals(1, textRead.getTextStart());
//	}
	
	@Test
	public void canPrintDNARead() throws InvalidGenomicCoordsException, IOException, InvalidColourException {
		
		GenomicCoords gc= new GenomicCoords("chr7:5566778-5566798", 80, samSeqDict, fastaFile);
		SAMRecord rec= new SAMRecord(null);
		rec.setAlignmentStart(5566780);
		rec.setCigarString("24M");
		rec.setReadNegativeStrandFlag(true);
		rec.setReadBases("AACCGGTTAACCGGTTAACCGGTT".getBytes());
		
		TextRead textRead= new TextRead(rec, gc, false);
		assertEquals(textRead.getTextStart(), 3);
		assertEquals(textRead.getTextEnd(), 21);
		assertEquals("AACCGGTTAACCGGTTAAC".length(), textRead.getTextEnd() - textRead.getTextStart() + 1);
		
		System.out.println(textRead.getPrintableTextRead(false, true, false));
		
		assertEquals("a,ccgg,t,acc,gtt,ac", textRead.getPrintableTextRead(false, true, false));
		assertEquals("a,ccgg,t,uccmgtt,ac", textRead.getPrintableTextRead(true, true, false));


		gc= new GenomicCoords("chr7:5566778-5566798", 80, samSeqDict, null);
		textRead= new TextRead(rec, gc, false);
		assertEquals("aaccggttaaccggttaac", textRead.getPrintableTextRead(false, true, false));

		gc= new GenomicCoords("chr7:5566780-5566782", 80, samSeqDict, fastaFile);
		textRead= new TextRead(rec, gc, false);
		System.out.println(textRead.getPrintableTextRead(true, true, false));
		System.out.println(textRead.getPrintableTextRead(true, false, false));		
	}

	@Test
	public void canPrintFormattedRead() throws InvalidGenomicCoordsException, IOException, InvalidColourException {
		
		GenomicCoords gc= new GenomicCoords("chr7:5566778-5566798", 80, null, null);
		SAMRecord rec= new SAMRecord(null);
		rec.setAlignmentStart(5566780);
		rec.setCigarString("24M");
		rec.setReadBases("AACCGGTTAACCGGTTAACCGGTT".getBytes());
		TextRead textRead= new TextRead(rec, gc, false);
		
		rec.setSecondOfPairFlag(false);
		assertTrue( ! textRead.getPrintableTextRead(true, false, false).contains("4;")); // '4': Underline

		rec.setSecondOfPairFlag(true);
		assertTrue( textRead.getPrintableTextRead(true, false, false).contains("4;")); // '4': Underline
	}
	
	@Test
	public void canPrintSquashedRead() throws InvalidGenomicCoordsException, IOException, InvalidColourException{
		GenomicCoords gc= new GenomicCoords("chr7:5566778-5566978", 80, samSeqDict, fastaFile);
		
		SAMRecord rec= new SAMRecord(null);
		rec.setAlignmentStart(5566780);
		rec.setCigarString("24M");
		rec.setReadBases("AACCGGTTAACCGGTTAACCGGTT".getBytes());
		
		TextRead textRead= new TextRead(rec, gc, false);
		assertEquals(2, textRead.getTextStart());
		assertEquals(11, textRead.getTextEnd());
		assertEquals(textRead.getTextEnd() - textRead.getTextStart() + 1, ">>>>>>>>>>".length());
		assertEquals(textRead.getPrintableTextRead(false, true, false), ">>>>>>>>>>");
		// System.out.println(textRead);
	}

	@Test
	public void canPrintReadWithSkippedBases() throws InvalidGenomicCoordsException, IOException, InvalidColourException{
		
		GenomicCoords gc= new GenomicCoords("chr7:1-800", 80, samSeqDict, fastaFile);
		
		SAMRecord rec= new SAMRecord(null);
		rec.setAlignmentStart(1);
		rec.setCigarString("70M200N50M200N10M");
		
		TextRead textRead= new TextRead(rec, gc, false);
		System.out.println(textRead.getPrintableTextRead(false, true, false));
		assertTrue(textRead.getPrintableTextRead(false, true, false).startsWith(">"));
		assertTrue(textRead.getPrintableTextRead(false, true, false).endsWith(">"));
		assertTrue(textRead.getPrintableTextRead(false, true, false).contains("_"));

		// Very large skipped region: Only show the skipped part.
		rec.setCigarString("1M5000N1M");
		textRead= new TextRead(rec, gc, false);
		System.out.println(textRead.getPrintableTextRead(false, true, false));
		assertTrue("_", textRead.getPrintableTextRead(false, true, false).startsWith("_"));
		assertTrue("_", textRead.getPrintableTextRead(false, true, false).endsWith("_"));
		assertEquals("", textRead.getPrintableTextRead(false, true, false).replaceAll("_", ""));
	}

//	@Test
//	public void canPrintWithReadName() throws InvalidGenomicCoordsException, IOException, InvalidColourException{
//
//		GenomicCoords gc= new GenomicCoords("chr7:5566778-5566978", 80, samSeqDict, fastaFile);
//		
//		SAMRecord rec= new SAMRecord(null);
//		rec.setAlignmentStart(5566780);
//		rec.setCigarString("24M");
//		rec.setReadBases("AACCGGTTAACCGGTTAACCGGTT".getBytes());
//		rec.setReadName("Read1");
//		
//		TextRead textRead= new TextRead(rec, gc, false);
//		assertEquals("Read1/>>>>", textRead.getPrintableTextRead(false, true, true));
//		
//		rec.setReadName("VeryLongReadNameMoreThanReadSequence");
//		textRead= new TextRead(rec, gc, false);
//		assertEquals("VeryLongRe", textRead.getPrintableTextRead(false, true, true));
//		
//		rec.setReadName("VeryLongRead");
//		textRead= new TextRead(rec, gc, false);
//		assertEquals("VeryLongRe", textRead.getPrintableTextRead(false, true, true));
//	}
}
