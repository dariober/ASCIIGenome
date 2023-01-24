package tracks;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.IOException;
import java.sql.SQLException;
import java.util.List;
import java.util.regex.Pattern;

import org.junit.Before;
import org.junit.Test;

import com.google.common.base.Splitter;

import coloring.Config;
import exceptions.InvalidColourException;
import exceptions.InvalidCommandLineException;
import exceptions.InvalidConfigException;
import exceptions.InvalidGenomicCoordsException;
import exceptions.InvalidRecordException;
import samTextViewer.GenomicCoords;
import samTextViewer.Utils;

public class TrackTest {

	@Before
	public void setConfig() throws IOException, InvalidConfigException{
		new Config(null);
	}
	
//	@Test
//	public void canSetVariantReadsFilterBasedOnScreenPercent() throws InvalidCommandLineException, InvalidGenomicCoordsException, IOException, ClassNotFoundException, InvalidRecordException, SQLException{
//		
//		GenomicCoords gc= new GenomicCoords("chr7:3-101", 1000, null, "test_data/chr7.fa");
//		
//		Track t1= new TrackReads("test_data/ds051.actb.bam", gc);
//
//		t1.setVariantReadInInterval("0.0", true);
//		assertEquals(3, t1.getFeatureFilter().getVariantFrom());
//		assertEquals(3, t1.getFeatureFilter().getVariantTo());
//		
//		t1.setVariantReadInInterval("1.0", true);
//		assertEquals(101, t1.getFeatureFilter().getVariantFrom());
//		assertEquals(101, t1.getFeatureFilter().getVariantTo());
//		
//		//     (101-3+1) * 0.5 = 49.5
//		// We round to 49.5 to 50. So: 
//		//     3 + 50 - 1 = 52
//		t1.setVariantReadInInterval("0.5", true);
//		assertEquals(52, t1.getFeatureFilter().getVariantFrom());
//		
//		gc= new GenomicCoords("chr7:3-102", 1000, null, "test_data/chr7.fa");
//		t1= new TrackReads("test_data/ds051.actb.bam", gc);
//		t1.setVariantReadInInterval("0.5", true);
//		assertEquals(52, t1.getFeatureFilter().getVariantFrom());
//		
//		// Must fail
//		boolean pass= false;
//		try{
//			t1.setVariantReadInInterval("foo", true);
//		} catch(NumberFormatException e){
//			pass= true;
//		}
//		assertTrue(pass);
//	}
	
	@Test
	public void canConcatTitleAndTrack() throws ClassNotFoundException, IOException, InvalidGenomicCoordsException, InvalidRecordException, SQLException, InvalidConfigException, InvalidColourException{
		new Config(null);
		GenomicCoords gc= new GenomicCoords("1:735171-2045891", 80, null, null);
		TrackIntervalFeature tif= new TrackIntervalFeature("test_data/CHD.exon.2010_03.sites.vcf", gc);
		tif.setNoFormat(true);
		tif.setTrackTag("title.bed");
		String[] lines= tif.concatTitleAndTrack().split("\n");
		assertEquals(3, lines.length);
		assertTrue(lines[0].startsWith("title.bed"));
		assertTrue(lines[0].trim().endsWith("A"));
		System.err.println(tif.concatTitleAndTrack());
		
		// Do not put on different lines
		tif.setTrackTag("LoooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooongTitle.bed");
		lines= tif.concatTitleAndTrack().split("\n");
		assertEquals(4, lines.length);
		System.err.println(tif.concatTitleAndTrack());
		
	}
	
	@Test
	public void canConcatTitleAndTrackWithNoFeatures() throws ClassNotFoundException, IOException, InvalidGenomicCoordsException, InvalidRecordException, SQLException, InvalidConfigException, InvalidColourException{
		new Config(null);
		GenomicCoords gc= new GenomicCoords("1:1-1000", 80, null, null);
		TrackIntervalFeature tif= new TrackIntervalFeature("test_data/CHD.exon.2010_03.sites.vcf", gc);
		tif.setNoFormat(true);
		tif.setTrackTag("title.bed");
		String[] lines= tif.concatTitleAndTrack().split("\n");
		assertEquals(1, lines.length);
		assertTrue( ! tif.concatTitleAndTrack().contains("\n"));
	}
	
	@Test
	public void canParsePrintableLinesWithSystemCommand() throws ClassNotFoundException, IOException, InvalidGenomicCoordsException, InvalidRecordException, SQLException, InvalidColourException, InvalidCommandLineException{

		GenomicCoords gc= new GenomicCoords("chr1:1-100000", 80, null, null);
		TrackIntervalFeature tif= new TrackIntervalFeature("test_data/hg19_genes_head.gtf", gc);
		tif.setNoFormat(true);
		tif.setPrintMode(PrintRawLine.FULL);
		tif.setSystemCommandForPrint("grep WASH7P | sort -k5,5nr | sed 's/WASH/wash/'");
		String out= tif.printLines();
		// Same as `awk '$4 <= 100000' hg19_genes_head.gtf | grep WASH7P | sort -k5,5nr | wc -l`
		assertEquals(11, Splitter.on("\n").omitEmptyStrings().splitToList(out).size()); // 11 lines grepped.
		assertTrue(out.contains("wash"));
	}

	@Test
	public void canParsePrintableLinesWithInvalidCommand() throws ClassNotFoundException, IOException, InvalidGenomicCoordsException, InvalidRecordException, SQLException, InvalidColourException, InvalidCommandLineException{

		// Invalid command: Empty output. But note that no exception is thrown!
		GenomicCoords gc = new GenomicCoords("chr1:1-100000", 80, null, null);
		TrackIntervalFeature tif = new TrackIntervalFeature("test_data/hg19_genes_head.gtf", gc);
		tif.setSystemCommandForPrint("foo");
		assertEquals("", tif.printLines());
	}
	
	@Test
	public void canHighlightLinesByIndex() throws ClassNotFoundException, IOException, InvalidGenomicCoordsException, InvalidRecordException, SQLException, InvalidColourException, InvalidCommandLineException{
		GenomicCoords gc= new GenomicCoords("1:1-400000", 80, null, null);
		TrackIntervalFeature tif= new TrackIntervalFeature("test_data/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf", gc);
		tif.setPrintMode(PrintRawLine.FULL);
		tif.setHighlightPattern(Pattern.compile("$2, $4"));
		String out= tif.printLines();
		assertTrue(out.contains("\033[7m200000\033[27m"));
		assertTrue(out.contains("\033[7mA\033[27m"));
		assertTrue(out.contains(" ALU_umary_ALU_2 ")); // Not highlighted

		tif.setHighlightPattern(Pattern.compile("$0, $100, 3"));
		out= tif.printLines();
		assertTrue(! out.contains("\033[7m"));
	}
	
	@Test
	public void canHighlightPrintableLines() throws ClassNotFoundException, IOException, InvalidGenomicCoordsException, InvalidRecordException, SQLException, InvalidColourException, InvalidCommandLineException{

		GenomicCoords gc= new GenomicCoords("1:1-400000", 80, null, null);
		TrackIntervalFeature tif= new TrackIntervalFeature("test_data/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf", gc);
		tif.setPrintMode(PrintRawLine.FULL);
		tif.setHighlightPattern(Pattern.compile("200000"));
		String out= tif.printLines();
		assertTrue(out.contains("\033[7m200000\033[27m"));
		
		// Multiple patterns
		tif.setHighlightPattern(Pattern.compile("ALU"));
		out= tif.printLines();
		assertTrue(out.contains("\033[7mALU\033[27m"));

		// Pattern start of line
		tif.setHighlightPattern(Pattern.compile("1"));
		tif.setPrintMode(PrintRawLine.FULL);
		out= tif.printLines();
		assertTrue(out.startsWith("\033[38;5;0;48;5;231m\033[7m1[27m"));
		
		// Pattern EOL
		tif.setHighlightPattern(Pattern.compile("1\\|1"));
		out= tif.printLines();
		assertTrue(out.trim().endsWith("1|1\033[27m"));
		
		// No pattern
		tif.setHighlightPattern(Pattern.compile("FOOBAR"));
		out= tif.printLines();
		assertTrue( ! out.contains("\033[7m"));
		
		// Reset
		tif.setHighlightPattern(Pattern.compile(""));
		out= tif.printLines();
		assertTrue( ! out.contains("\033[7m"));
		
		// Highlight FORMAT tag and associated values
		tif.setHighlightPattern(Pattern.compile("GT"));
		out= tif.printLines();
		assertTrue(out.contains("\033[7mGT\033[27m"));
		assertTrue(out.contains("\033[7m.|.\033[27m"));
		assertTrue(out.contains("\033[7m0|0\033[27m"));
		
		// Turn the string-table to back to a list and ensure that all fields are there
		out = Utils.stripAnsiCodes(out);
		List<String> xlist= Splitter.on("\t").omitEmptyStrings().splitToList(out.replaceAll(" \\| |\n", "\t"));
		assertEquals(24, xlist.size());
	}
	
	@Test
	public void canHighlightFormatInVcfLines() throws ClassNotFoundException, IOException, InvalidGenomicCoordsException, InvalidRecordException, SQLException, InvalidColourException, InvalidCommandLineException{

		GenomicCoords gc= new GenomicCoords("1:1-20000000", 80, null, null);
		TrackIntervalFeature tif= new TrackIntervalFeature("test_data/info_formats.vcf.gz", gc);
		tif.setPrintMode(PrintRawLine.FULL);

		// Highlight FORMAT tag and associated values
		tif.setHighlightPattern(Pattern.compile("GT"));
		String out= tif.printLines();
		assertTrue(out.contains("\033[7mGT\033[27m"));
		assertTrue(out.contains("\033[7m1/0\033[27m"));
		assertTrue(out.contains("\033[7m1|2\033[27m"));

//		// Tag GTFF is not highlighted because the match is partial
//		assertTrue(out.contains("\033[7m1|2\033[27m:100:0,0"));
//		
//		// This should match both:
//		tif.setHighlightPattern(Pattern.compile("GT|GTFF"));
//		// TODO: assertTrue(out.contains("\033[7m1|2\033[27m:100:0,0"));
		
		// Turn the string-table to back to a list and ensure that all fields are there
		out = Utils.stripAnsiCodes(out);
		List<String> xlist= Splitter.on("\t").omitEmptyStrings().splitToList(out.replaceAll(" \\| |\n", "\t"));
		assertEquals(33, xlist.size());
	}
	
//	@Test
//	public void canParsePrintableLinesWithSystemCommandBcftools() throws ClassNotFoundException, IOException, InvalidGenomicCoordsException, InvalidRecordException, SQLException, InvalidColourException, InvalidCommandLineException{
//
//		GenomicCoords gc= new GenomicCoords("1:1-400000", 80, null, null);
//		TrackIntervalFeature tif= new TrackIntervalFeature("test_data/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf", gc);
//		tif.setNoFormat(true);
//		tif.setPrintMode(PrintRawLine.FULL);
//		// Use full path to bcftools only for testing because Eclipse doesn't look into user's PATH
//		tif.setSystemCommandForPrint("/usr/local/bin/bcftools query -f '%CHROM\t%POS\t%AF\t[%GT:]\n'");
//		tif.setPrintNumDecimals(-1);
//		String out= tif.printLines();
//		assertEquals(2, out.split("\n").length);
//		assertTrue(out.contains("0.00698882"));
//	}
	
	@Test
	public void canParsePrintableLinesWithNoFeatures() throws ClassNotFoundException, IOException, InvalidGenomicCoordsException, InvalidRecordException, SQLException, InvalidColourException, InvalidCommandLineException{

		// Test region with no features
		GenomicCoords gc = new GenomicCoords("chr10:1-100000", 80, null, null);
		TrackIntervalFeature tif = new TrackIntervalFeature("test_data/hg19_genes_head.gtf", gc);
		tif.setNoFormat(true);
		tif.setPrintMode(PrintRawLine.FULL);
		tif.setSystemCommandForPrint("head");
		assertEquals("", tif.printLines());
	}

	@Test
	public void canParsePrintBAM() throws InvalidGenomicCoordsException, IOException, InvalidColourException, InvalidCommandLineException, ClassNotFoundException, InvalidRecordException, SQLException, InvalidConfigException{
		new Config(null);
		// BAM 
		GenomicCoords gc= new GenomicCoords("chr7:5566733-5566903", 80, null, null);
		TrackReads tif= new TrackReads("test_data/ds051.short.bam", gc);
		tif.setNoFormat(true);
		tif.setPrintMode(PrintRawLine.FULL);
		tif.setSystemCommandForPrint("grep NCNNTCCC");
		assertEquals(2, tif.printLines().split("\n").length);
	}

	@Test
	public void canPrintBamWithClippedSeq() throws InvalidGenomicCoordsException, IOException, InvalidColourException, InvalidCommandLineException, ClassNotFoundException, InvalidRecordException, SQLException, InvalidConfigException{
		
		new Config(null);
		GenomicCoords gc= new GenomicCoords("chr7:5566733-5566903", 160, null, null);
		TrackReads tif= new TrackReads("test_data/ds051.actb.bam", gc);
		tif.setNoFormat(true);
		tif.setPrintMode(PrintRawLine.CLIP);
		assertTrue(tif.printLines().contains(" CTCAT[+"));
		assertTrue(tif.printLines().contains(" CCCFF[+"));
	}

	@Test
	public void canPrintVcfWithClippedSeq() throws InvalidGenomicCoordsException, IOException, InvalidColourException, InvalidCommandLineException, ClassNotFoundException, InvalidRecordException, SQLException, InvalidConfigException{
		new Config(null);
		GenomicCoords gc= new GenomicCoords("1:1019492-1019672", 160, null, null);
		TrackIntervalFeature tif= new TrackIntervalFeature("test_data/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf", gc);
		tif.setNoFormat(true);
		tif.setPrintMode(PrintRawLine.CLIP);
		assertTrue(tif.printLines().contains(" GTCAC["));
	}
	
	@Test
	public void printIsNotResetAfterExec() throws InvalidGenomicCoordsException, IOException, InvalidColourException, InvalidCommandLineException, ClassNotFoundException, InvalidRecordException, SQLException, InvalidConfigException{
		new Config(null);
		// BAM 
		GenomicCoords gc= new GenomicCoords("chr7:5566733-5566903", 80, null, null);
		TrackReads tif= new TrackReads("test_data/ds051.short.bam", gc);
		tif.setNoFormat(true);
		tif.setPrintMode(PrintRawLine.FULL);
		tif.setSystemCommandForPrint("grep NCNNTCCC");
		assertEquals(2, tif.printLines().split("\n").length);
		
		// Call printLines again on new coordinates: 
		// The sys command is still on
		gc= new GenomicCoords("chr7:5566733-5566904", 80, null, null);
		tif.setGc(gc);
		System.err.println(tif.printLines());
		assertEquals(2, tif.printLines().split("\n").length);
		
		// Now turn it off:
		tif.setSystemCommandForPrint("");
		assertEquals(22, tif.printLines().split("\n").length);
	}

	
}
