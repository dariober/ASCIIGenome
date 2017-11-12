package tracks;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.IOException;
import java.sql.SQLException;

import org.junit.Test;

import com.google.common.base.Splitter;

import coloring.Config;
import exceptions.InvalidColourException;
import exceptions.InvalidCommandLineException;
import exceptions.InvalidConfigException;
import exceptions.InvalidGenomicCoordsException;
import exceptions.InvalidRecordException;
import samTextViewer.GenomicCoords;

public class TrackTest {

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
		tif.setSystemCommandForPrint("grep WASH7P | sort -k5,5nr");
		String out= tif.printLines();
		// Same as `awk '$4 <= 100000' hg19_genes_head.gtf | grep WASH7P | sort -k5,5nr | wc -l`
		assertEquals(11, Splitter.on("\n").omitEmptyStrings().splitToList(out).size()); // 11 lines grepped.
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
	public void canParsePrintableLinesWithSystemCommandBcftools() throws ClassNotFoundException, IOException, InvalidGenomicCoordsException, InvalidRecordException, SQLException, InvalidColourException, InvalidCommandLineException{

		GenomicCoords gc= new GenomicCoords("1:1-400000", 80, null, null);
		TrackIntervalFeature tif= new TrackIntervalFeature("test_data/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf", gc);
		tif.setNoFormat(true);
		tif.setPrintMode(PrintRawLine.FULL);
		// Use full path to bcftools only for testing because Eclipse doesn't look into user's PATH
		tif.setSystemCommandForPrint("/usr/local/bin/bcftools query -f '%CHROM\t%POS\t%AF\t[%GT:]\n'");
		String out= tif.printLines();
		assertEquals(2, out.split("\n").length);
		assertTrue(out.contains("0.00698882"));
	}
	
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
		assertEquals(2, tif.printLines().split("\n").length);
		
		// Now turn it off:
		tif.setSystemCommandForPrint("");
		assertEquals(22, tif.printLines().split("\n").length);
	}

	
}
