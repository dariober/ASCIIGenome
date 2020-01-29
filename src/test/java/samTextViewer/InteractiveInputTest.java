package samTextViewer;

import static org.junit.Assert.*;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.sql.SQLException;
import java.util.ArrayList;

import org.junit.Test;

import com.google.common.base.Splitter;

import coloring.Config;
import exceptions.InvalidCommandLineException;
import exceptions.InvalidConfigException;
import exceptions.InvalidGenomicCoordsException;
import exceptions.InvalidRecordException;
import jline.console.ConsoleReader;
import tracks.TrackSet;

public class InteractiveInputTest {
	
	private TrackProcessor gimmeTrackProcessor(String region, int terminalWidth) throws InvalidGenomicCoordsException, IOException, ClassNotFoundException, InvalidRecordException, SQLException {
		GenomicCoords gc= new GenomicCoords(region, terminalWidth, null, "test_data/chr7.fa");
		GenomicCoordsHistory gch= new GenomicCoordsHistory();
		gch.add(gc);
		TrackSet trackSet= new TrackSet(new ArrayList<String>(), gc);
		TrackProcessor proc= new TrackProcessor(trackSet, gch); 
		return proc;
	}
	
	@Test 
	public void canMovePositionByColumn() throws InvalidGenomicCoordsException, IOException, InvalidRecordException, ClassNotFoundException, SQLException, InvalidCommandLineException, InvalidConfigException {

		new Config(null);
		TrackProcessor proc;
		
		InteractiveInput ip= new InteractiveInput(new ConsoleReader());

		proc= this.gimmeTrackProcessor("chr7:1001-1800", 80);
		GenomicCoords gc2= ip.processInput("[", proc, 1).getGenomicCoordsHistory().current();
		assertEquals(1011, (int)gc2.getFrom());
		assertEquals(1810, (int)gc2.getTo());
		
		proc= this.gimmeTrackProcessor("chr7:1001-1800", 80);
		gc2= ip.processInput("]", proc, 1).getGenomicCoordsHistory().current();
		assertEquals(1001-10, (int)gc2.getFrom());
		assertEquals(1800-10, (int)gc2.getTo());

		proc= this.gimmeTrackProcessor("chr7:1001-1800", 80);
		gc2= ip.processInput("[ 20", proc, 1).getGenomicCoordsHistory().current();
		assertEquals(1001 + (20*10), (int)gc2.getFrom());
		assertEquals(1800 + (20*10), (int)gc2.getTo());

		proc= this.gimmeTrackProcessor("chr7:1001-1800", 80);
		gc2= ip.processInput("[20", proc, 1).getGenomicCoordsHistory().current();
		assertEquals(1001 + (20*10), (int)gc2.getFrom());
		assertEquals(1800 + (20*10), (int)gc2.getTo());
		
		proc= this.gimmeTrackProcessor("chr7:1001-1800", 80);
		gc2= ip.processInput("]] 3", proc, 1).getGenomicCoordsHistory().current();
		assertEquals(1001 - (6 * 10), (int)gc2.getFrom());
		assertEquals(1800 - (6 * 10), (int)gc2.getTo());
	
		proc= this.gimmeTrackProcessor("chr7:1001-1800", 80);
		gc2= ip.processInput("] 0", proc, 1).getGenomicCoordsHistory().current();
		assertEquals(1001, (int)gc2.getFrom());
		assertEquals(1800, (int)gc2.getTo());
		
		// Test left bound
		proc= this.gimmeTrackProcessor("chr7:1001-1800", 80);
		gc2= ip.processInput("]] 30000", proc, 1).getGenomicCoordsHistory().current();
		assertEquals(1, (int)gc2.getFrom());
		assertEquals(800, (int)gc2.getTo());

		// Test right bound
		proc= this.gimmeTrackProcessor("chr7:1001-1800", 80);
		gc2= ip.processInput("[ 30000000", proc, 1).getGenomicCoordsHistory().current();
		assertEquals(159138663, (int)gc2.getTo());
		assertEquals(159138663-800+1, (int)gc2.getFrom());
		
		// Invalid input
		ip.processInput("[]", proc, 1).getGenomicCoordsHistory().current();
		assertEquals(ExitCode.ERROR, ip.getInteractiveInputExitCode());

		// Ensure this is fine
		ip.processInput("]", proc, 1).getGenomicCoordsHistory().current();
		assertEquals(ExitCode.CLEAN, ip.getInteractiveInputExitCode());
		
		// Another invalid input
		ip.processInput("] foo", proc, 1).getGenomicCoordsHistory().current();
		assertEquals(ExitCode.ERROR, ip.getInteractiveInputExitCode());
	}
	
	
	@Test
	public void canPrintHelp() throws InvalidGenomicCoordsException, IOException, ClassNotFoundException, InvalidRecordException, SQLException, InvalidCommandLineException {
		
		InteractiveInput ip= new InteractiveInput(null);
		
		GenomicCoords gc= new GenomicCoords("chr7:1-100", 80, null, null);
		GenomicCoordsHistory gch= new GenomicCoordsHistory();
		gch.add(gc);
		TrackProcessor proc= new TrackProcessor(null, gch); 

		// Send the output of the help from stderr to baos
		ByteArrayOutputStream baos = new ByteArrayOutputStream();
		System.setErr(new PrintStream(baos));
		
		// Various ways of getting general help
		ip.processInput("-h", proc, 1);
		String H1= baos.toString(); baos.reset(); 
		assertTrue(H1.contains("show this help"));
		
		System.out.println(H1);
		
		ip.processInput("h", proc, 1);
		String H2= baos.toString(); baos.reset(); 
		
		ip.processInput("help", proc, 1);
		String H3= baos.toString(); baos.reset(); 
		
		ip.processInput("  ?", proc, 1);
		String H4= baos.toString(); baos.reset(); 
		
		assertEquals(H1, H2); 
		assertEquals(H1, H3);
		assertEquals(H1, H4);
		
		// Various way of getting command help
		ip.processInput("next -h", proc, 1);
		String h1= baos.toString(); baos.reset(); 
		assertTrue(h1.contains("Move to the next"));
		
		ip.processInput("?next", proc, 1);
		String h2= baos.toString(); baos.reset(); 
		
		ip.processInput("  help  next ", proc, 1);
		String h3= baos.toString(); baos.reset(); 

		assertEquals(h1, h2); 
		assertEquals(h1, h3);
	}
	
	@Test
	public void canGoToRegion() throws InvalidGenomicCoordsException, IOException {
		
		InteractiveInput ip= new InteractiveInput(null);
		
		GenomicCoords gc= new GenomicCoords("chr1:1-1000", 80, null, null);
		String region= ip.gotoOnCurrentChrom(Splitter.on(" ").splitToList("1 100"), gc);
		assertEquals("chr1:1-100", region);
		
		region= ip.gotoOnCurrentChrom(Splitter.on(" ").splitToList("1 20 30 100"), gc);
		assertEquals("chr1:1-100", region);
		
		region= ip.gotoOnCurrentChrom(Splitter.on(" ").splitToList("2"), gc);
		assertEquals("chr1:2-81", region);
	}

	@Test
	public void canGoToRegionAndCenter() throws InvalidGenomicCoordsException, IOException {
		
		InteractiveInput ip= new InteractiveInput(null);
		
		GenomicCoords gc= new GenomicCoords("chr1:1-20", 20, null, null);
		
		String region= ip.gotoOnCurrentChrom(Splitter.on(" ").splitToList("16 c"), gc);
		assertEquals("chr1:6-25", region);

		gc= new GenomicCoords("chr1:1-21", 20, null, null);
		
		region= ip.gotoOnCurrentChrom(Splitter.on(" ").splitToList("16 c"), gc);
		assertEquals("chr1:6-26", region);
	}
	
}
