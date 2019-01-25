package samTextViewer;

import static org.junit.Assert.*;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.sql.SQLException;
import java.util.ArrayList;
import org.junit.Test;

import com.google.common.base.Splitter;

import exceptions.InvalidCommandLineException;
import exceptions.InvalidGenomicCoordsException;
import exceptions.InvalidRecordException;

public class InteractiveInputTest {

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
