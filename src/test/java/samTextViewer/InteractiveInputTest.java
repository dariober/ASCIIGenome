package samTextViewer;

import static org.junit.Assert.*;

import java.io.IOException;

import org.junit.Test;

import com.google.common.base.Splitter;

import exceptions.InvalidGenomicCoordsException;

public class InteractiveInputTest {

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
