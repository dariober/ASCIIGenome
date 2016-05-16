package samTextViewer;

import static org.junit.Assert.*;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.junit.Test;

import exceptions.InvalidGenomicCoordsException;

public class GenomicCoordsHistoryTest {

	@Test
	public void canMoveBackAndForthInHistory() throws InvalidGenomicCoordsException, IOException{

		GenomicCoordsHistory gch= new GenomicCoordsHistory();
		GenomicCoords g1= new GenomicCoords("chr7:1-100", null, 100, null);
		GenomicCoords g2= new GenomicCoords("chr7:2-100", null, 100, null);
		GenomicCoords g3= new GenomicCoords("chr7:3-100", null, 100, null);
		GenomicCoords g4= new GenomicCoords("chr7:4-100", null, 100, null);
		GenomicCoords g5= new GenomicCoords("chr7:5-100", null, 100, null);
		gch.add(g1);
		gch.add(g2);
		gch.add(g3);
		gch.add(g4);
		gch.add(g5);
		
		assertTrue(g5.equalCoords(gch.current()));
		gch.previous();
		assertTrue(g4.equalCoords(gch.current()));
		gch.previous();
		assertTrue(g3.equalCoords(gch.current()));
		gch.previous();
		gch.previous();
		gch.previous();
		gch.previous();
		gch.previous();
		assertTrue(g1.equalCoords(gch.current()));
		
		gch.next();
		assertTrue(g2.equalCoords(gch.current()));
		gch.next();
		gch.next();
		gch.next();
		gch.next();
		gch.next();
		gch.next();
		assertTrue(g5.equalCoords(gch.current()));
		
		// Move back and ensure adding an item moves the tracker to last position
		gch.previous();
		gch.previous();
		gch.previous();
		GenomicCoords g6= new GenomicCoords("chr7:6-100", null, 100, null);
		gch.add(g6);
		assertTrue(g6.equalCoords(gch.current()));
		
		// Adding same pos more than once has no effect:
		GenomicCoords g7= new GenomicCoords("chr7:6-100", null, 100, null);
		gch.add(g7);
		GenomicCoords g8= new GenomicCoords("chr7:6-100", null, 100, null);
		gch.add(g8);
		// System.out.println(gch.getHistory());
		
		GenomicCoordsHistory gchNull= new GenomicCoordsHistory();
		assertNull(gchNull.current());
		//assertNull(gchNull.previous());
		//assertNull(gchNull.next());
	
	}

	@Test
	public void getAndPutItems() throws InvalidGenomicCoordsException, IOException{
		
		GenomicCoords g1= new GenomicCoords("chr7:1-100", null, 100, null);
		
		GenomicCoordsHistory gch= new GenomicCoordsHistory();
		gch.add(g1);

		// Get an item, change it, put it back as a different one by *cloning*:
		GenomicCoords g1z= (GenomicCoords) gch.current().clone();
		g1z.zoomOut();
		gch.add(g1z);
		assertEquals(2, gch.getHistory().size());
		// System.out.println(gch.getHistory());
	}
	
}
