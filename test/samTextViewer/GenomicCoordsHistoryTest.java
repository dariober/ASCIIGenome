package samTextViewer;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.junit.Test;

import exceptions.InvalidCommandLineException;
import exceptions.InvalidGenomicCoordsException;

public class GenomicCoordsHistoryTest {

	@Test
	public void canReadHistoricPositions() throws InvalidGenomicCoordsException, IOException{

		// This is the initial coordinates. Set by the user or reading the input files.
		// The historic positions will be checked against this coordinates.
		GenomicCoords checkGc= new GenomicCoords("chr1:1-10000", 80, null, "test_data/chr7.fa");

		GenomicCoordsHistory gch= new GenomicCoordsHistory();
		gch.readHistory(new File("test_data/asciigenome_history"), checkGc);
		assertEquals(2, gch.getHistory().size());
		assertEquals("chr7:10-1000", gch.getHistory().get(0).toStringRegion());
		
		// Initial genomic coords has no sequence dict, so everything in the history
		// file is loaded:
		checkGc= new GenomicCoords("chr1:1-10000", 80, null, null);

		gch= new GenomicCoordsHistory();
		gch.readHistory(new File("test_data/asciigenome_history"), checkGc);
		assertEquals(4, gch.getHistory().size());
		assertEquals("bar:10-1000", gch.getHistory().get(0).toStringRegion());

		// Behavior with missing or file
		gch= new GenomicCoordsHistory();
		gch.readHistory(new File("test_data/nonsense"), checkGc);
	}
	
	@Test
	public void canGetHistoryAsStringForHistoryFile() throws InvalidGenomicCoordsException, IOException{

		GenomicCoordsHistory gch= new GenomicCoordsHistory();
		// History is mixture of positions read from file, some of which invalid, and
		// positions visited in this session.
		gch.readHistory(new File("test_data/asciigenome_history"), new GenomicCoords("chr7:1-100", 80, null, "test_data/chr7.fa"));
		GenomicCoords g1= new GenomicCoords("chr7:1-100", 80, null, "test_data/chr7.fa");
		GenomicCoords g2= new GenomicCoords("chr7:2-100", 80, null, "test_data/chr7.fa");
		GenomicCoords g3= new GenomicCoords("chr7:3-100", 80, null, "test_data/chr7.fa");
		GenomicCoords g4= new GenomicCoords("chr7:4-100", 80, null, "test_data/chr7.fa");
		GenomicCoords g5= new GenomicCoords("chr7:5-100", 80, null, "test_data/chr7.fa");
		gch.add(g1);
		gch.add(g2);
		gch.add(g3);
		gch.add(g4);
		gch.add(g5);
		
		int maxPos= 100;
		List<String> hist= gch.prepareHistoryForHistoryFile(new File("test_data/asciigenome_history"), maxPos);
		assertEquals(9, hist.size()); // 2 valid pos from history + 2 invalid + 5 from current session
		assertEquals("## pos ##bar:10-1000", hist.get(0));
		assertEquals("## pos ##chr7:5-100", hist.get(hist.size()-1));

		maxPos= 2;
		hist= gch.prepareHistoryForHistoryFile(new File("test_data/asciigenome_history"), maxPos);
		assertEquals(maxPos, hist.size());
		assertEquals("## pos ##chr7:4-100", hist.get(0));
		assertEquals("## pos ##chr7:5-100", hist.get(1));

		maxPos= 0; // Don't write history
		hist= gch.prepareHistoryForHistoryFile(new File("test_data/asciigenome_history"), maxPos);
		assertEquals(0, hist.size());
		
		// Behave nicely with missing file
		hist= gch.prepareHistoryForHistoryFile(new File("test_data/nonsense"), 100);

	}
	
	@Test
	public void canMoveBackAndForthInHistory() throws InvalidGenomicCoordsException, IOException{

		GenomicCoordsHistory gch= new GenomicCoordsHistory();
		GenomicCoords g1= new GenomicCoords("chr7:1-100", 80, null, null);
		GenomicCoords g2= new GenomicCoords("chr7:2-100", 80, null, null);
		GenomicCoords g3= new GenomicCoords("chr7:3-100", 80, null, null);
		GenomicCoords g4= new GenomicCoords("chr7:4-100", 80, null, null);
		GenomicCoords g5= new GenomicCoords("chr7:5-100", 80, null, null);
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
		GenomicCoords g6= new GenomicCoords("chr7:6-100", 80, null, null);
		gch.add(g6);
		assertTrue(g6.equalCoords(gch.current()));
		
		// Adding same pos more than once has no effect:
		GenomicCoords g7= new GenomicCoords("chr7:6-100", 80, null, null);
		gch.add(g7);
		GenomicCoords g8= new GenomicCoords("chr7:6-100", 80, null, null);
		gch.add(g8);
		// System.out.println(gch.getHistory());
		
		GenomicCoordsHistory gchNull= new GenomicCoordsHistory();
		assertNull(gchNull.current());
		//assertNull(gchNull.previous());
		//assertNull(gchNull.next());
	
	}

	
	@Test
	public void canMoveBackAndForthInHistoryWithInvalidPosition() throws InvalidGenomicCoordsException, IOException, InvalidCommandLineException{

		GenomicCoordsHistory gch= new GenomicCoordsHistory();
		gch.setGenome(Utils.tokenize("test_data/ds051.actb.bam", "\t"));
		
		GenomicCoords g1= new GenomicCoords("foo:1-100", 80, null, null);
		GenomicCoords g2= new GenomicCoords("chr7:1-100", 80, null, null);
		gch.add(g1);
		gch.add(g2);
		
		gch.previous();
				
	}

	
	@Test
	public void getAndPutItems() throws InvalidGenomicCoordsException, IOException{
		
		GenomicCoords g1= new GenomicCoords("chr7:1-100", 80, null, null);
		
		GenomicCoordsHistory gch= new GenomicCoordsHistory();
		gch.add(g1);

		// Get an item, change it, put it back as a different one by *cloning*:
		GenomicCoords g1z= (GenomicCoords) gch.current().clone();
		g1z.zoomOut();
		gch.add(g1z);
		assertEquals(2, gch.getHistory().size());
		// System.out.println(gch.getHistory());
	}
	
	@Test
	public void canSetGenome() throws InvalidGenomicCoordsException, IOException, InvalidCommandLineException{
		
		GenomicCoordsHistory gch= new GenomicCoordsHistory();
		gch.add(new GenomicCoords("chr7:1-100", 80, null, null));
		List<String> cmdInput= new ArrayList<String>();
		cmdInput.add("test_data/chr7.fa");
		gch.setGenome(cmdInput);
		assertEquals("test_data/chr7.fa", gch.current().getFastaFile());
		assertTrue(gch.current().getSamSeqDict().toString().length() > 10);
		
		gch= new GenomicCoordsHistory();
		gch.add(new GenomicCoords("chr7:1-100", 80, null, null));
		gch.add(new GenomicCoords("chr7:1-1000", 80, null, null));
		gch.add(new GenomicCoords("chr7:1-10000", 80, null, null));
		cmdInput= new ArrayList<String>();
		cmdInput.add("hg19");
		gch.setGenome(cmdInput);
		assertTrue(gch.current().getSamSeqDict().toString().length() > 10);
		assertTrue(gch.getHistory().get(0).getSamSeqDict().toString().length() > 10);

	}

	@Test
	public void canSetGenomeWithInvalidInput() throws InvalidCommandLineException, IOException, InvalidGenomicCoordsException{
		
		GenomicCoordsHistory gch= new GenomicCoordsHistory();
		gch.add(new GenomicCoords("chr7:1-100", 80, null, null));

		// Input string is neither not an exiting file or a genome tag:
		List<String> cmdInput= new ArrayList<String>();
		cmdInput.add("test_data/foo.fa");
		gch.setGenome(cmdInput);
		
	}
	
	@Test
	public void canSetGenomeFromInvalidPosition() throws InvalidGenomicCoordsException, IOException, InvalidCommandLineException{
		
		GenomicCoordsHistory gch= new GenomicCoordsHistory();
		gch.add(new GenomicCoords("nonexisting:1-100", 80, null, null));
		gch.add(new GenomicCoords("nonexisting:2-100", 80, null, null));
		gch.add(new GenomicCoords("nonexisting:3-100", 80, null, null));
		
		List<String> cmdInput= new ArrayList<String>();
		cmdInput.add("test_data/chr7.fa");
		gch.setGenome(cmdInput);
		
		assertEquals("chr7", gch.current().getChrom());
		assertEquals("test_data/chr7.fa", gch.current().getFastaFile());
		assertEquals(1, gch.getHistory().size()); // Invalid positions have been removed.
	}
}
