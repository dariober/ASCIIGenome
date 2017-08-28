package samTextViewer;

import static org.junit.Assert.*;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import org.junit.Test;

public class ASCIIGenomeHistoryTest {

	@Test
	public void canReadHistory() throws IOException {
		ASCIIGenomeHistory ag= new ASCIIGenomeHistory("test_data/asciigenome.yaml");
		assertEquals(8, ag.getFiles().size());
		assertEquals(4, ag.getPositions().size());
		assertEquals(6, ag.getCommandHistory().size());
		
		ag= new ASCIIGenomeHistory("non-existing.yaml");
		assertEquals(0, ag.getFiles().size());

		ag= new ASCIIGenomeHistory();
		assertNotNull(ag.getFileName());
	}

	@Test
	public void canPrepareEmptyObject() throws IOException{
		ASCIIGenomeHistory ag= new ASCIIGenomeHistory(null);
		assertEquals(0, ag.getFiles().size());
	}
	
	@Test
	public void canWriteHistory() throws IOException {
		ASCIIGenomeHistory ag= new ASCIIGenomeHistory("test_data/asciigenome.yaml");
		
		File out= new File("tmp.yaml");
		out.deleteOnExit();
		ag.write(out);
		assertTrue(out.length() > 100);
	}

}
