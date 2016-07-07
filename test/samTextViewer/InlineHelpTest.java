package samTextViewer;

import static org.junit.Assert.*;

import org.junit.Test;

public class InlineHelpTest {

	@Test
	public void compileHelp(){
		System.out.println(InlineHelp.getHelp());
	}
	
	@Test
	public void allPlaceholdersReplaced(){
		// This is a dumb check that all placeholders are substituted
		String fmtHelp= InlineHelp.getHelp();
		assertTrue(!fmtHelp.contains("${"));
	}
	
}
