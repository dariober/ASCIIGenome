package coloring;

import static org.junit.Assert.*;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;

//import org.fusesource.jansi.Ansi.Color;
import java.awt.Color;
import org.junit.Test;

import com.google.common.base.Splitter;

public class PngTest {

	@Test
	public void canGetColorFromStringStartingWithAnsiSequence(){
		
		// To convert RGB to HSB:
		// System.out.println( Arrays.toString(Color.RGBtoHSB(255, 182, 193, null)));
		
		Png png= new Png(new File("foo"));
		assertEquals("java.awt.Color[r=255,g=0,b=0]", png.ansiCodesToFgColor("\n[31mFOOBAR").toString());
		assertEquals("java.awt.Color[r=255,g=0,b=0]", png.ansiCodesToFgColor("[30;32;31;100mFOOBAR").toString());
		// No ansi escape: default
		assertEquals("java.awt.Color[r=0,g=0,b=0]", png.ansiCodesToFgColor("[31;FOOBAR").toString());
		// No foreground color: Default
		assertEquals("java.awt.Color[r=0,g=0,b=0]", png.ansiCodesToFgColor("[100;120mFOOBAR").toString());
		assertEquals("java.awt.Color[r=255,g=0,b=0]", png.ansiCodesToFgColor("[0;31mds051.actb.bam#1; ").toString());
		
	}
	
	@Test
	public void canPrintPngFromAnsiFile() throws IOException {
		
		File tmp= new File("test_data/ansicolor.png");
		
		Png png= new Png(new File("test_data/ansicolor.txt"));
		png.convert(tmp);
		System.out.println(tmp.getAbsolutePath());
		assertTrue(tmp.length() > 50000); // Check file size is about right.
	}

}
