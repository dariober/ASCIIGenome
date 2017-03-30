package coloring;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.io.IOException;

import org.apache.commons.io.FileUtils;
import org.junit.Test;

import com.itextpdf.text.DocumentException;
import com.itextpdf.text.pdf.PdfReader;

import exceptions.InvalidColourException;
import exceptions.InvalidConfigException;

public class PdfTest {

	@Test
	public void canProcessAllLines() throws IOException, InvalidColourException, DocumentException{

		Pdf withoutNewline= new Pdf("foo\nbar\nbaz");
		Pdf withNewline= new Pdf("foo\nbar\nbaz\n");
		File tmp1= new File("deleteme1.pdf");
		File tmp2= new File("deleteme2.pdf");
		tmp1.deleteOnExit();
		tmp2.deleteOnExit();
		withoutNewline.convert(tmp1, 10, false);
		withNewline.convert(tmp2, 10, false);
		assertEquals(tmp1.length(), tmp2.length());
	}
	
	@Test
	public void getColorFromAnsi() throws IOException, InvalidColourException, InvalidConfigException{

		new Config(null);
		
		// Dummy pdf object. It doesn't matter how you create it.
		String ansiInput= FileUtils.readFileToString(new File("test_data/ansicolor.txt"));
		Pdf pdf= new Pdf(ansiInput);


		// Expected: xterm 244 = grey with rgb= 128
		String x= "[38;5;244;48;5;15m FOO";
		assertEquals("java.awt.Color[r=128,g=128,b=128]", pdf.xterm256ToColor(x, false).toString());
		
		x= "[48;5;15;38;5;244m FOO";
		assertEquals("java.awt.Color[r=128,g=128,b=128]", pdf.xterm256ToColor(x, false).toString());
		
		// Colour for foreground not given, default to Config:
		x= "[48;5;15m FOO";
		assertEquals("java.awt.Color[r=0,g=0,b=0]", pdf.xterm256ToColor(x, false).toString());
		// Background as in string
		assertEquals("java.awt.Color[r=255,g=255,b=255]", pdf.xterm256ToColor(x, true).toString());
	}
	
	@Test
	public void canPrintPdfFromAnsiFile() throws DocumentException, IOException, InvalidColourException{
		
		String ansiInput= FileUtils.readFileToString(new File("test_data/ansicolor.txt"));
		
		File tmp= new File("deleteme.pdf");
		tmp.deleteOnExit();
		Pdf pdf= new Pdf(ansiInput);
		pdf.convert(tmp, 10, false);
		assertTrue(tmp.length() > 1000); // Check file size is about right.
	
		PdfReader p = new PdfReader(tmp.getAbsolutePath());
		assertEquals(1, p.getNumberOfPages());
		
	}
	
	@Test 
	public void canAppendToPdf() throws IOException, DocumentException, InvalidColourException{

		File tmp= new File("append.pdf");
		tmp.delete();
		tmp.deleteOnExit();
		String ansiInput= FileUtils.readFileToString(new File("test_data/ansicolor.txt"));
		Pdf pdf= new Pdf(ansiInput);

		// Append to file not yet created:
		pdf.convert(tmp, 10, true);

		PdfReader p = new PdfReader(tmp.getAbsolutePath());
		assertEquals(1, p.getNumberOfPages());
		
		// Append another file 
		pdf.convert(tmp, 10, true);
		p = new PdfReader(tmp.getAbsolutePath());
		assertEquals(2, p.getNumberOfPages()); // We got two pages right?
	}
}
