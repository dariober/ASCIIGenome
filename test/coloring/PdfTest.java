package coloring;

import static org.junit.Assert.*;

import java.io.File;
import java.io.IOException;

import org.junit.Test;

import com.itextpdf.text.DocumentException;

public class PdfTest {

	@Test
	public void canPrintPdfFromAnsiFile() throws DocumentException, IOException{
		
		File tmp= new File("deleteme.pdf");
		tmp.deleteOnExit();
		Pdf pdf= new Pdf(new File("test_data/ansicolor.txt"));
		pdf.convert(tmp, 10);
		assertTrue(tmp.length() > 1000); // Check file size is about right.
	
	}
}
