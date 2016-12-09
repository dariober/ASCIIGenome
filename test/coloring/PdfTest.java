package coloring;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;

import org.apache.commons.io.FileUtils;
import org.junit.Test;

import com.itextpdf.text.Document;
import com.itextpdf.text.DocumentException;
import com.itextpdf.text.pdf.PdfCopy;
import com.itextpdf.text.pdf.PdfReader;
import com.itextpdf.text.pdf.PdfSmartCopy;

public class PdfTest {

	@Test
	public void testAppend() throws DocumentException, IOException{
		
		String existing= "test_data/deleteme.1.pdf";
		String toAppend= "test_data/deleteme.2.pdf";
		
        Document document = new Document();
        FileOutputStream outputStream = new FileOutputStream(existing);
        PdfCopy copy = new PdfSmartCopy(document, outputStream);
        document.open();

        PdfReader reader = new PdfReader(toAppend);
        copy.addDocument(reader);
        reader.close();
        document.close();
		
	}
	
	@Test
	public void canPrintPdfFromAnsiFile() throws DocumentException, IOException{
		
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
	public void canAppendToPdf() throws IOException, DocumentException{

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
