package coloring;

import java.awt.Color;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.nio.file.StandardCopyOption;
import java.util.ArrayList;
import java.util.List;

import com.google.common.base.Splitter;
import com.itextpdf.text.BaseColor;
import com.itextpdf.text.Chunk;
import com.itextpdf.text.Document;
import com.itextpdf.text.DocumentException;
import com.itextpdf.text.Font;
import com.itextpdf.text.Paragraph;
import com.itextpdf.text.Rectangle;
import com.itextpdf.text.pdf.PdfCopy;
import com.itextpdf.text.pdf.PdfReader;
import com.itextpdf.text.pdf.PdfSmartCopy;
import com.itextpdf.text.pdf.PdfWriter;

public class Pdf {

	String ansiInput;
	/** Max height and width of the text to print. I think it's in "point" unit.   
	 * */
	int maxWidth= 0;
	int maxHeight= 0;
	
	/* C o n s t r u c t o r */
	
	public Pdf(String ansiInput) {
		// Possibly some validation about input file?
		this.ansiInput= ansiInput;
	}

	/* M e t h o d s */

	public void convert(File pdfOut, float fontSize, boolean append) throws IOException, DocumentException {

		// First we write to tmp file then we copy to given destination, possibly appending
		File tmpPdf= File.createTempFile(pdfOut.getName(), ".pdf");
		// tmpPdf.deleteOnExit();
		
		List<Paragraph> pdfLines= this.ansiFileToPdfParagraphs(fontSize);
        
		Document document = new Document(new Rectangle((float) (this.getMaxWidth() * 1.01), (float) (this.getMaxHeight())), 5f, 0f, 0f, 0f);
        PdfWriter.getInstance(document, new FileOutputStream(tmpPdf));
        document.open();

		for(Paragraph line : pdfLines){
	        document.add(line);
        }
        document.close();
        
        if(append){
        	this.appendPdf(tmpPdf, pdfOut);
        } else {
        	Files.move(Paths.get(tmpPdf.getAbsolutePath()), Paths.get(pdfOut.getAbsolutePath()), StandardCopyOption.REPLACE_EXISTING);
        }
	}

	private void appendPdf(File inFile, File outputFile) throws DocumentException, IOException {

		if(inFile.equals(outputFile)){
			// Consistent with Linux `cat test.txt >> test.txt`
			System.err.println("Cannot append pdf to itself. Input file is output file: " + inFile.getAbsolutePath());
			throw new RuntimeException(); 
		}
		
		if( ! outputFile.exists()){
			Files.move(Paths.get(inFile.getAbsolutePath()), Paths.get(outputFile.getAbsolutePath()));	
			return;
		}
		
		// It's silly to create a list of files first but we have a solution for concatenting lists of files.
		List<File> listOfPdfFiles= new ArrayList<File>();
		listOfPdfFiles.add(outputFile);
		listOfPdfFiles.add(inFile);
		
		File tmpPdf= File.createTempFile(outputFile.getName(), ".pdf");
		// tmpPdf.deleteOnExit();
		
		System.err.println(listOfPdfFiles);
		System.err.println(tmpPdf);
		
		// concatenatePdfs(listOfPdfFiles, tmpPdf);
		// Files.move(Paths.get(tmpPdf.getAbsolutePath()), Paths.get(outputFile.getAbsolutePath()), StandardCopyOption.REPLACE_EXISTING);
	}

	/* Concatenate pdf files. From 
	 * http://stackoverflow.com/questions/23062345/function-that-can-use-itext-to-concatenate-merge-pdfs-together-causing-some 
	 * **/
	public void concatenatePdfs(List<File> listOfPdfFiles, File outputFile) throws DocumentException, IOException {
        Document document = new Document();
        FileOutputStream outputStream = new FileOutputStream(outputFile);
        PdfCopy copy = new PdfSmartCopy(document, outputStream);
        document.open();

        for (File inFile : listOfPdfFiles) {
            PdfReader reader = new PdfReader(inFile.getAbsolutePath());
            copy.addDocument(reader);
            reader.close();
        }
        document.close();
	}
	
	/** Read ansi formatted file line by lone and convert each line to 
	 * a iText paragraph. 
	 * @throws IOException 
	 * */
	private List<Paragraph> ansiFileToPdfParagraphs(float fontSize) throws IOException{
		
		List<Paragraph> paraList= new ArrayList<Paragraph>();

		byte[] encoded = this.ansiInput.getBytes(); // Files.readAllBytes(Paths.get(this.ansiInput.toURI()));
		String str= new String(encoded);

		// Each string in this list is a sequence of a characters sharing the same formatting.
		List<String> ansiList= Splitter.on("\033").splitToList(str);  

		Paragraph pdfLine= new Paragraph();
		int nChars= 0;
		int currentMax= 0;
		for(String xv : ansiList){
			if(xv.equals("[48;5;231m")){
				continue;
			}
			BaseColor fgBaseCol= new BaseColor(this.ansiCodesToFgColor(xv).getRGB());
			BaseColor bgBaseCol= new BaseColor(this.ansiCodesToBgColor(xv).getRGB());

			if(this.extractAnsiCodes(xv).size() != 0){
				// This string begins with ansi sequence and the color has been extracted.
				// So remove the ansi sequence at the beginnig and the end, we don't need them anymore
				xv= xv.replaceAll("^.+?m", "");
			}
			
			for(int i= 0; i < xv.length(); i++){
				
				char c= xv.charAt(i);
				nChars += 1;
				Chunk chunk= new Chunk(c, new Font(Font.FontFamily.COURIER, fontSize, Font.NORMAL, fgBaseCol));
				chunk.setBackground(bgBaseCol);
				pdfLine.add(chunk);
				
				if(c == '\n'){
					// Newline found: Start a new parapgrah
					paraList.add(pdfLine);
					if(nChars > currentMax){
						currentMax= nChars;
					}
					nChars= 0;
					pdfLine= new Paragraph();
				}
			}
		}
		Font f= new Font(Font.FontFamily.COURIER, fontSize);
		Chunk chunk= new Chunk("X", f);
		this.setMaxWidth((int) ((currentMax) * chunk.getWidthPoint()));
		for(Paragraph p : paraList){
			p.setSpacingBefore(-3);
			p.setSpacingAfter(-3);
		}
		this.setMaxHeight((int) Math.rint((paraList.size() + 1) * fontSize));
		return paraList;
	}
	
	/**String x is expected to start with the ansi escape sequence. From the ansi sequence
	 * extract the (last) color for the foreground color.
	 * If x doesn't start with the ansi escape return a default color (black?). 
	 * */
	private Color ansiCodesToFgColor(String x){

		List<Integer> ansiCodes= this.extractAnsiCodes(x);

		Color col= Color.BLACK;
		for(int n : ansiCodes){
			if(ColoredChar.ansiForegroundColourToGraphicsColor(n) != null){
				// Code is a valid foreground color so update the current Color
				col= ColoredChar.ansiForegroundColourToGraphicsColor(n);
			}
		}
		return col;
	}
	
	/**Get the background color from the ansi sequence 
	 * */
	private Color ansiCodesToBgColor(String x){

		List<Integer> ansiCodes= this.extractAnsiCodes(x);

		Color col= Color.WHITE;
		for(int n : ansiCodes){
			if(ColoredChar.ansiBackgroundColourToGraphicsColor(n) != null){
				// Code is a valid background color so update the current Color
				col= ColoredChar.ansiBackgroundColourToGraphicsColor(n);
			}
		}
		return col;
	}
	
	/** Return a list the ansi codes found at the beginning of String x. 
	 * Not the the escape \033 has been already removed from x. 
	 *  extractAnsiCodes("[31;32;33mFOOBAR") -> [31, 32, 33]
	 *  extractAnsiCodes("FOOBAR") -> []
	 * */
	private List<Integer> extractAnsiCodes(String x){
		List<Integer> ansiCodes= new ArrayList<Integer>();
		x= x.replaceAll("\\s", "");
		if( ! (x.startsWith("[") || x.contains("m"))){
			// Beginning of x cannot be an ansi sequence 
			return ansiCodes;
		}
		List<String> digits = Splitter.on(";")
				.omitEmptyStrings()
				.splitToList(x.replaceAll("^\\[", "").replaceAll("m.*", ""));
		
		for(String d : digits){
			try{
				Integer.parseInt(d);
			} catch(NumberFormatException e){
				// Something between [ and m is not a digit. This is not an ansi sequence
				return new ArrayList<Integer>();
			}
			int n= Integer.parseInt(d);
			if(n < 0 || n > 999){
				return new ArrayList<Integer>();
			}
			ansiCodes.add(n);
		}
		return ansiCodes;
	}

//	private File getAnsiInput() {
//		return ansiInput;
//	}
//
//	private void setAnsiInput(File ansiInput) {
//		this.ansiInput = ansiInput;
//	}

	private int getMaxWidth() {
		return maxWidth;
	}

	private void setMaxWidth(int maxWidth) {
		this.maxWidth = maxWidth;
	}

	private int getMaxHeight() {
		return maxHeight;
	}

	private void setMaxHeight(int maxHeight) {
		this.maxHeight = maxHeight;
	}

	
}
