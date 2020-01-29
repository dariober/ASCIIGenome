package coloring;

import java.awt.Color;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.nio.file.StandardCopyOption;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
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

import exceptions.InvalidColourException;
import samTextViewer.Utils;

public class Pdf {

	String ansiInput;
	/** Max height and width of the text to print. I think it's in "point" unit.   
	 * */
	int maxWidth= 0;
	int maxHeight= 0;
	
	/* C o n s t r u c t o r */
	
	public Pdf(String ansiInput) {
		// Possibly some validation about input file?
		// Input string must end with newline otherwise the last line is not going to be printed to pdf.
		if( ! ansiInput.endsWith("\n")){
			ansiInput += '\n';
		}
		this.ansiInput= ansiInput;
	}

	/* M e t h o d s */

	public void convert(File pdfOut, float fontSize, boolean append) throws IOException, DocumentException, InvalidColourException {

		// First we write to tmp file then we copy to given destination, possibly appending
		File tmpPdf= Utils.createTempFile(pdfOut.getName(), ".pdf", true);
		
		List<Paragraph> pdfLines= this.ansiFileToPdfParagraphs(fontSize);
        
		Rectangle pageSize = new Rectangle((float) (this.getMaxWidth() * 1.01), (float) (this.getMaxHeight()));

		int background256= Config.get256Color(ConfigKey.background);
		Color pageColor= Xterm256.xterm256ToColor(background256);
		
		pageSize.setBackgroundColor(new BaseColor(pageColor.getRed(), pageColor.getGreen(), pageColor.getBlue()));
		Document document = new Document(pageSize, 5f, 0f, 0f, 0f);
		//Document document = new Document(new Rectangle((float) (this.getMaxWidth() * 1.01), (float) (this.getMaxHeight())), 5f, 0f, 0f, 0f);
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

	/** Append pdf file `in` to pdf file `dest`
	 * @throws IOException 
	 * @throws DocumentException 
	 * */
	private void appendPdf(File in, File dest) throws IOException, DocumentException{
		
		if( ! dest.exists()){
			Files.move(Paths.get(in.getAbsolutePath()), Paths.get(dest.getAbsolutePath()));
			return;
		}
		
		File template= Utils.createTempFile(".template.", ".pdf", true);
		
		Files.copy(Paths.get(dest.getAbsolutePath()), Paths.get(template.getAbsolutePath()), StandardCopyOption.REPLACE_EXISTING);

		Document document = new Document();
        FileOutputStream outputStream = new FileOutputStream(template);
        PdfCopy copy = new PdfSmartCopy(document, outputStream);
        document.open();

        PdfReader reader0 = new PdfReader(dest.getAbsolutePath());
        copy.addDocument(reader0);
        reader0.close();
        
        PdfReader reader = new PdfReader(in.getAbsolutePath());
        copy.addDocument(reader);
        reader.close();

        document.close();
        Files.move(Paths.get(template.getAbsolutePath()), Paths.get(dest.getAbsolutePath()), StandardCopyOption.REPLACE_EXISTING);
		
	}
	
	/** Read ansi formatted file line by lone and convert each line to 
	 * a iText paragraph. 
	 * @throws IOException 
	 * @throws InvalidColourException 
	 * */
	private List<Paragraph> ansiFileToPdfParagraphs(float fontSize) throws IOException, InvalidColourException{
		
		List<Paragraph> paraList= new ArrayList<Paragraph>();

		byte[] encoded = this.ansiInput.getBytes(); // Files.readAllBytes(Paths.get(this.ansiInput.toURI()));
		String str= new String(encoded);

		// Each string in this list is a sequence of a characters sharing the same formatting.
		List<String> ansiList= Splitter.on("\033").splitToList(str);  

		Paragraph pdfLine= new Paragraph();
		int nChars= 0;
		int currentMax= 0;
		for(String xv : ansiList){
			
			BaseColor bgBaseCol;
			BaseColor fgBaseCol;
			List<Integer> extractAnsi = this.extractAnsiCodes(xv);

			if(extractAnsi.contains(7)) {
				// Invert colour if ansi formatting says so
				bgBaseCol= new BaseColor(this.xterm256ToColor(xv, false).getRGB());
				fgBaseCol= new BaseColor(this.xterm256ToColor(xv, true).getRGB());
			} else {
				fgBaseCol= new BaseColor(this.xterm256ToColor(xv, false).getRGB());
				bgBaseCol= new BaseColor(this.xterm256ToColor(xv, true).getRGB());
			}
			
			if(extractAnsi.size() != 0){
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

	/** Parse the string x to get the colour for foreground or background. 
	 * If the input string doesn't contain the escape sequence for fore or back ground,
	 * use the colour from Config. 
	 * This method should be private. It is protected only for unit test.
	 * */
	protected Color xterm256ToColor(String x, boolean isBackground) throws InvalidColourException{

		List<Integer> ansiCodes= this.extractAnsiCodes(x);
		
		Color configBg= Xterm256.xterm256ToColor(Config.get256Color(ConfigKey.background));
		Color configFg= Xterm256.xterm256ToColor(Config.get256Color(ConfigKey.foreground));
		
		Color col= null;
		int xtag= -1;
		if(isBackground){
			col= configBg; // Color.WHITE; // Default background
			xtag= Collections.indexOfSubList(ansiCodes, Arrays.asList(new Integer[] {48, 5}));
		} else {
			col= configFg; // Color.BLACK; // Default foreground
			xtag= Collections.indexOfSubList(ansiCodes, Arrays.asList(new Integer[] {38, 5}));
		}
		
		if(xtag != -1 && ansiCodes.size() > 2){
			col= Xterm256.xterm256ToColor(ansiCodes.get(xtag + 2));
		}
		
		return col;
	}
	
	/** Return a list the ansi codes found at the beginning of String x. 
	 * Not the the escape \033 has been already removed from x. 
	 *  extractAnsiCodes("[31;32;33mFOOBAR") -> [31, 32, 33]
	 *  extractAnsiCodes("FOOBAR") -> []
	 *  engineer 
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
