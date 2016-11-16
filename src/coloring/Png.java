package coloring;

import java.awt.Color;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;

import javax.imageio.ImageIO;

import com.google.common.base.Splitter;

public class Png {

	File ansiInput;
	
	/* C o n s t r u c t o r */
	
	public Png(File ansiInput) {
		// Possibly some validation about input file?
		this.ansiInput= ansiInput;
	}

	/* M e t h o d s */

	public void convert(File pngOut, int fontSize) throws IOException {
		
		List<ColoredChar> coloredList= this.ansiFileToColoredChar();
		BufferedImage img= this.colorListToGraphics(coloredList, fontSize);
		ImageIO.write(img, "png", pngOut);
	}

	private BufferedImage colorListToGraphics(List<ColoredChar> coloredList, int fontSize) {

		BufferedImage img = new BufferedImage(1, 1, BufferedImage.TYPE_INT_ARGB);
        Graphics2D g2d = img.createGraphics();

        Font font= new Font("Courier", Font.BOLD, fontSize);
        
    	/* Width of the image and number of lines */
    	int width= -1;
    	int height= 0;
    	for(ColoredChar x : coloredList){
    		if(x.getX() > width){
    			width= x.getX();
    		}
    		if(x.getY() > height){
    			height = x.getY();
    		}
    	}
        g2d.setFont(font);
        FontMetrics fm = g2d.getFontMetrics();
    	g2d.dispose();
        
        width= (int) Math.round(width * fm.charWidth('X') * 1.01);
        height= (int) Math.round(height * fm.getAscent() * 1.5);
        
        img = new BufferedImage(width, height, BufferedImage.TYPE_INT_ARGB);
        g2d = img.createGraphics();
        g2d.setColor(Color.WHITE); // Set background colour by filling a rect
        g2d.fillRect(0, 0, width * fm.charWidth('X'), height);

        //g2d.setRenderingHint(RenderingHints.KEY_ALPHA_INTERPOLATION, RenderingHints.VALUE_ALPHA_INTERPOLATION_QUALITY);
        //g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
        //g2d.setRenderingHint(RenderingHints.KEY_COLOR_RENDERING, RenderingHints.VALUE_COLOR_RENDER_QUALITY);
        //g2d.setRenderingHint(RenderingHints.KEY_DITHERING, RenderingHints.VALUE_DITHER_ENABLE);
        //g2d.setRenderingHint(RenderingHints.KEY_FRACTIONALMETRICS, RenderingHints.VALUE_FRACTIONALMETRICS_ON);
        //g2d.setRenderingHint(RenderingHints.KEY_INTERPOLATION, RenderingHints.VALUE_INTERPOLATION_BILINEAR);
        //g2d.setRenderingHint(RenderingHints.KEY_RENDERING, RenderingHints.VALUE_RENDER_QUALITY);
        //g2d.setRenderingHint(RenderingHints.KEY_STROKE_CONTROL, RenderingHints.VALUE_STROKE_PURE);
        g2d.setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_ON);
        g2d.setFont(font);
        fm = g2d.getFontMetrics();
        g2d.setColor(Color.BLACK);
        for(ColoredChar x : coloredList){
        	
        	g2d.setColor(x.getBackgroundColor());
        	g2d.fillRect(
	        	x.getX() * fm.charWidth('X'),
	        	(int)Math.round(((x.getY() * fm.getHeight()) - fm.getHeight() * 0.8 )),
	        	fm.charWidth('X'),
	        	(int)Math.round(fm.getHeight() * 1)
	        );
	        
        	g2d.setColor(x.getForegroundColor());
        	g2d.drawString(Character.toString(x.getChar()), (float)x.getX() * fm.charWidth('X'), (float)(x.getY() * fm.getHeight()));
        }
        g2d.dispose();
        return img;
	}

	/** Convert the input text containing ansi formatting to a list of individual characters. 
	 * @throws IOException 
	 * */
	private List<ColoredChar> ansiFileToColoredChar() throws IOException {
		
		byte[] encoded = Files.readAllBytes(Paths.get(this.ansiInput.toURI()));
		String str= new String(encoded);
		List<String> ansiList= Splitter.on("\033").splitToList(str);

		List<ColoredChar> coloredList= new ArrayList<ColoredChar>();
		int x= 1;
		int y= 1;
		for(String xv : ansiList){
			if(xv.equals("[48;5;231m")){
				continue;
			}
			Color fgColor= this.ansiCodesToFgColor(xv);
			Color bgColor= this.ansiCodesToBgColor(xv);
			if(this.extractAnsiCodes(xv).size() != 0){
				// This string begins with ansi sequence and the color has been extracted.
				// So remove the ansi sequence at the beginnig and the end, we don't need them anymore
				xv= xv.replaceAll("^.+?m", "");
			}
			for(int i= 0; i < xv.length(); i++){
				char c= xv.charAt(i);
				if(c == '\n'){
					y += 1;
					x= 0;
				}
				coloredList.add(new ColoredChar(c, fgColor, bgColor, x, y));
				x += 1;
			}
		}
			
		return coloredList;
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
	
	/**String x is expected to start with the ansi escape sequence. From the ansi sequence
	 * extract the (last) color for the foreground color.
	 * If x doesn't start with the ansi escape return a default color (black?). 
	 * */
	protected Color ansiCodesToFgColor(String x){

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
	protected Color ansiCodesToBgColor(String x){

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
	
}
