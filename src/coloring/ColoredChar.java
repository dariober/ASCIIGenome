package coloring;

import java.awt.Color;

public class ColoredChar {
	
	private char xchar;
	private Color foregroundColor= Color.BLACK; // The default color
	private Color backgroundColor= Color.WHITE; // The default color
	private int x= 0; // X and Y coordinates of this character. In screen coordinates.
	private int y= 0;
	
	public ColoredChar(char xchar, Color foregroundColor, Color backgroundColor, int x, int y){
		this.xchar= xchar;
		this.foregroundColor= foregroundColor;
		this.backgroundColor= backgroundColor;
		this.x= x;
		this.y= y;
		
	}
	
	protected Color getForegroundColor() {
		return this.foregroundColor;
		//return ansiColourToGraphicsColor(this.ansiFgColorCode);
	}
	protected Color getBackgroundColor() {
		return this.backgroundColor;
	}
	/** Input to ColorChar is the ansi color code. Here you translate this code 
	 * to the corresponding Color obj. 
	 * */
	protected static Color ansiForegroundColourToGraphicsColor(int ansiColor){
		
		/* Useful resources:
		Convert RGB to HSB:
		System.out.println( Arrays.toString(Color.RGBtoHSB(255, 225, 0, null)));
		*/
		
		if(ansiColor == 30){ return Color.getHSBColor(0, 0, 0); } // black
		if(ansiColor == 31){ return Color.getHSBColor(0, 1, 1); } // red
		if(ansiColor == 32){ return Color.getHSBColor((float)0.3333, 1, 1); } // green
		if(ansiColor == 33){ return Color.getHSBColor((float)0.14705883, 1, 1); } // yellow
		if(ansiColor == 34){ return Color.getHSBColor((float)0.6667, 1, 1); } // blue
		
		if(ansiColor == 37){ return Color.getHSBColor(0, 0, (float)0.7); } // light grey
		if(ansiColor == 90){ return Color.getHSBColor(0, 0, (float)0.4); } // dark grey
		if(ansiColor == 97){ return Color.getHSBColor(0, 0, 1); } // white

		return null;
	}

	protected static Color ansiBackgroundColourToGraphicsColor(int ansiColor){
		
		if(ansiColor == 40){ return Color.getHSBColor(0, 0, 0); } // black
		if(ansiColor == 41){ return Color.getHSBColor(0, 1, 1); } // red
		if(ansiColor == 42){ return Color.getHSBColor((float)0.3333, 1, 1); } // green
		if(ansiColor == 43){ return Color.getHSBColor((float)0.14705883, 1, 1); } // yellow
		if(ansiColor == 44){ return Color.getHSBColor((float)0.6667, 1, 1); } // blue
		
		if(ansiColor == 47){ return Color.getHSBColor(0, 0, (float)0.7); } // light grey
		if(ansiColor == 100){ return Color.getHSBColor(0, 0, (float)0.4); } // dark grey
		if(ansiColor == 107){ return Color.getHSBColor(0, 0, 1); } // white
		
		if(ansiColor == 147){ return Color.getHSBColor((float)0.5409, (float)0.2478, (float)0.9019); } // 147: Sort of lightblue
		if(ansiColor == 225){ return Color.getHSBColor((float)0.9748858, (float)0.28627452, (float)1.0); } // 225: Sort of light pink

		return null;
	}
		
	protected char getChar() {
		return this.xchar;
	}
	protected int getX() {
		return this.x;
	}
	protected int getY() {
		return this.y;
	}

	@Override
	public String toString(){
		return "char: " + this.xchar + " x: " + this.x + " y: " + this.y + " color: " + this.foregroundColor;
	}
}
