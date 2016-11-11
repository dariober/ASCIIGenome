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
		if(ansiColor == 31){ return Color.getHSBColor((float)(352.0/360.0), (float)0.76, (float)0.57); } // red
		if(ansiColor == 91){ return Color.getHSBColor((float)(356.0/360.0), (float)0.90, (float)0.87); } // light_red
		if(ansiColor == 32){ return Color.getHSBColor((float)(120.0/360.0), (float)0.71, (float)0.67); } // green
		if(ansiColor == 92){ return Color.getHSBColor((float)(117.0/360.0), (float)0.82, (float)0.85); } // light_green
		if(ansiColor == 33){ return Color.getHSBColor((float)(62.0/360.0), (float)0.41, (float)0.74); } // yellow
		if(ansiColor == 93){ return Color.getHSBColor((float)(62.0/360.0), (float)0.76, (float)0.91); } // light_yellow
		if(ansiColor == 34){ return Color.getHSBColor((float)(242.0/360.0), (float)0.68, (float)0.64); } // blue
		if(ansiColor == 94){ return Color.getHSBColor((float)(239.0/360.0), (float)0.88, (float)1.0); } // light_blue
		if(ansiColor == 35){ return Color.getHSBColor((float)(299.0/360.0), (float)0.68, (float)0.72); } // magenta
		if(ansiColor == 95){ return Color.getHSBColor((float)(299.0/360.0), (float)0.81, (float)0.91); } // light_magenta
		if(ansiColor == 36){ return Color.getHSBColor((float)(185.0/360.0), (float)0.71, (float)0.69); } // cyan
		if(ansiColor == 96){ return Color.getHSBColor((float)(185.0/360.0), (float)0.71, (float)0.90); } // light_cyan
		
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
		
		if(ansiColor == 147){ return Color.getHSBColor((float)(242.0/360.0), (float)0.39, 1); } // 147: Sort of lightblue
		if(ansiColor == 225){ return Color.getHSBColor((float)(299.0/360.0), (float)0.20, 1); } // 225: Sort of light pink

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
