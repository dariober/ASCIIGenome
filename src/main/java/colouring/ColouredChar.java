package colouring;

import java.awt.Color;

public class ColouredChar {

  private char xchar;
  private Color foregroundColour = Color.BLACK; // The default colour
  private Color backgroundColour = Color.WHITE; // The default colour
  private int x = 0; // X and Y coordinates of this character. In screen coordinates.
  private int y = 0;

  public ColouredChar(char xchar, Color foregroundColour, Color backgroundColour, int x, int y) {
    this.xchar = xchar;
    this.foregroundColour = foregroundColour;
    this.backgroundColour = backgroundColour;
    this.x = x;
    this.y = y;
  }

  protected Color getForegroundColour() {
    return this.foregroundColour;
    // return ansiColourToGraphicsColour(this.ansiFgColourCode);
  }

  protected Color getBackgroundColour() {
    return this.backgroundColour;
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
  public String toString() {
    return "char: "
        + this.xchar
        + " x: "
        + this.x
        + " y: "
        + this.y
        + " colour: "
        + this.foregroundColour;
  }
}
