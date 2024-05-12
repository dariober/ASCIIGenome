package tracks;

import coloring.Config;
import coloring.ConfigKey;
import coloring.Xterm256;
import com.google.common.base.Joiner;
import com.google.common.base.Splitter;
import com.google.common.base.Strings;
import exceptions.InvalidColourException;
import java.util.ArrayList;
import java.util.List;

class TrackHeader {
  private String headerText = null;
  private String color = null;
  private boolean isBold = true;
  // private HeaderAlignment headerAlignment = HeaderAlignment.CENTER;
  private double headerAlignmentPct = 0.5;
  private int terminalWidth = 0;

  public TrackHeader(String header) {
    this.setHeaderText(header);
  }

  /*   M E T H O D S   */

  protected String format(boolean isNoFormat) throws InvalidColourException {

    if (this.headerText == null) {
      return "";
    }

    String text = this.padHeaderText() + "\n";

    if (isNoFormat) {
      return text;
    } else {
      int colourCode = Config.get256Color(ConfigKey.foreground);
      if (this.color != null) {
        new Xterm256();
        colourCode = Xterm256.colorNameToXterm256(this.color);
      }
      String bold = "";
      if (this.isBold) {
        bold = "1;";
      }
      return "\033["
          + bold
          + "48;5;"
          + Config.get256Color(ConfigKey.background)
          + ";38;5;"
          + colourCode
          + "m"
          + text
          + "\033[0m";
    }
  }

  private String padHeaderText() {
    List<String> lines = Splitter.on("\n").splitToList(this.headerText);
    List<String> padded = new ArrayList<String>();
    for (String x : lines) {
      int gap = (int) this.terminalWidth - x.length();
      if (gap > 0) {
        int leftPad = (int) Math.floor(this.getHeaderAlignmentPct() * gap);
        x = Strings.repeat(" ", leftPad) + x;
      }
      padded.add(x);
    }
    return Joiner.on("\n").join(padded);
  }

  protected boolean isBold() {
    return isBold;
  }

  protected void setBold(boolean isBold) {
    this.isBold = isBold;
  }

  protected void setHeaderText(String header) {
    if (header == null) {
      this.headerText = header;
    } else {
      String current = this.headerText;
      if (current != null) {
        header = header.replace("{-}", current);
      }
      this.headerText = header;
    }
  }

  protected void setColor(String color) {
    this.color = color;
  }

  protected int getTerminalWidth() {
    return terminalWidth;
  }

  protected void setTerminalWidth(int terminalWidth) {
    this.terminalWidth = terminalWidth;
  }

  protected double getHeaderAlignmentPct() {
    return headerAlignmentPct;
  }

  protected void setHeaderAlignmentPct(double headerAlignmentPct) {
    if (headerAlignmentPct < 0) {
      headerAlignmentPct = 0;
    }
    if (headerAlignmentPct > 1) {
      headerAlignmentPct = 1;
    }
    this.headerAlignmentPct = headerAlignmentPct;
  }
}
