package utils;

import htsjdk.tribble.Feature;

public class GtfLine implements Feature {
  private String tokens[];

  public GtfLine(String tokens[]) {
    this.tokens = tokens;
  }

  @Override
  @Deprecated
  public final String getChr() {
    return getContig();
  }

  @Override
  public String getContig() {
    return tokens[0];
  }

  @Override
  public int getStart() {
    return Integer.parseInt(tokens[3]);
  }

  @Override
  public int getEnd() {
    return (tokens.length < 3 ? getStart() : Integer.parseInt(tokens[4]));
  }

  public String get(int index) {
    return (index < tokens.length ? tokens[index] : null);
  }

  public String join(final CharSequence delimiter) {
    StringBuilder sb = new StringBuilder();
    for (String x : this.tokens) {
      sb.append(x);
      sb.append(delimiter);
    }
    return sb.toString().replace("\t$", "");
    // return String.join(delimiter, this.tokens);
  }

  public String join() {
    return join("\t");
  }

  public int getColumnCount() {
    return tokens.length;
  }

  public static boolean isBedHeader(String line) {
    return line.startsWith("#") || line.startsWith("track") || line.startsWith("browser");
  }
}
