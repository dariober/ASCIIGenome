package utils;

public class CsvFormat {

  private final int chromColIndex;
  private final int startColIndex;
  private final int endColIndex;
  private int scoreColIndex;
  private final Boolean isZeroBased;
  private final int numHeaderLinesToSkip;
  private final char metaCharacter;
  private final char columnSeparator;

  /** NB: Column indexes are 0-based! */
  public CsvFormat(
      int chromColIndex,
      int startColIndex,
      int endColIndex,
      int scoreColIndex,
      Boolean isZeroBased,
      int numHeaderLinesToSkip,
      char metaCharacter,
      char columnSeparator) {

    this.chromColIndex = chromColIndex;
    this.endColIndex = endColIndex;
    this.startColIndex = startColIndex;
    this.scoreColIndex = scoreColIndex;
    this.isZeroBased = isZeroBased;
    this.numHeaderLinesToSkip = numHeaderLinesToSkip;
    this.metaCharacter = metaCharacter;
    this.columnSeparator = columnSeparator;
  }

  public int getChromColIndex() {
    return chromColIndex;
  }

  public int getStartColIndex() {
    return startColIndex;
  }

  public int getEndColIndex() {
    if (endColIndex < 0) {
      return startColIndex;
    }
    return endColIndex;
  }

  public Boolean isZeroBased() {
    return isZeroBased;
  }

  public char getMetaCharacter() {
    return metaCharacter;
  }

  public int getNumHeaderLinesToSkip() {
    return numHeaderLinesToSkip;
  }

  public int getScoreColIndex() {
    return scoreColIndex;
  }

  public void setScoreColIndex(int scoreColIndex) {
    this.scoreColIndex = scoreColIndex;
  }

  public char getColumnSeparator() {
    return columnSeparator;
  }
}
