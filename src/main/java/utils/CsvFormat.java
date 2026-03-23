package utils;

import java.lang.reflect.Field;

public class CsvFormat {

  private final int chromColIndex;
  private final int startColIndex;
  private final int endColIndex;
  private int scoreColIndex;
  private final Boolean isZeroBased;
  private int numHeaderLinesToSkip;
  private final char metaCharacter;
  private char columnSeparator;

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

  public void setNumHeaderLinesToSkip(int numHeaderLinesToSkip) {
    this.numHeaderLinesToSkip = numHeaderLinesToSkip;
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

  public void setColumnSeparator(char columnSeparator) {
    this.columnSeparator = columnSeparator;
  }

  @Override
  public String toString() {
    StringBuilder sb = new StringBuilder();
    sb.append(getClass().getSimpleName()).append(" {");

    Field[] fields = getClass().getDeclaredFields(); // Get all fields, including private
    for (Field field : fields) {
      try {
        field.setAccessible(true); // To access private fields
        sb.append(field.getName()) // Field name
            .append("=")
            .append(field.get(this)) // Field value
            .append(", ");
      } catch (IllegalAccessException e) {
        e.printStackTrace(); // Handle the exception if a field is not accessible
      }
    }
    // Remove the trailing comma and space, if there are any fields
    if (fields.length > 0) {
      sb.setLength(sb.length() - 2);
    }
    sb.append("}");
    return sb.toString();
  }
}
