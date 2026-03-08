package tracks;

import com.google.common.base.Splitter;
import java.util.List;
import utils.CsvFormat;

public class QuantitativeFeature extends IntervalFeature {

  private double score;

  /* C o n s t r u c t o r s */

  public QuantitativeFeature(String line, CsvFormat csvFormat) {
    if (csvFormat == null) {
      csvFormat = new CsvFormat(0, 1, 2, 3, true, 0, '#', '\t');
    }

    this.setRaw(line);
    List<String> row = Splitter.on(csvFormat.getColumnSeparator()).splitToList(line);
    if (row.size() < 2) {
      throw new RuntimeException("Invalid line:\n" + row);
    }

    this.setChrom(row.get(csvFormat.getChromColIndex()).trim());
    this.setFrom(Integer.parseInt(row.get(csvFormat.getStartColIndex())));
    int from = this.getFrom();
    from += csvFormat.isZeroBased() ? 1 : 0;
    this.setFrom(from);
    this.setTo(
        csvFormat.getEndColIndex() > 0
            ? Integer.parseInt(row.get(csvFormat.getEndColIndex()))
            : this.getFrom());
    try {
      this.score = Double.parseDouble(row.get(csvFormat.getScoreColIndex()));
    } catch (NumberFormatException e) {
      this.score = Double.NaN;
    } catch (IndexOutOfBoundsException e) {
      throw new RuntimeException(
          "Invalid index for score column: " + (csvFormat.getScoreColIndex() + 1));
    }
  }

  /* M e t h o d s */

  public double getScore() {
    return score;
  }
}
