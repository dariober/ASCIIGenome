package tracks;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import exceptions.InvalidGenomicCoordsException;
import org.junit.Test;
import utils.CsvFormat;

public class QuantitativeFeatureTest {

  @Test
  public void canParseScorecolumn() throws InvalidGenomicCoordsException {

    String line = "chr1 0 1 9 8 FOO 10".replaceAll(" ", "\t");

    CsvFormat csv = new CsvFormat(0, 1, 2, 4, true, 0, '#', '\t');

    // Default column indexes for scores
    QuantitativeFeature ift = new QuantitativeFeature(line, csv);
    assertEquals(8, ift.getScore(), 0.001);

    csv.setScoreColIndex(6);
    ift = new QuantitativeFeature(line, csv);
    assertEquals(10, ift.getScore(), 0.001);

    csv.setScoreColIndex(5);
    ift = new QuantitativeFeature(line, csv);
    assertEquals(Float.NaN, ift.getScore(), 0.001);
  }

  @Test
  public void canConstructFromCsv() throws InvalidGenomicCoordsException {
    String line = "10,chr1,9.9";
    CsvFormat csv = new CsvFormat(1, 0, -1, 2, false, 0, '#', ',');
    QuantitativeFeature x = new QuantitativeFeature(line, csv);
    assertEquals("chr1", x.getChrom());
    assertEquals(10, x.getFrom());
    assertEquals(10, x.getTo());
    assertEquals(9.9, x.getScore(), 0.0001);

    line = "9|chr1|20|foo";
    csv = new CsvFormat(1, 0, 2, 3, true, 0, '#', '|');
    x = new QuantitativeFeature(line, csv);
    assertEquals("chr1", x.getChrom());
    assertEquals(10, x.getFrom());
    assertEquals(20, x.getTo());
    assertEquals(Double.NaN, x.getScore(), 0.0001);

    line = "chr1\t0\t1";
    csv = new CsvFormat(0, 1, 2, -1, true, 0, '#', '\t');
    boolean pass = false;
    try {
      new QuantitativeFeature(line, csv);
    } catch (RuntimeException e) {
      assertEquals("Invalid index for score column: 0", e.getMessage().trim());
      pass = true;
    }
    assertTrue(pass);
  }
}
