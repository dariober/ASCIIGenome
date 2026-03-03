package tracks;

import exceptions.InvalidGenomicCoordsException;
import exceptions.InvalidRecordException;
import htsjdk.samtools.util.FileExtensions;
import htsjdk.tribble.index.tabix.TabixFormat;
import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.List;
import samTextViewer.GenomicCoords;
import samTextViewer.Utils;
import sortBgzipIndex.MakeTabixIndex;
import utils.CsvFormat;
import utils.FlexibleTabixReader;

public class TrackCSV extends TrackBedgraph {

  private final CsvFormat csvFormat;

  public TrackCSV(String filename, GenomicCoords gc, CsvFormat csvFormat)
      throws IOException,
          SQLException,
          ClassNotFoundException,
          InvalidGenomicCoordsException,
          InvalidRecordException {

    this.csvFormat = csvFormat;
    this.columnSeparator = csvFormat.getColumnSeparator();
    this.scoreColIdx = csvFormat.getScoreColIndex() + 1;
    this.setFilename(filename);
    this.setTrackFormat(TrackFormat.BEDGRAPH);

    if (!Utils.hasTabixIndex(filename)) {
      String suffix = new File(filename).getName();
      if (!suffix.endsWith(".gz")) {
        suffix += ".gz";
      }
      String tmpWorkFile =
          Utils.createTempFile(".asciigenome.", "." + suffix, true).getAbsolutePath();
      new File(tmpWorkFile + FileExtensions.TABIX_INDEX).deleteOnExit();
      this.setWorkFilename(tmpWorkFile);

      TabixFormat tabixFormat =
          new TabixFormat(
              this.csvFormat.isZeroBased() ? TabixFormat.ZERO_BASED : TabixFormat.GENERIC_FLAGS,
              this.csvFormat.getChromColIndex() + 1,
              this.csvFormat.getStartColIndex() + 1,
              this.csvFormat.getEndColIndex() + 1,
              this.csvFormat.getMetaCharacter(),
              this.csvFormat.getNumHeaderLinesToSkip());
      new MakeTabixIndex(
          filename, new File(this.getWorkFilename()), tabixFormat, csvFormat.getColumnSeparator());
      this.tabixReader = this.getTabixReader(this.getWorkFilename());
    } else { // This means the input is tabix indexed.
      this.setWorkFilename(filename);
      this.tabixReader = new FlexibleTabixReader(this.getWorkFilename());
    }
    this.setGc(gc);
  }

  @Override
  protected List<IntervalFeature> getFeaturesInInterval(String chrom, int from, int to)
      throws IOException, InvalidGenomicCoordsException {
    List<IntervalFeature> xFeatures = new ArrayList<>();
    TabixBigBedReader reader = this.getReader();
    TabixBigBedIterator qry = reader.query(chrom, from - 1, to);
    while (true) {
      String line = qry.next();
      if (line == null) {
        break;
      }
      IntervalFeature intervalFeature = new IntervalFeature(line, this.csvFormat);
      xFeatures.add(intervalFeature);
    }
    this.removeInvisibleFeatures(xFeatures);
    return xFeatures;
  }
}
