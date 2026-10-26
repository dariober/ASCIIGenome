package tracks;

import exceptions.InvalidGenomicCoordsException;
import exceptions.InvalidRecordException;
import org.checkerframework.checker.units.qual.C;
import samTextViewer.GenomicCoords;
import samTextViewer.Utils;
import utils.CsvFormat;
import utils.CsvPresets;

import java.io.IOException;
import java.sql.SQLException;

public class TrackGemma extends TrackBedgraph {

//  protected DataTransformation dataTransformation = DataTransformation.IDENTITY;
//  protected DataAggregationMethod dataAggregationMethod = DataAggregationMethod.MEAN;
//  protected int scoreColIdx = 12;

  public TrackGemma(String filename, GenomicCoords gc)
          throws SQLException,
          InvalidGenomicCoordsException,
          IOException,
          ClassNotFoundException,
          InvalidRecordException {
    this.csvFormat = CsvPresets.get(TrackFormat.GEMMA);
    this.scoreColIdx = this.csvFormat.getScoreColIndex() + 1;
    this.dataTransformation = DataTransformation.MINUS_LOG10;
    this.dataAggregationMethod = DataAggregationMethod.MAX;

    this.setFilename(filename);
    this.setTrackFormat(TrackFormat.GEMMA);

    if (Utils.hasTabixIndex(filename)) {
      this.setWorkFilename(filename);
    } else {
      this.sortAndIndex(filename);
    }
    this.tabixReader = this.getTabixReader(this.getWorkFilename());
    this.tabixReader.setColumnSeparator(this.csvFormat.getColumnSeparator());
    this.setGc(gc);
  }
}
