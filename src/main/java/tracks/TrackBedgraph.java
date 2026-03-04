package tracks;

import colouring.Config;
import colouring.ConfigKey;
import colouring.Xterm256;
import com.google.common.base.Joiner;
import exceptions.InvalidColourException;
import exceptions.InvalidGenomicCoordsException;
import exceptions.InvalidRecordException;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.List;
import org.apache.commons.lang3.NotImplementedException;
import org.apache.commons.lang3.StringUtils;
import samTextViewer.GenomicCoords;
import samTextViewer.Utils;
import utils.CsvFormat;

public class TrackBedgraph extends TrackIntervalFeature {

  private DataTransformation dataTransformation = DataTransformation.IDENTITY;
  private DataAggregationMethod dataAggregationMethod = DataAggregationMethod.MEAN;

  public TrackBedgraph(String filename, GenomicCoords gc)
      throws SQLException,
          InvalidGenomicCoordsException,
          IOException,
          ClassNotFoundException,
          InvalidRecordException {
    this(filename, gc, null);
  }

  public TrackBedgraph(String filename, GenomicCoords gc, CsvFormat csvFormat)
      throws ClassNotFoundException,
          IOException,
          InvalidGenomicCoordsException,
          InvalidRecordException,
          SQLException {
    this.csvFormat = csvFormat;
    this.scoreColIdx = csvFormat == null ? 4 : csvFormat.getScoreColIndex() + 1;

    this.setFilename(filename);
    this.setTrackFormat(TrackFormat.BEDGRAPH);

    if (Utils.hasTabixIndex(filename)) {
      this.setWorkFilename(filename);
    } else {
      this.sortAndIndex(filename);
    }
    this.tabixReader = this.getTabixReader(this.getWorkFilename());
    this.tabixReader.setColumnSeparator(csvFormat == null ? '\t' : csvFormat.getColumnSeparator());
    this.setGc(gc);
  }

  public TrackBedgraph() {}

  /* ----------- METHODS ----------- */

  /** NB: index here is 1-based
   * */
  protected void setScoreColIdx(int scoreColIdx)
      throws ClassNotFoundException,
          IOException,
          InvalidGenomicCoordsException,
          InvalidRecordException,
          SQLException {
    if (this.scoreColIdx != scoreColIdx) {
      this.scoreColIdx = scoreColIdx;
      if (this.csvFormat != null) {
        this.csvFormat.setScoreColIndex(scoreColIdx - 1); // CsvFormat is 0-based
      }
      this.update();
    }
  }

  protected void setDataAggregationMethod(DataAggregationMethod dataAggregationMethod)
      throws ClassNotFoundException,
          IOException,
          InvalidGenomicCoordsException,
          InvalidRecordException,
          SQLException {
    if (this.dataAggregationMethod != dataAggregationMethod) {
      this.dataAggregationMethod = dataAggregationMethod;
      this.update();
    }
  }

  protected DataAggregationMethod getDataAggregationMethod() {
    return this.dataAggregationMethod;
  }

  protected void setDataTransformation(DataTransformation dataTransformation)
      throws SQLException,
          InvalidGenomicCoordsException,
          IOException,
          InvalidRecordException,
          ClassNotFoundException {
    if (this.dataTransformation != dataTransformation) {
      this.dataTransformation = dataTransformation;
      this.update();
    }
  }

  protected DataTransformation getDataTransformation() {
    return this.dataTransformation;
  }

  /** Get values for bedgraph */
  private void bedGraphToScores(List<IntervalFeature> intervalFeatureList)
      throws IOException, InvalidGenomicCoordsException {

    List<ScreenWiggleLocusInfo> screenWigLocInfoList = new ArrayList<ScreenWiggleLocusInfo>();
    for (int i = 0; i < getGc().getUserWindowSize(); i++) {
      screenWigLocInfoList.add(new ScreenWiggleLocusInfo());
    }

    for (IntervalFeature ift : intervalFeatureList) {
      ift.mapToScreen(this.getGc().getMapping());
      for (int i = ift.getScreenFrom(); i <= ift.getScreenTo(); i++) {
        screenWigLocInfoList.get(i).increment(ift.getScore(), this.getDataTransformation());
      }
    }

    List<Double> screenScores = this.aggregateScreenScores(screenWigLocInfoList);
    this.setScreenScores(screenScores);
  }

  private List<Double> aggregateScreenScores(List<ScreenWiggleLocusInfo> screenWigLocInfoList) {
    List<Double> screenScores = new ArrayList<>();
    for (ScreenWiggleLocusInfo x : screenWigLocInfoList) {
      if (this.getDataAggregationMethod() == DataAggregationMethod.MEAN) {
        screenScores.add(x.getMeanScore());
      } else if (this.getDataAggregationMethod() == DataAggregationMethod.MAX) {
        screenScores.add(x.getMax());
      } else if (this.getDataAggregationMethod() == DataAggregationMethod.MIN) {
        screenScores.add(x.getMin());
      } else {
        throw new NotImplementedException(
            "Data '" + this.getDataAggregationMethod() + "' aggregation method not implemented");
      }
    }
    return screenScores;
  }

  @Override
  public void update()
      throws IOException,
          InvalidRecordException,
          InvalidGenomicCoordsException,
          ClassNotFoundException,
          SQLException {

    this.setFeatureList(
        this.getFeaturesInInterval(
            this.getGc().getChrom(), this.getGc().getFrom(), this.getGc().getTo()));
    this.setDataTransformation(this.dataTransformation);
    this.bedGraphToScores(this.getFeatureList());
  }

  @Override
  public String printToScreen() throws InvalidColourException {

    if (this.getyMaxLines() == 0) {
      return "";
    }

    TextProfile textProfile =
        new TextProfile(
            this.getScreenScores(), this.getyMaxLines(), this.getYLimitMin(), this.getYLimitMax());

    ArrayList<String> lineStrings = new ArrayList<String>();
    for (int i = (textProfile.getProfile().size() - 1); i >= 0; i--) {
      List<String> xl = textProfile.getProfile().get(i);
      lineStrings.add(StringUtils.join(xl, ""));
    }

    String printable = Joiner.on("\n").join(lineStrings);
    if (!this.isNoFormat()) {
      new Xterm256();
      printable =
          "\033[48;5;"
              + Config.get256Colour(ConfigKey.background)
              + ";38;5;"
              + Xterm256.colourNameToXterm256(this.getTitleColour())
              + "m"
              + printable;
    }
    return printable;
  }

  @Override
  public String getTitle()
      throws InvalidColourException, InvalidGenomicCoordsException, IOException {

    if (this.isHideTitle()) {
      return "";
    }

    Double[] range = Utils.range(this.getScreenScores());
    String[] rounded = Utils.roundToSignificantDigits(range[0], range[1], 2);

    String ymin = this.getYLimitMin().isNaN() ? "auto" : this.getYLimitMin().toString();
    String ymax = this.getYLimitMax().isNaN() ? "auto" : this.getYLimitMax().toString();

    String xtitle =
        this.getTrackTag()
            + "; ylim["
            + ymin
            + " "
            + ymax
            + "]"
            + "; range["
            + rounded[0]
            + " "
            + rounded[1]
            + "]";

    String filters = this.getTitleForActiveFilters();
    xtitle += filters;

    return this.formatTitle(xtitle) + "\n";
  }
}
