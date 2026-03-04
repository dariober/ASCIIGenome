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
import java.util.Map;
import org.apache.commons.lang3.NotImplementedException;
import org.apache.commons.lang3.StringUtils;
import samTextViewer.GenomicCoords;
import samTextViewer.Utils;
import utils.CsvFormat;

public class TrackBedgraph extends AbstractTrackFeature<QuantitativeFeature> {

  private DataTransformation dataTransformation = DataTransformation.IDENTITY;
  private DataAggregationMethod dataAggregationMethod = DataAggregationMethod.MEAN;
  protected int scoreColIdx = -1;

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
    if (csvFormat == null) {
      this.csvFormat = new CsvFormat(0, 1, 2, 3, true, 0, '#', '\t');
    } else {
      this.csvFormat = csvFormat;
    }
    this.scoreColIdx = this.csvFormat.getScoreColIndex() + 1;

    this.setFilename(filename);
    this.setTrackFormat(TrackFormat.BEDGRAPH);

    if (Utils.hasTabixIndex(filename)) {
      this.setWorkFilename(filename);
    } else {
      this.sortAndIndex(filename);
    }
    this.tabixReader = this.getTabixReader(this.getWorkFilename());
    this.tabixReader.setColumnSeparator(this.csvFormat.getColumnSeparator());
    this.setGc(gc);
  }

  public TrackBedgraph() {}

  /* ----------- METHODS ----------- */

  @Override
  protected List<QuantitativeFeature> getFeaturesInInterval(String chrom, int from, int to)
      throws IOException, InvalidGenomicCoordsException {
    List<QuantitativeFeature> xFeatures = new ArrayList<>();
    TabixBigBedIterator qry = this.getReader().query(chrom, from - 1, to);
    while (true) {
      String line = qry.next();
      if (line == null) {
        break;
      }
      xFeatures.add(new QuantitativeFeature(line, this.csvFormat));
    }
    this.removeInvisibleFeatures(xFeatures);
    return xFeatures;
  }

  protected int getScoreColIdx() {
    return this.scoreColIdx;
  }

  @Override
  protected QuantitativeFeature createFeature(String line) throws InvalidGenomicCoordsException {
    return new QuantitativeFeature(line, this.getCsvFormat());
  }

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
  private void bedGraphToScores(List<QuantitativeFeature> quantitativeFeatureList)
      throws IOException, InvalidGenomicCoordsException {

    List<ScreenWiggleLocusInfo> screenWigLocInfoList = new ArrayList<ScreenWiggleLocusInfo>();
    for (int i = 0; i < getGc().getUserWindowSize(); i++) {
      screenWigLocInfoList.add(new ScreenWiggleLocusInfo());
    }

    for (QuantitativeFeature ift : quantitativeFeatureList) {
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

    List<QuantitativeFeature> newFeatures = this.getFeaturesInInterval(
        this.getGc().getChrom(), this.getGc().getFrom(), this.getGc().getTo());
    this.setFeatureList(newFeatures);
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

  // This is a bad but for now let's get on with it
  // >---->
  @Override
  protected Map<String, List<QuantitativeFeature>> groupByGFFAttribute() {
    return Map.of();
  }

  @Override
  protected Map<String, List<QuantitativeFeature>> groupByGTFAttribute() {
    return Map.of();
  }

  @Override
  protected QuantitativeFeature collapseGFFTranscript(List<QuantitativeFeature> features, List<Double> mapToScreen) {
    return null;
  }
  // <----<
}
