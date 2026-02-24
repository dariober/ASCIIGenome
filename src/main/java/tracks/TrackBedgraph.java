package tracks;

import colouring.Config;
import colouring.ConfigKey;
import colouring.Xterm256;
import com.google.common.base.Joiner;
import com.google.common.base.Splitter;
import com.google.common.collect.Lists;
import exceptions.InvalidColourException;
import exceptions.InvalidGenomicCoordsException;
import exceptions.InvalidRecordException;
import htsjdk.samtools.util.FileExtensions;
import htsjdk.tribble.readers.TabixReader;
import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.List;
import org.apache.commons.lang3.NotImplementedException;
import org.apache.commons.lang3.StringUtils;
import org.apache.commons.lang3.math.NumberUtils;
import samTextViewer.GenomicCoords;
import samTextViewer.Utils;
import sortBgzipIndex.MakeTabixIndex;

public class TrackBedgraph extends TrackIntervalFeature {

  private DataTransformation dataTransformation = DataTransformation.IDENTITY;
  private DataAggregationMethod dataAggregationMethod = DataAggregationMethod.MEAN;

  public TrackBedgraph(String filename, GenomicCoords gc)
      throws ClassNotFoundException,
          IOException,
          InvalidGenomicCoordsException,
          InvalidRecordException,
          SQLException {
    this.scoreColIdx = 4;
    this.setFilename(filename);
    this.setTrackFormat(Utils.getFileTypeFromName(filename));

    if (!Utils.hasTabixIndex(filename)) {
      String suffix = new File(filename).getName();
      if (!suffix.endsWith(".gz")) {
        suffix += ".gz";
      }
      String tmpWorkFile =
          Utils.createTempFile(".asciigenome.", "." + suffix, true).getAbsolutePath();
      new File(tmpWorkFile + FileExtensions.TABIX_INDEX).deleteOnExit();
      this.setWorkFilename(tmpWorkFile);
      new MakeTabixIndex(
          filename,
          new File(this.getWorkFilename()),
          Utils.trackFormatToTabixFormat(this.getTrackFormat()));
      this.tabixReader = this.getTabixReader(this.getWorkFilename());
    } else { // This means the input is tabix indexed.
      this.setWorkFilename(filename);
      this.tabixReader = new TabixReader(this.getWorkFilename());
    }
    this.setGc(gc);

    this.setDataTransformation(this.dataTransformation);
  }

  public TrackBedgraph() {}

  /* ----------- METHODS ----------- */

  protected void setScoreColIdx(int scoreColIdx)
      throws ClassNotFoundException,
      IOException,
      InvalidGenomicCoordsException,
      InvalidRecordException,
      SQLException {
    if (this.scoreColIdx != scoreColIdx) {
      this.scoreColIdx = scoreColIdx;
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

  protected void setDataTransformation(DataTransformation dataTransformation) {
    this.dataTransformation = dataTransformation;
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

    List<Float> screenScores = this.aggregateScreenScores(screenWigLocInfoList);
    this.setScreenScores(screenScores);
  }

  private List<Float> aggregateScreenScores(List<ScreenWiggleLocusInfo> screenWigLocInfoList) {
    List<Float> screenScores = new ArrayList<>();
    for (ScreenWiggleLocusInfo x : screenWigLocInfoList) {
      if (this.getDataAggregationMethod() == DataAggregationMethod.MEAN) {
        screenScores.add(x.getMeanScore());
      } else if (this.getDataAggregationMethod() == DataAggregationMethod.MAX) {
        screenScores.add(x.getMax());
      } else if (this.getDataAggregationMethod() == DataAggregationMethod.MIN) {
        screenScores.add(x.getMin());
      } else {
        throw new NotImplementedException("Data '" + this.getDataAggregationMethod() + "' aggregation method not implemented");
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

    if (this.getScoreColIdx() < 4) {
      System.err.println(
          "Invalid index for bedgraph column of data value. Expected >=4. Got "
              + this.getScoreColIdx());
      this.scoreColIdx = 4;
      throw new InvalidRecordException();
    }
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

    Float[] range = Utils.range(this.getScreenScores());
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
