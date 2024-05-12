package tracks;

import coloring.Config;
import coloring.ConfigKey;
import coloring.Xterm256;
import com.google.common.base.Joiner;
import exceptions.InvalidColourException;
import exceptions.InvalidGenomicCoordsException;
import exceptions.InvalidRecordException;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.List;
import org.apache.commons.lang3.StringUtils;
import org.broad.igv.bbfile.BBFileReader;
import org.broad.igv.bbfile.BigWigIterator;
import org.broad.igv.bbfile.WigItem;
import org.broad.igv.tdf.TDFDataset;
import org.broad.igv.tdf.TDFGroup;
import org.broad.igv.tdf.TDFReader;
import org.broad.igv.tdf.TDFTile;
import org.broad.igv.util.ResourceLocator;
import samTextViewer.GenomicCoords;
import samTextViewer.Utils;

/** Process wiggle file formats. Mostly using IGV classes. bigBed, bigWig, */
public class TrackWiggles extends Track {

  private List<ScreenWiggleLocusInfo> screenWiggleLocusInfoList;
  private BBFileReader bigWigReader;

  /* C o n s t r u c t o r s */
  protected TrackWiggles(String filename, GenomicCoords gc, TrackFormat trackFormat)
      throws IOException,
          InvalidRecordException,
          InvalidGenomicCoordsException,
          ClassNotFoundException,
          SQLException {
    this.setFilename(filename);
    this.setWorkFilename(filename);
    this.setTrackFormat(trackFormat);

    if (this.getTrackFormat().equals(TrackFormat.BIGWIG)) {
      this.setTrackFormat(TrackFormat.BIGWIG);
      this.bigWigReader = new BBFileReader(this.getWorkFilename()); // or url for remote access.
      if (!this.bigWigReader.getBBFileHeader().isBigWig()) {
        throw new RuntimeException("Invalid file type " + this.getWorkFilename());
      }
    }
    this.setGc(gc);
  }

  public TrackWiggles(String filename, GenomicCoords gc)
      throws IOException,
          InvalidRecordException,
          InvalidGenomicCoordsException,
          ClassNotFoundException,
          SQLException {
    this(filename, gc, Utils.getFileTypeFromName(filename));
  }

  /*  M e t h o d s  */

  @Override
  public void close() {
    if (this.bigWigReader != null) {
      this.bigWigReader.close();
    }
  }

  @Override
  public void update()
      throws IOException,
          InvalidRecordException,
          InvalidGenomicCoordsException,
          ClassNotFoundException,
          SQLException {

    if (this.getTrackFormat().equals(TrackFormat.BIGWIG)) {

      this.bigWigToScores(this.bigWigReader);

    } else if (this.getTrackFormat().equals(TrackFormat.TDF)) {

      this.updateTDF();

    } else {
      throw new RuntimeException(
          "Invalid format not recognized for "
              + this.getWorkFilename()
              + ": "
              + this.getTrackFormat());
    }
  }

  private void updateTDF() throws InvalidGenomicCoordsException, IOException {

    this.screenWiggleLocusInfoList =
        this.tdfRangeToScreen(
            this.getWorkFilename(),
            this.getGc().getChrom(),
            this.getGc().getFrom(),
            this.getGc().getTo(),
            this.getGc().getMapping());

    List<Float> screenScores = new ArrayList<Float>();
    for (ScreenWiggleLocusInfo x : screenWiggleLocusInfoList) {
      screenScores.add((float) x.getMeanScore());
    }
    if (this.isRpm()) {
      screenScores = this.normalizeToRpm(screenScores);
    }
    this.setScreenScores(screenScores);
  }

  /**
   * Fetch data in tdf file in given range and puts it in a list of ScreenWiggleLocusInfo. a Adapted
   * from dumpRange. Really it should implement iterator.
   *
   * @param genomeToScreenMapping Typically from GenomicCoords.getMapping()
   * @author berald01
   */
  private List<ScreenWiggleLocusInfo> tdfRangeToScreen(
      String ibfFile,
      String chrom,
      int startLocation,
      int endLocation,
      List<Double> genomeToScreenMapping) {

    List<ScreenWiggleLocusInfo> screenWiggleLocusInfoList = new ArrayList<ScreenWiggleLocusInfo>();
    for (int i = 0; i < genomeToScreenMapping.size(); i++) {
      screenWiggleLocusInfoList.add(new ScreenWiggleLocusInfo());
    }

    TDFReader reader = TDFReader.getReader(ibfFile);

    for (String dsName : reader.getDatasetNames()) {

      String[] tokens = dsName.split("/");
      String chrName = tokens[1];
      if (!chrName.equals(chrom) || !dsName.contains("raw")) { // Not the right chrom or track
        continue;
      }

      TDFDataset ds = reader.getDataset(dsName);

      int tileWidth = ds.getTileWidth();
      int startTile = startLocation / tileWidth;
      int endTile = endLocation / tileWidth;

      for (int tileNumber = startTile; tileNumber <= endTile; tileNumber++) {
        TDFTile tile = reader.readTile(ds, tileNumber);
        if (tile == null) {
          // System.out.println("Null tile: " + dsName + " [" + tileNumber + "]");
        } else {
          int nTracks = reader.getTrackNames().length;
          if (nTracks > 1) {
            throw new RuntimeException("More than one track found in tdf file " + ibfFile);
          }
          int nBins = tile.getSize();
          if (nBins > 0) {
            for (int b = 0; b < nBins; b++) {
              int start = tile.getStartPosition(b);
              int end = tile.getEndPosition(b);
              if (start > endLocation) {
                break;
              }
              if (end >= startLocation) {
                int tileStartPos = tile.getStartPosition(b);
                float tileValue = tile.getValue(0, b);
                int idx =
                    Utils.getIndexOfclosestValue(
                        tileStartPos + 1,
                        genomeToScreenMapping); // Where should this position be mapped on screen?
                screenWiggleLocusInfoList.get(idx).increment(tileValue);
              }
            } // End process bins in this tile
          }
        } // End process this tile
      } // End iter tiles
    } // End iter datasets names
    reader.close();
    return screenWiggleLocusInfoList;
  }

  @Override
  protected void updateToRPM() {
    if (this.getTrackFormat().equals(TrackFormat.TDF)) {
      // Re-run update only for track types that can be converted to RPM
      try {
        this.update();
      } catch (ClassNotFoundException
          | IOException
          | InvalidRecordException
          | InvalidGenomicCoordsException
          | SQLException e) {
        e.printStackTrace();
      }
    }
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
              + Config.get256Color(ConfigKey.background)
              + ";38;5;"
              + Xterm256.colorNameToXterm256(this.getTitleColour())
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

    // xtitle= Utils.padEndMultiLine(xtitle, this.getGc().getUserWindowSize());
    return this.formatTitle(xtitle) + "\n";
  }

  /**
   * Populate object using bigWig data
   *
   * @throws IOException
   * @throws InvalidGenomicCoordsException
   */
  private void bigWigToScores(BBFileReader reader)
      throws InvalidGenomicCoordsException, IOException {

    // List of length equal to screen size. Each inner map contains info about the screen locus
    List<ScreenWiggleLocusInfo> screenWigLocInfoList = new ArrayList<ScreenWiggleLocusInfo>();
    for (int i = 0; i < getGc().getUserWindowSize(); i++) {
      screenWigLocInfoList.add(new ScreenWiggleLocusInfo());
    }

    BigWigIterator iter =
        reader.getBigWigIterator(
            getGc().getChrom(), getGc().getFrom(), getGc().getChrom(), getGc().getTo(), false);
    while (iter.hasNext()) {
      WigItem bw = iter.next();
      for (int i = bw.getStartBase(); i <= bw.getEndBase(); i++) {
        int idx =
            Utils.getIndexOfclosestValue(
                i, this.getGc().getMapping()); // Where should this position be mapped on screen?
        screenWigLocInfoList.get(idx).increment(bw.getWigValue());
      }
    }
    List<Float> screenScores = new ArrayList<Float>();
    for (ScreenWiggleLocusInfo x : screenWigLocInfoList) {
      screenScores.add((float) x.getMeanScore());
    }
    this.setScreenScores(screenScores);
  }

  private List<Float> normalizeToRpm(List<Float> screenScores) {
    ArrayList<Float> rpmed = new ArrayList<Float>();
    String x = this.getAttributesFromTDF("totalCount");
    if (x == null) {
      System.err.println("Warning: Cannot get total counts for " + this.getFilename());
      return screenScores;
    }
    Integer totalCount = Integer.parseInt(x);
    for (int i = 0; i < screenScores.size(); i++) {
      rpmed.add((float) (screenScores.get(i) / totalCount * 1000000.0));
    }
    return rpmed;
  }

  private String getAttributesFromTDF(String attr) {

    String path = this.getWorkFilename();

    try {
      ResourceLocator resourceLocator = new ResourceLocator(path);
      TDFReader reader = new TDFReader(resourceLocator);
      TDFGroup rootGroup = reader.getGroup("/");
      return rootGroup.getAttribute(attr);
    } catch (Exception e) {
      return null;
    }
  }

  @Override
  public ArrayList<String> getChromosomeNames() {

    if (this.getTrackFormat().equals(TrackFormat.TDF)) {

      ResourceLocator resourceLocator = new ResourceLocator(this.getWorkFilename());
      TDFReader reader = new TDFReader(resourceLocator);
      ArrayList<String> chroms = new ArrayList<String>(reader.getChromosomeNames());
      if (chroms.get(0).equals("All")) {
        chroms.remove(0);
      }
      return chroms;
      // chroms.addAll();
    }
    if (this.getTrackFormat().equals(TrackFormat.BIGWIG)) {
      return this.bigWigReader.getChromosomeNames();
    }
    return null;
  }

  @Override
  public String printLines() {
    return "";
  }

  @Override
  protected List<String> getRecordsAsStrings() {
    return new ArrayList<String>();
  }

  @Override
  public void setAwk(String awk)
      throws ClassNotFoundException,
          IOException,
          InvalidGenomicCoordsException,
          InvalidRecordException,
          SQLException {
    //
  }

  @Override
  public String getAwk() {
    return "";
  }

  @Override
  protected String getTitleForActiveFilters() {
    return "";
  }

  @Override
  public void reload()
      throws InvalidGenomicCoordsException,
          IOException,
          ClassNotFoundException,
          InvalidRecordException,
          SQLException {
    this.update();
  }

  @Override
  public void setFeatureName(String gtfAttributeForName) {}
}
