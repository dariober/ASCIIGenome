package samTextViewer;

import coloring.Config;
import coloring.ConfigKey;
import coloring.Pdf;
import com.google.common.base.Joiner;
import com.google.common.base.Splitter;
import com.itextpdf.text.DocumentException;
import exceptions.InvalidColourException;
import exceptions.InvalidCommandLineException;
import exceptions.InvalidGenomicCoordsException;
import exceptions.InvalidRecordException;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import org.apache.commons.lang3.StringUtils;
import tracks.Track;
import tracks.TrackSet;

/** Process a TrackSet given the necessary elements */
public class TrackProcessor {

  private TrackSet trackSet;
  private boolean noFormat = false;
  private GenomicCoordsHistory genomicCoordsHistory;
  private String snapshotFile = null;
  private boolean appendToSnapshotFile = false;
  private boolean stripAnsi = true;
  private boolean showGruler = true;
  private boolean showCruler = true;

  /* C O N S T R U C T O R S */

  public TrackProcessor(TrackSet trackSet, GenomicCoordsHistory genomicCoordsHistory)
      throws IOException, InvalidRecordException {
    this.trackSet = trackSet;
    this.genomicCoordsHistory = genomicCoordsHistory;
  }

  /* M E T H O D S */

  public void iterateTracks()
      throws IOException,
          InvalidGenomicCoordsException,
          InvalidRecordException,
          ClassNotFoundException,
          SQLException,
          InvalidCommandLineException,
          DocumentException,
          InvalidColourException {

    final GenomicCoords currentGC = this.genomicCoordsHistory.current();

    StringBuilder outputString = new StringBuilder();

    currentGC.getChromIdeogram(20, this.noFormat);

    if (currentGC.getChromIdeogram(20, this.noFormat) != null) {
      outputString.append(currentGC.getChromIdeogram(20, this.noFormat) + "\n");
    }

    // Update tracks to new genomic coords
    for (Track track : trackSet.getTrackList()) {
      if (!track.getGc().equalCoordsAndWindowSize(currentGC)
          && track.getyMaxLines() > 0
          && !track.isHideTrack()) {
        track.setGc(currentGC);
      }
    }
    // Set new y limits as required. This step has to come after the positioning to new coordinates
    // because
    // we may need to autoscale to global min or max.
    this.getTrackSet().setAutoYLimits();

    // Visualise as required
    for (Track track : trackSet.getTrackList()) {

      track.setNoFormat(this.noFormat);
      if (track.getyMaxLines() > 0 && !track.isHideTrack()) {
        outputString.append(track.concatTitleAndTrack() + "\n");
        outputString.append(track.printLines());
      }
    }

    // Ruler and sequence
    // ------------------
    outputString.append(currentGC.printableRefSeq(noFormat));
    if (this.isShowGruler()) {
      outputString.append(currentGC.printableGenomicRuler(10, noFormat) + "\n");
    }
    if (this.isShowCruler()) {
      outputString.append(currentGC.printablePercentRuler(10, noFormat) + "\n");
    }

    // Position, memory, etc
    // ---------------------
    String footer = this.getFooter(currentGC);
    if (!noFormat) {
      outputString.append("\033[48;5;");
      outputString.append(Config.get256Color(ConfigKey.background));
      outputString.append(";38;5;");
      outputString.append(Config.get256Color(ConfigKey.footer));
      outputString.append("m");
      outputString.append(footer);
      outputString.append("\033[38;5;");
      outputString.append(Config.get256Color(ConfigKey.foreground));
      outputString.append("m");
    } else {
      outputString.append(footer);
    }

    String printable = this.fillUpLines(outputString.toString());

    // Print to screen
    System.out.println(printable);

    // Optionally save to file
    // -----------------------
    if (this.snapshotFile != null && this.snapshotFile.endsWith(".pdf")) {
      (new Pdf(printable)).convert(new File(this.snapshotFile), 10, this.appendToSnapshotFile);

    } else if (this.snapshotFile != null) {
      BufferedWriter wr =
          new BufferedWriter(
              new FileWriter(new File(this.snapshotFile), this.appendToSnapshotFile));
      if (this.getStripAnsi()) {
        wr.write(Utils.stripAnsiCodes(printable));
      } else {
        wr.write(printable);
      }
      wr.write("\n-------8<-------------[ cut here ]----------------------\n\n");
      wr.close();
    }
    this.snapshotFile = null;
  }

  /**
   * String screenshot is the string almost ready to be printed to screen or file. What is missing
   * is to fill up lines with whitespaces until the end of the screen so that terminals like tmux do
   * not show lines of mixed colours.
   *
   * @throws IOException
   * @throws InvalidGenomicCoordsException
   */
  private String fillUpLines(String screenshot) throws InvalidGenomicCoordsException, IOException {
    StringBuilder full = new StringBuilder();
    Iterator<String> iter = Splitter.on("\n").split(screenshot).iterator();
    int screenSize = this.getWindowSize();
    while (iter.hasNext()) {
      String line = iter.next();
      int npad = screenSize - Utils.stripAnsiCodes(line).length();
      if (npad < 0) {
        npad = 0; // Nothing to be done.
      }
      String filler = StringUtils.repeat(' ', npad);
      full.append(line);
      full.append(filler);
      if (iter.hasNext()) {
        full.append("\n");
      }
    }
    return full.toString();
  }

  private String getFooter(GenomicCoords currentGC)
      throws InvalidGenomicCoordsException, IOException {

    List<String> footList = new ArrayList<String>();
    footList.add(currentGC.toString());
    footList.add(Math.rint(currentGC.getBpPerScreenColumn() * 10d) / 10d + " bp/char");
    String footer = Joiner.on("; ").join(footList);

    // Add midpoint marker
    int mid = (int) Math.rint(this.getWindowSize() / 2.0);
    int npad = mid - footer.length();
    if (npad > 1) {
      String pad = StringUtils.repeat(" ", npad - 1);
      footer += pad;
      footer += "/\\";
    }
    return footer;
  }

  protected TrackSet getTrackSet() {
    return trackSet;
  }

  protected void setTrackSet(TrackSet trackSet) {
    this.trackSet = trackSet;
  }

  protected boolean isNoFormat() {
    return noFormat;
  }

  protected void setNoFormat(boolean noFormat) {
    this.noFormat = noFormat;
  }

  protected GenomicCoordsHistory getGenomicCoordsHistory() {
    return genomicCoordsHistory;
  }

  protected void setGenomicCoordsHistory(GenomicCoordsHistory genomicCoordsHistory) {
    this.genomicCoordsHistory = genomicCoordsHistory;
  }

  protected String getSnapshotFile() {
    return snapshotFile;
  }

  protected void setSnapshotFile(String snapshotFile) {
    this.snapshotFile = snapshotFile;
  }

  protected int getWindowSize() throws InvalidGenomicCoordsException, IOException {
    return this.genomicCoordsHistory.current().getUserWindowSize();
  }

  protected boolean isAppendToSnapshotFile() {
    return appendToSnapshotFile;
  }

  protected void setAppendToSnapshotFile(boolean appendToSnapshotFile) {
    this.appendToSnapshotFile = appendToSnapshotFile;
  }

  protected void setStripAnsi(boolean stripAnsi) {
    this.stripAnsi = stripAnsi;
  }

  protected boolean getStripAnsi() {
    return this.stripAnsi;
  }

  protected boolean isShowGruler() {
    return showGruler;
  }

  protected void setShowGruler(boolean showGruler) {
    this.showGruler = showGruler;
  }

  protected boolean isShowCruler() {
    return showCruler;
  }

  protected void setShowCruler(boolean showCruler) {
    this.showCruler = showCruler;
  }
}
