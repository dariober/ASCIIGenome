package tracks;

import colouring.Xterm256;
import com.google.common.collect.Lists;
import exceptions.InvalidColourException;
import exceptions.InvalidCommandLineException;
import exceptions.InvalidGenomicCoordsException;
import exceptions.InvalidRecordException;
import htsjdk.samtools.util.FileExtensions;
import htsjdk.tribble.readers.TabixReader;
import htsjdk.variant.vcf.VCFCodec;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.regex.Pattern;
import org.apache.commons.lang3.StringUtils;
import org.broad.igv.bbfile.BBFileReader;
import samTextViewer.GenomicCoords;
import samTextViewer.Utils;
import sortBgzipIndex.MakeTabixIndex;

public class TrackIntervalFeature extends AbstractTrackIntervalFeature<IntervalFeature> {

  /** For GTF/GFF data: Use this attribute to get the feature names */
  protected int scoreColIdx = -1;
  private List<Argument> colourForRegex = null;
  private String gtfAttributeForName = null;
  private int bedFieldForName = 3; // 0-based!

  /* C o n s t r u c t o r */

  public TrackIntervalFeature(final String filename, GenomicCoords gc)
      throws IOException,
          InvalidGenomicCoordsException,
          ClassNotFoundException,
          InvalidRecordException,
          SQLException {

    this.setFilename(filename);
    this.setTrackFormat(Utils.getFileTypeFromName(filename));

    if (this.getTrackFormat().equals(TrackFormat.BIGBED)
        || this.getTrackFormat().equals(TrackFormat.BED)) {
      this.scoreColIdx =
          5; // Don't use setScoreColIdx because it call update and the GenomicCoordinates object is
      // null
    } else if (this.getTrackFormat().equals(TrackFormat.BEDGRAPH)) {
      this.scoreColIdx = 4;
    }
    if (this.getTrackFormat().equals(TrackFormat.BIGBED)) {
      this.bigBedReader = new BBFileReader(filename); // or url for remote access.
      if (!this.bigBedReader.getBBFileHeader().isBigBed()) {
        throw new RuntimeException("File " + filename + " is not bigBed.");
      }
      this.setWorkFilename(filename);
    } else if (!Utils.hasTabixIndex(filename)) {
      // Tabix index not found for this file. Sort and index input to tmp.
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
  }

  protected TrackIntervalFeature(GenomicCoords gc) {}

  public TrackIntervalFeature() {}

  /* M e t h o d s */
  @Override
  protected IntervalFeature mergeFeatures(
          IntervalFeature a,
          IntervalFeature b,
          boolean screenCoords)
          throws InvalidGenomicCoordsException, InvalidColourException {

    IntervalFeature x =
            new IntervalFeature(
                    a.getChrom() + "\t" +
                            (a.getFrom() - 1) + "\t" +
                            b.getTo(),
                    TrackFormat.BED,
                    null,
                    -1
            );

    if (screenCoords) {
      x.setScreenFrom(a.getScreenFrom());
      x.setScreenTo(b.getScreenTo());
    }

    if (a.getStrand() == b.getStrand()) {
      x.setStrand(a.getStrand());
    }

    return x;
  }

  @Override
  protected IntervalFeature createFeature(String line)
          throws InvalidGenomicCoordsException {
    return new IntervalFeature(
            line,
            getTrackFormat(),
            null,
            getScoreColIdx()
    );
  }

  protected List<IntervalFeature> getFeaturesInInterval(String chrom, int from, int to)
          throws IOException, InvalidGenomicCoordsException {
    List<IntervalFeature> xFeatures = new ArrayList<>();
    TabixBigBedIterator qry = this.getReader().query(chrom, from - 1, to);
    while (true) {
      String q = qry.next();
      if (q == null) {
        break;
      }
      IntervalFeature intervalFeature =
              new IntervalFeature(q, this.getTrackFormat(), null, this.getScoreColIdx());
      xFeatures.add(intervalFeature);
    }
    this.removeInvisibleFeatures(xFeatures);
    return xFeatures;
  }

}
