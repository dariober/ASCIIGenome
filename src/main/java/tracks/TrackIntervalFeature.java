package tracks;

import coloring.Xterm256;
import com.google.common.base.Joiner;
import com.google.common.base.Splitter;
import com.google.common.collect.Lists;
import exceptions.InvalidColourException;
import exceptions.InvalidCommandLineException;
import exceptions.InvalidGenomicCoordsException;
import exceptions.InvalidRecordException;
import htsjdk.samtools.util.FileExtensions;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.tribble.readers.TabixReader;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import java.io.File;
import java.io.IOException;
import java.net.URL;
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

public class TrackIntervalFeature extends Track {

  protected List<IntervalFeature> intervalFeatureList = new ArrayList<IntervalFeature>();

  /** For GTF/GFF data: Use this attribute to get the feature names */
  protected TabixReader tabixReader; // Leave *protected* for TrackBookmark to work

  private BBFileReader bigBedReader;
  protected int scoreColIdx = -1;
  private List<Argument> colorForRegex = null;
  private VCFCodec vcfCodec;
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

  protected TrackIntervalFeature() {}

  /* M e t h o d s */

  @Override
  public void close() {
    if (this.tabixReader != null) {
      this.tabixReader.close();
    }
    if (this.bigBedReader != null) {
      this.bigBedReader.close();
    }
  }

  @Override
  /**
   * Collect features mapping to the current genomic coordinates and update the list of interval
   * features for the track. Also update the mapping of features to the terminal screen. This method
   * should be called only from within methods that change the features being displayed. E.g.
   * setGc(), which changes the coordinates, or setHideRegex() & setShowRegex() which change the
   * visible features. update() should not change anything other than the list of features and the
   * mapping.
   */
  public void update()
      throws IOException,
          InvalidGenomicCoordsException,
          ClassNotFoundException,
          InvalidRecordException,
          SQLException {

    this.intervalFeatureList =
        this.getFeaturesInInterval(
            this.getGc().getChrom(), this.getGc().getFrom(), this.getGc().getTo());

    for (IntervalFeature ift : this.intervalFeatureList) {
      ift.mapToScreen(this.getGc().getMapping());
    }
  }

  protected List<IntervalFeature> getFeaturesInInterval(String chrom, int from, int to)
      throws IOException, InvalidGenomicCoordsException {

    if (from < 1) {
      System.err.println("from < 1: " + from + "; resetting to 1.");
      from = 1;
    }

    if (from > to || to < 1) {
      System.err.println(
          "Invalid coordinates: from: "
              + from
              + "; to: "
              + to
              + "; Resetting to initial 1-"
              + Integer.MAX_VALUE);
      from = 1;
      to = Integer.MAX_VALUE;
      throw new InvalidGenomicCoordsException();
    }

    List<IntervalFeature> xFeatures = new ArrayList<IntervalFeature>();

    if (this.getTrackFormat().equals(TrackFormat.VCF)) {
      xFeatures = getFeaturesInVCFInterval(chrom, from, to);
    } else {
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
    }
    this.removeInvisibleFeatures(xFeatures);
    return xFeatures;
  }

  private List<IntervalFeature> getFeaturesInVCFInterval(String chrom, int from, int to)
      throws IOException, InvalidGenomicCoordsException {

    // Get header if not set yet
    if (this.getVcfHeader() == null) {
      if (Utils.urlFileExists(this.getFilename())) {
        URL url = new URL(this.getFilename());
        AbstractFeatureReader<VariantContext, LineIterator> reader =
            AbstractFeatureReader.getFeatureReader(url.toExternalForm(), new VCFCodec(), false);
        this.setVcfHeader((VCFHeader) reader.getHeader());
      } else {
        VCFFileReader reader = new VCFFileReader(new File(this.getWorkFilename()));
        this.setVcfHeader(reader.getFileHeader());
        reader.close();
      }
    }

    // Collect feature
    List<IntervalFeature> xFeatures = new ArrayList<IntervalFeature>();
    TabixBigBedIterator qry = this.getReader().query(chrom, from - 1, to);
    while (true) {
      String q = qry.next();
      if (q == null) {
        break;
      }
      IntervalFeature intervalFeature =
          new IntervalFeature(q, TrackFormat.VCF, this.getVCFCodec(), this.getScoreColIdx());
      xFeatures.add(intervalFeature);
    }
    return xFeatures;
  }

  protected void removeInvisibleFeatures(List<IntervalFeature> iftList)
      throws InvalidGenomicCoordsException, IOException {

    for (int i = 0; i < iftList.size(); i++) {

      IntervalFeature tr = iftList.get(i);
      String x = tr.getRaw();

      boolean showIt = true;
      if (this.getShowRegex() != null
          && !this.getShowRegex().pattern().equals(Filter.DEFAULT_SHOW_REGEX.getValue())) {
        showIt = this.getShowRegex().matcher(x).find();
      }

      boolean hideIt = false;
      if (this.getHideRegex() != null
          && !this.getHideRegex().pattern().isEmpty()
          && !this.getHideRegex().pattern().equals(Filter.DEFAULT_HIDE_REGEX.getValue())) {
        hideIt = this.getHideRegex().matcher(x).find();
      }
      if (!showIt || hideIt) {
        iftList.set(i, null);
      }
    }
    iftList.removeIf(Objects::isNull);

    if (!this.getAwk().isEmpty() && !this.getAwk().equals(Filter.DEFAULT_AWK)) {
      String[] rawLines = new String[iftList.size()];
      for (int i = 0; i < iftList.size(); i++) {
        rawLines[i] = iftList.get(i).getRaw();
      }
      boolean[] passAwk;

      // Awk
      try {
        passAwk = Utils.passAwkFilter(rawLines, this.getAwk());
      } catch (Exception e) {
        System.err.print(
            Utils.padEndMultiLine(
                "Error processing awk script.", this.getGc().getUserWindowSize()));
        try {
          this.setAwk("");
        } catch (ClassNotFoundException | InvalidRecordException | SQLException e1) {
          e1.printStackTrace();
        }
        throw new IOException();
      }
      for (int i = 0; i < passAwk.length; i++) {
        if (!passAwk[i]) {
          iftList.set(i, null);
        }
      }
      iftList.removeIf(Objects::isNull);
    }
  }

  /**
   * Return the coordinates of the next feature so that the start coincide with the start of the
   * feature and the end is the start + windowSize.
   */
  public GenomicCoords coordsOfNextFeature(GenomicCoords currentGc, boolean getPrevious)
      throws InvalidGenomicCoordsException, IOException {

    IntervalFeature nextFeature;
    if (getPrevious) {
      nextFeature = this.getPreviousFeature(currentGc.getChrom(), currentGc.getFrom());
    } else {
      nextFeature = this.getNextFeature(currentGc.getChrom(), currentGc.getTo());
    }
    if (nextFeature == null) {
      return currentGc;
    }
    GenomicCoords nextGc =
        new GenomicCoords(
            Utils.coordinatesToString(
                nextFeature.getChrom(),
                nextFeature.getFrom(),
                nextFeature.getFrom() + currentGc.getGenomicWindowSize() - 1),
            currentGc.getUserWindowSize(),
            currentGc.getSamSeqDict(),
            currentGc.getFastaFile());
    return nextGc;
  }

  private IntervalFeature getNextFeature(String startChrom, int from)
      throws IOException, InvalidGenomicCoordsException {

    IntervalFeature next = getNextFeatureOnChrom(startChrom, from);
    if (next != null) {
      return next;
    } // There is no feature left on the starting chrom at the given position.

    // Search the remaining chroms in order starting from position 1:
    List<String> chroms = this.getChromListStartingAt(startChrom);
    if (chroms.contains(startChrom)) {
      // The start chroms is searched last. The `if` statement controls for the case where
      // the startChrom is not present at all in this track.
      chroms.remove(startChrom);
      chroms.add(startChrom);
    }
    for (String chrom : chroms) {
      next =
          this.getNextFeatureOnChrom(
              chrom, 0); // Use 0 so if the next feature starts at the beginning of the chrom,
      // i.e. at start=1, it is not missed. See issue #50
      if (next != null) {
        return next;
      }
    }
    return null;
  }

  /**
   * Get the next feature on chrom after "from" position or null if no feature found. This function
   * should be used only by getNextFeature() which is more general since it searches all the
   * chromosomes in turn.
   *
   * @throws IOException
   * @throws InvalidGenomicCoordsException
   */
  private IntervalFeature getNextFeatureOnChrom(String chrom, int from)
      throws IOException, InvalidGenomicCoordsException {

    int qend = (from - 1) < 0 ? 0 : (from - 1);

    TabixBigBedIterator iter = this.getReader().query(chrom, qend, Integer.MAX_VALUE);

    // We accumulate here a number of feature to test for visibility and we return the first one
    // that passes. Testing one by one would be slow.
    List<IntervalFeature> buffer = new ArrayList<IntervalFeature>();

    while (true) {
      String line = iter.next();
      if (buffer.size() < 1000 && line != null) {
        IntervalFeature x =
            new IntervalFeature(
                line, this.getTrackFormat(), this.getVCFCodec(), this.getScoreColIdx());
        buffer.add(x);
        continue;
      }
      if (buffer.size() > 0) {
        this.removeInvisibleFeatures(buffer);
        for (int i = 0; i < buffer.size(); i++) {
          if (buffer.get(i).getFrom() > from) {
            return buffer.get(i);
          }
        }
        buffer.clear();
      }
      if (line == null) {
        return null;
      }
      //            IntervalFeature x= new IntervalFeature(line, this.getTrackFormat(),
      // this.getVCFCodec(), this.getScoreColIdx());
      //            if(x.getFrom() > from && this.featureIsVisible(x.getRaw())){
      //                return x;
      //            }
    }
  }

  private IntervalFeature getPreviousFeature(String startChrom, int pos)
      throws IOException, InvalidGenomicCoordsException {

    IntervalFeature prev = getPreviousFeatureOnChrom(startChrom, pos);
    if (prev != null) {
      return prev;
    }
    // There is no feature left on the starting chrom at the given position.

    // Search the remaining chroms in order starting from last:
    List<String> chroms = Lists.reverse(this.getChromListStartingAt(startChrom));
    for (String chrom : chroms) {
      prev = getPreviousFeatureOnChrom(chrom, Integer.MAX_VALUE);
      if (prev != null) {
        return prev;
      }
    }
    return null;
  }

  /**
   * Experimental method Find the previous feature relative to the given chrom-pos position. Read
   * the chunk of file *before* the pos position and return the last feature found. I.e. find the
   * feature immediately to to left of pos and not overlapping pos.
   */
  private IntervalFeature getPreviousFeatureOnChrom(String chrom, int pos)
      throws IOException, InvalidGenomicCoordsException {

    int chunkFrom = Math.max(0, pos - 1000000);
    int chunkTo = pos - 1; // -1 because we don't include the current position
    IntervalFeature last = null;

    while (chunkTo > 0) {

      TabixBigBedIterator iter = this.getReader().query(chrom, chunkFrom, chunkTo);

      List<IntervalFeature> buffer = new ArrayList<IntervalFeature>();

      while (true) {
        // Find the last feature in this chunk where end coordinate is less then pos
        String line = iter.next();
        if (buffer.size() < 1000 && line != null) {
          IntervalFeature x =
              new IntervalFeature(
                  line, this.getTrackFormat(), this.getVCFCodec(), this.getScoreColIdx());
          buffer.add(x);
          continue;
        }
        if (buffer.size() > 0) {
          this.removeInvisibleFeatures(buffer);
          for (int i = 0; i < buffer.size(); i++) {
            if (buffer.get(i).getTo() < pos) {
              last = buffer.get(i);
            }
          }
          buffer.clear();
        }
        if (line == null) {
          break;
        }
        //                line= line.trim();
        //                IntervalFeature candidate= new IntervalFeature(line,
        // this.getTrackFormat(), this.getVCFCodec(), this.getScoreColIdx());
        //                if(candidate.getTo() < pos && this.featureIsVisible(line)){
        //                    // This is a candidate feature but we don't know yet of it's the last
        // one
        //                    last= new IntervalFeature(line, this.getTrackFormat(),
        // this.getVCFCodec(), this.getScoreColIdx());
        //                }
      }
      if (last != null) {
        break; // The last feature is not null and valid. Stop looking
      } else {
        // Move to previous chunk;
        chunkTo = chunkFrom;
        chunkFrom = Math.max(0, chunkFrom - 1000000);
      }
    }
    if (last == null) {
      return null;
    }
    List<IntervalFeature> xLast = new ArrayList<IntervalFeature>();
    xLast.add(last);
    this.removeInvisibleFeatures(xLast);
    if (last.getTo() > pos || xLast.size() == 0) {
      return null;
    }
    return last;
  }

  protected GenomicCoords startEndOfNextFeature(GenomicCoords currentGc, boolean getPrevious)
      throws InvalidGenomicCoordsException, IOException {
    IntervalFeature nextFeature;
    if (getPrevious) {
      nextFeature = getPreviousFeature(currentGc.getChrom(), currentGc.getFrom());
    } else {
      nextFeature = getNextFeature(currentGc.getChrom(), currentGc.getTo());
    }
    if (nextFeature == null) {
      return currentGc;
    }
    GenomicCoords nextGc =
        new GenomicCoords(
            Utils.coordinatesToString(
                nextFeature.getChrom(), nextFeature.getFrom(), nextFeature.getTo()),
            currentGc.getUserWindowSize(),
            currentGc.getSamSeqDict(),
            currentGc.getFastaFile());
    return nextGc;
  }

  public GenomicCoords findNextMatch(GenomicCoords currentGc, Pattern pattern)
      throws IOException, InvalidGenomicCoordsException {

    IntervalFeature nextFeature =
        findNextRegexInGenome(pattern, currentGc.getChrom(), currentGc.getTo());
    if (nextFeature == null) {
      return currentGc;
    }
    GenomicCoords nextGc =
        new GenomicCoords(
            Utils.coordinatesToString(
                nextFeature.getChrom(),
                nextFeature.getFrom(),
                nextFeature.getFrom() + currentGc.getGenomicWindowSize() - 1),
            currentGc.getUserWindowSize(),
            currentGc.getSamSeqDict(),
            currentGc.getFastaFile());
    return nextGc;
  }

  /**
   * Execute findAllChromRegexInGenome() and return the extreme coordinates of the matched features
   */
  protected GenomicCoords genomicCoordsAllChromMatchInGenome(
      Pattern pattern, GenomicCoords currentGc) throws IOException, InvalidGenomicCoordsException {

    List<IntervalFeature> matchedFeatures = this.findAllChromMatchInGenome(pattern, currentGc);

    if (matchedFeatures.size() == 0) {
      return currentGc;
    }

    // Now get the coords of the first and last feature matched.
    String chrom = matchedFeatures.get(0).getChrom();
    int startFrom = matchedFeatures.get(0).getFrom();
    int endTo = matchedFeatures.get(matchedFeatures.size() - 1).getTo();
    GenomicCoords allMatchesGc =
        new GenomicCoords(
            Utils.coordinatesToString(chrom, startFrom, endTo),
            currentGc.getUserWindowSize(),
            currentGc.getSamSeqDict(),
            currentGc.getFastaFile());
    return allMatchesGc;
  }

  /**
   * Find all the feature matching regex. Only the feature on one chromosome are returned and this
   * chromsome is the first one to have a match. The search starts from the beginning of the current
   * chrom and if nothing is found continues to the other chroms.
   *
   * @throws InvalidGenomicCoordsException
   */
  private List<IntervalFeature> findAllChromMatchInGenome(Pattern pattern, GenomicCoords currentGc)
      throws IOException, InvalidGenomicCoordsException {

    // Accumulate features here
    List<IntervalFeature> matchedFeatures = new ArrayList<IntervalFeature>();

    // We start search from input chrom
    List<String> chromSearchOrder = null;
    chromSearchOrder = getChromListStartingAt(currentGc.getChrom());

    chromSearchOrder.add(currentGc.getChrom());
    for (String curChrom : chromSearchOrder) {

      TabixBigBedIterator iter = this.getReader().query(curChrom, 0, Integer.MAX_VALUE);

      while (true) {
        String line = iter.next();
        if (line == null) {
          break;
        }
        boolean matched = pattern.matcher(line).find();
        if (matched) {
          IntervalFeature x =
              new IntervalFeature(
                  line, this.getTrackFormat(), this.getVCFCodec(), this.getScoreColIdx());
          matchedFeatures.add(x);
          // if(this.featureIsVisible(x.getRaw())){
          //    matchedFeatures.add(x);
          // }
        }
      }
      // Now you have all the features matching regex on this chromosome. We need to exclude those
      // that are not visible
      this.removeInvisibleFeatures(matchedFeatures);

      if (matchedFeatures.size() > 0) {
        // At least one feature matching regex found on this chrom.
        // Check we are at the same position as the beginning. if so, continue to other chroms
        if (matchedFeatures.get(0).getChrom().equals(currentGc.getChrom())
            && matchedFeatures.get(0).getFrom() == currentGc.getFrom()
            && matchedFeatures.get(matchedFeatures.size() - 1).getTo() == currentGc.getTo()) {
          // Discard results and keep searching other chroms.
          matchedFeatures = new ArrayList<IntervalFeature>();
        } else {
          break;
        }
      }
    } // Loop chrom
    return matchedFeatures;
  }

  @Override
  public String printToScreen() throws InvalidGenomicCoordsException, InvalidColourException {

    for (IntervalFeature x : this.getIntervalFeatureList()) {
      if (this.getTrackFormat().equals(TrackFormat.GFF)
          || this.getTrackFormat().equals(TrackFormat.GTF)) {
        x.setGtfAttributeForName(this.gtfAttributeForName);
      } else if (this.getTrackFormat().equals(TrackFormat.BED)
          || this.getTrackFormat().equals(TrackFormat.BIGBED)) {
        x.setBedFieldName(Integer.valueOf(this.bedFieldForName));
      }
    }

    List<String> printable = new ArrayList<String>();
    int nLines = 0;
    try {
      for (List<IntervalFeature> listToPrint : this.stackFeatures()) {
        nLines++;
        if (nLines > this.yMaxLines) {
          // Limit the number of lines in output
          break;
        }
        printable.add(this.printToScreenOneLine(listToPrint));
      }
    } catch (Exception e) {
      e.printStackTrace();
    }

    // Genotype matrix
    if (this.getTrackFormat().equals(TrackFormat.VCF)) {
      try {
        String gtm =
            this.getGenotypeMatrix()
                .printToScreen(
                    this.isNoFormat(),
                    this.intervalFeatureList,
                    this.getGc().getUserWindowSize(),
                    this.getVcfHeader());
        printable.add(gtm);
      } catch (InvalidColourException | IOException e) {
        e.printStackTrace();
      }
    }

    return StringUtils.join(printable, "\n").replaceAll("\n$", "");
  }

  /**
   * Put in the same list reads that will go in the same line of text. This method separates
   * features touching or overlapping each other, useful for visualization. Each item of the output
   * list going on its own line.
   *
   * @param space Space between text feature that touch each other on screen. Use 1 to have at least
   *     one space so that distinct features never look merged with adjacent ones. 0 will not put
   *     any space and will give more compact view. See also TrackReads.stackReads();
   * @throws InvalidGenomicCoordsException
   * @throws IOException
   * @throws InvalidColourException
   */
  private List<List<IntervalFeature>> stackFeatures()
      throws InvalidGenomicCoordsException, IOException, InvalidColourException {

    List<IntervalFeature> intervals;
    List<IntervalFeature> flatListOfTx = this.flatListOfPrintableFeatures();
    if (this.getFeatureDisplayMode().equals(FeatureDisplayMode.COLLAPSED)) {
      intervals = Utils.mergeIntervalFeatures(flatListOfTx, false);
    } else if (this.getFeatureDisplayMode().equals(FeatureDisplayMode.ONELINE)) {
      intervals = Utils.mergeIntervalFeatures(flatListOfTx, true);
    } else {
      intervals = flatListOfTx;
    }

    // Make a copy of the IntervalFeature list. Items will be popped out as they are
    // added to individual lines.
    List<IntervalFeature> flatList = new ArrayList<IntervalFeature>(intervals);

    List<List<IntervalFeature>> listOfLines = new ArrayList<List<IntervalFeature>>();
    if (flatList.size() == 0) {
      return listOfLines;
    }
    List<IntervalFeature> line = new ArrayList<IntervalFeature>();
    line.add(flatList.get(0));
    flatList.remove(0);
    listOfLines.add(line);

    while (true) {
      ArrayList<IntervalFeature> trToRemove = new ArrayList<IntervalFeature>();
      // Find a read in input whose start is greater then end of current
      for (int i = 0; i < flatList.size(); i++) {
        IntervalFeature intervalFeature = flatList.get(i);
        // int gap= 1; // Add a space between book-end features
        if (intervalFeature.getScreenFrom()
            > line.get(line.size() - 1).getScreenTo()
                + this.getGap()) { // +2 because we want some space between adjacent reads
          listOfLines.get(listOfLines.size() - 1).add(intervalFeature); // Append to the last line.
          trToRemove.add(intervalFeature);
        }
      } // At the end of the loop you have put in line as many reads as you can.
      for (IntervalFeature intervalFeature : trToRemove) {
        flatList.remove(flatList.indexOf(intervalFeature));
      }
      // Create a new line, add the first IntervalFeature in list
      if (flatList.size() > 0) {
        line = new ArrayList<IntervalFeature>();
        line.add(flatList.get(0));
        listOfLines.add(line);
        flatList.remove(0);
      } else {
        break;
      }
    }
    return listOfLines;
  }

  /**
   * Return a string of a single line of (typically de-stacked) reads
   *
   * @throws IOException
   * @throws InvalidGenomicCoordsException
   * @throws InvalidColourException
   */
  private String printToScreenOneLine(List<IntervalFeature> listToPrint)
      throws InvalidGenomicCoordsException, IOException, InvalidColourException {

    List<String> printable =
        new ArrayList<
            String>(); // Each item in this list occupies a character space in the terminal.
    // NB: Each item is String not char because it might contain the ansi formatting.
    for (int i = 0; i < this.getGc().getMapping().size(); i++) { // First create empty line
      printable.add(" ");
    }
    for (IntervalFeature intervalFeature : listToPrint) {
      if (intervalFeature.getScreenFrom() == -1) {
        continue; // See test canProcessIndelAtWindowBoundary for how this can happen
      }
      List<FeatureChar> text = intervalFeature.getIdeogram(false, false);

      int i = 0;
      for (int j = intervalFeature.getScreenFrom(); j <= intervalFeature.getScreenTo(); j++) {
        printable.set(j, text.get(i).format(this.isNoFormat()));
        i++;
      }
    }
    return StringUtils.join(printable, "");
  }

  /**
   * List where the original records have been grouped into transcripts. If there are no
   * transcripts, just return the input feature(s) as it is. TODO: Process here also squash, merge
   * and gap?
   *
   * @throws IOException
   * @throws InvalidGenomicCoordsException
   * @throws InvalidColourException
   * @throws InvalidCommandLineException
   */
  private List<IntervalFeature> flatListOfPrintableFeatures()
      throws InvalidGenomicCoordsException, IOException, InvalidColourException {

    for (IntervalFeature x : this.getIntervalFeatureList()) {
      x.getIdeogram(true, false);
    }
    this.changeFeatureColor(this.getColorForRegex());

    List<IntervalFeature> flatList = new ArrayList<IntervalFeature>();

    List<Double> mapToScreen = this.getGc().getMapping();

    if (this.getTrackFormat().equals(TrackFormat.GFF)
        || this.getTrackFormat().equals(TrackFormat.GTF)) {

      Map<String, List<IntervalFeature>> tx;
      if (this.getTrackFormat().equals(TrackFormat.GFF)) {
        tx = this.groupByGFFAttribute();
      } else if (this.getTrackFormat().equals(TrackFormat.GTF)) {
        tx = this.groupByGTFAttribute();
      } else {
        throw new RuntimeException("This is not a GTF or GFF track!");
      }
      for (String txId : tx.keySet()) {
        if (txId.equals("_na_")) {
          // Features that are not part of transcript
          for (IntervalFeature x : tx.get(txId)) {
            x.getIdeogram(false, false);
            flatList.add(x);
          }
        } else {
          flatList.add(this.collapseGFFTranscript(tx.get(txId), mapToScreen));
        }
      }

    } else {
      for (IntervalFeature x : this.getIntervalFeatureList()) {
        flatList.add(x);
      }
    }

    Collections.sort(flatList);
    for (IntervalFeature x : flatList) {
      // Add name to the ideogram of each feature.
      x.getIdeogram(false, true);
    }
    return flatList;
  }

  protected String getUnformattedTitle() {

    String sq = "";
    if (this.getFeatureDisplayMode().equals(FeatureDisplayMode.COLLAPSED)) {
      sq = "; collapsed";
    }
    String gapped = "";
    if (this.getGap() == 0) {
      gapped = "; ungapped";
    }
    String title =
        this.getTrackTag()
            + ";"
            + " N: "
            + this.intervalFeatureList.size()
            + sq
            + gapped
            + this.getTitleForActiveFilters();
    return title;
  }

  @Override
  protected String getTitleForActiveFilters() {
    List<String> title = new ArrayList<String>();
    if (!this.getAwk().equals(Filter.DEFAULT_AWK.getValue())) {
      title.add("awk");
    }
    if (!this.getShowRegex().pattern().equals(Filter.DEFAULT_SHOW_REGEX.getValue())
        || !this.getHideRegex().pattern().equals(Filter.DEFAULT_HIDE_REGEX.getValue())) {
      title.add("grep");
    }
    if (title.size() > 0) {
      return "; filters: " + title.toString();
    } else {
      return "";
    }
  }

  @Override
  public String getTitle()
      throws InvalidColourException, InvalidGenomicCoordsException, IOException {

    if (this.isHideTitle()) {
      return "";
    }
    return this.formatTitle(this.getUnformattedTitle()) + "\n";
  }

  /**
   * Searching the current chrom starting at "from" to find the *next* feature matching the given
   * string. If not found, search the other chroms, if not found restart from the beginning of the
   * current chrom until the "from" position is reached.
   *
   * @throws InvalidGenomicCoordsException
   */
  protected IntervalFeature findNextRegexInGenome(Pattern pattern, String chrom, int from)
      throws IOException, InvalidGenomicCoordsException {

    int startingPoint = from - 1; // -1 because tabix.query from is 0 based (seems so at least)
    List<String> chromSearchOrder = this.getChromListStartingAt(chrom);
    chromSearchOrder.add(chrom);

    for (String curChrom : chromSearchOrder) {

      TabixBigBedIterator iter = this.getReader().query(curChrom, startingPoint, Integer.MAX_VALUE);
      List<IntervalFeature> buffer = new ArrayList<IntervalFeature>();

      while (true) {
        String line = iter.next();
        //                if(line == null) {
        //                    break;
        //                }
        if (line != null) {
          boolean matched = pattern.matcher(line).find();
          if (matched) {
            IntervalFeature x =
                new IntervalFeature(
                    line, this.getTrackFormat(), this.getVCFCodec(), this.getScoreColIdx());
            if (x.getFrom() > startingPoint) {
              buffer.add(x);
            }
            // IntervalFeature x= new IntervalFeature(line, this.getTrackFormat(),
            // this.getVCFCodec(), this.getScoreColIdx());
            // if(x.getFrom() > startingPoint && this.featureIsVisible(x.getRaw())){
            //    return x;
            // }
          }
        }
        if (line == null || buffer.size() > 10) {
          this.removeInvisibleFeatures(buffer);
          if (buffer.size() > 0) {
            return buffer.get(0);
          }
        }
        if (line == null) {
          break;
        }
      }
      startingPoint = 0;
    }
    return null; // Not found anywhere
  }

  private VCFCodec getVCFCodec() {
    if (this.getVcfHeader() == null) {
      return null;
    }
    if (this.vcfCodec == null) {
      VCFCodec vcfCodec = new VCFCodec();
      vcfCodec.setVCFHeader(this.getVcfHeader(), Utils.getVCFHeaderVersion(this.getVcfHeader()));
      this.vcfCodec = vcfCodec;
    }
    return this.vcfCodec;
  }

  /**
   * Return the list of known chroms with startChrom in first position. chroms: chr1 chr2 chr3 chr4
   * chr5 chr6 chr7 startChrom: chr3 return: chr3 chr4 chr5 chr6 chr7 chr1 chr2
   */
  private List<String> getChromListStartingAt(String startChrom) {

    List<String> chroms = this.getChromosomeNames();

    List<String> chromsStartingAt = new ArrayList<String>();

    int idx = chroms.indexOf(startChrom);
    if (idx == -1) {
      // If startChrom is not present at all in the bed/gtf file.
      chromsStartingAt.addAll(chroms);
      chromsStartingAt.add(0, startChrom);
    } else {
      chromsStartingAt.addAll(chroms.subList(idx, chroms.size()));
      chromsStartingAt.addAll(chroms.subList(0, idx));

      // Sanity check
      if (chroms.size() != chromsStartingAt.size()) {
        throw new RuntimeException(
            "Error reordering chroms. Expected "
                + chroms.size()
                + " chroms got "
                + chromsStartingAt.size());
      }
      if (!(chromsStartingAt.containsAll(chroms) && chroms.containsAll(chromsStartingAt))) {
        throw new RuntimeException("Error re-ordering chromosomes");
      }
    }
    if (!(chromsStartingAt.containsAll(chroms))) {
      throw new RuntimeException("Not all known chroms have been included");
    }
    return chromsStartingAt;
  }

  /**
   * Group the features in this genomic window by GFF attribute (typically a transcripts). Features
   * that don't have the attribute make each a length=1 list.
   */
  private Map<String, List<IntervalFeature>> groupByGFFAttribute() {

    // * First collect the IDs of the transcripts

    // Key is transcript ID e.g. ENST00001234.
    // Values is all the IntervalFeatures captured by this ID and part of a transcript.
    // I.e. their are in txFeature set,
    Map<String, List<IntervalFeature>> txIds = new LinkedHashMap<String, List<IntervalFeature>>();

    // This key:value is for records which are not part transcripts. E.g. features like "chromosome"
    // or rRNA.
    txIds.put("_na_", new ArrayList<IntervalFeature>());

    // Now populate the lists of values by assigning to each key the transcript records:
    for (IntervalFeature x : this.getIntervalFeatureList()) {

      if (FormatGTF.getTxSuperFeatures().contains(x.getFeature().toLowerCase())) {

        // Transcript feature. E.g.
        // 7 ensembl_havana mRNA 5527151 5530709 . - .
        // ID=transcript:ENST00000331789;Parent=gene:ENSG00000075624;Name=ACTB-...
        String txId = x.getGFFValueFromKey("ID");
        if (!txIds.containsKey(txId)) {
          txIds.put(txId, new ArrayList<IntervalFeature>());
        }
        txIds.get(txId).add(x);
      } else if (FormatGTF.getTxSubFeatures().contains(x.getFeature().toLowerCase())) {
        // Part of transcript, e.g:
        // 7 ensembl_havana exon 5527151 5527891 . - .
        // Parent=transcript:ENST00000331789;Name=ENSE00001902654;constitutive=0;ensembl_end_pha
        String txId = x.getGFFValueFromKey("Parent");
        if (txId == null) {
          txId = "_na_";
        }
        if (!txIds.containsKey(txId)) {
          txIds.put(txId, new ArrayList<IntervalFeature>());
        }
        txIds.get(txId).add(x);
      } else {
        // Not a transcript or part thereof. E.g.
        // 7 . biological_region 5529708 5529709 0.999 - . logic_name=eponine
        txIds.get("_na_").add(x);
      }
    }
    // We don't need to return the full Map, only the list of lists (groups) would suffice.
    // However, we need to separate the group of non-trascripts (_na_ key)
    return txIds;
  }

  /** Group the features in this genomic window by GTF attribute (typically a transcripts). */
  private Map<String, List<IntervalFeature>> groupByGTFAttribute() {

    // * First collect the IDs of the transcripts

    // Key is transcript ID e.g. ENST00001234.
    // Values is all the IntervalFeatures captured by this ID and part of a transcript.
    // I.e. their are in txFeature set,
    Map<String, List<IntervalFeature>> txIds = new LinkedHashMap<String, List<IntervalFeature>>();

    // This key:value is for records which are not part transcripts. E.g. features like "chromosome"
    // or rRNA.
    txIds.put("_na_", new ArrayList<IntervalFeature>());

    // Now populate the lists of values by assigning to each key the transcript records:
    for (IntervalFeature x : this.getIntervalFeatureList()) {

      if (FormatGTF.getTxSuperFeatures().contains(x.getFeature().toLowerCase())
          || FormatGTF.getTxSubFeatures().contains(x.getFeature().toLowerCase())) {
        // Transcript feature. E.g.
        // chr7 hg19_wgEncodeGencodeBasicV19 exon       5566782 5567522 0.000000 - . gene_id
        // "ENST00000331789.5"; transcript_id "ENST00000331789.5";
        String txId = x.getGFFValueFromKey("transcript_id");
        if (!txIds.containsKey(txId)) {
          txIds.put(txId, new ArrayList<IntervalFeature>());
        }
        txIds.get(txId).add(x);
      } else {
        // Not a transcript or part thereof. E.g.
        // 7 . biological_region 5529708 5529709 0.999 - . logic_name=eponine
        txIds.get("_na_").add(x);
      }
    }
    // We don't need to return the full Map, only the list of lists (groups) would suffice.
    // However, we need to separate the group of non-trascripts (_na_ key)
    return txIds;
  }

  /**
   * Collapse the list of features in a single IntervalFeature representing the transcript. The
   * elements of txFeatures are expected to represent the entire transcript, nothing more (e.g.
   * "chromosome"). The transcript may not be biologically complete as part of it may be outside the
   * current genomic coords.
   *
   * <p>mapToScreen: Mapping of genomic coordinates to screen coordinates. This could be obtained
   * inside this method but better to pass it from outside as it can take time to get the terminal
   * window size several times.
   *
   * @throws InvalidGenomicCoordsException
   * @throws InvalidColourException
   */
  private IntervalFeature collapseGFFTranscript(
      List<IntervalFeature> txFeatures, List<Double> mapToScreen)
      throws InvalidGenomicCoordsException, InvalidColourException {

    if (txFeatures.size() == 0) {
      System.err.println("Unexpected transcript: Length zero!");
      throw new RuntimeException();
    }

    // Collect the genomic and screen coordinates of this transcript
    int gFrom = Integer.MAX_VALUE;
    int gTo = 0;
    int screenFrom = Integer.MAX_VALUE;
    int screenTo = 0;
    for (IntervalFeature x : txFeatures) {

      if (x.getFrom() < gFrom) {
        gFrom = x.getFrom();
      }
      if (x.getTo() > gTo) {
        gTo = x.getTo();
      }
      if (x.getScreenFrom() < screenFrom) {
        screenFrom = x.getScreenFrom();
      }
      if (x.getScreenTo() > screenTo) {
        screenTo = x.getScreenTo();
      }
    }

    IntervalFeature transcript =
        new IntervalFeature(txFeatures.get(0).getChrom(), gFrom, gTo, TrackFormat.GFF);
    transcript.setStrand(txFeatures.get(0).getStrand());
    transcript.mapToScreen(mapToScreen);

    // Now we need to prepare the ideogram
    int txIdeogramSize = screenTo - screenFrom + 1;
    List<FeatureChar> ideogram = new ArrayList<FeatureChar>(txIdeogramSize);
    for (int i = 0; i < txIdeogramSize; i++) {
      FeatureChar c = new FeatureChar();
      c.setText('-');
      ideogram.add(c); // Default character to print. Typically this should apply to introns only.
    }

    for (String txSubType : FormatGTF.getTxSubFeatures()) {

      for (IntervalFeature subFeature : txFeatures) {

        if (subFeature.getFeature().toLowerCase().equals(txSubType)) {
          // Replace the featureChars in the novel transcript with the those from the individual
          // features
          // cccccccc                  <- subfeature#1
          //               eeee     <- subfeature#2
          // ---------------------- <- novel ideogram to be replaced
          List<FeatureChar> subFeatureIdeogram = subFeature.getIdeogram(false, false);
          int offset = subFeature.getScreenFrom() - screenFrom;
          for (FeatureChar x : subFeatureIdeogram) {
            ideogram.set(offset, x);
            offset++;
          }
        }
      }
    }
    // Now we get the name for this transcript
    String txName = "."; // Default: No name
    outerloop:
    for (String txSuperType : FormatGTF.getTxSuperFeatures()) {
      for (IntervalFeature x : txFeatures) {
        if (x.getFeature().toLowerCase().equals(txSuperType)) {
          txName = x.getName();
        }
        if (txName != null && !txName.isEmpty() && !txName.equals(".")) {
          break outerloop; // A name found, break
        }
      }
    }
    if (txName == null || txName.isEmpty() || txName.equals(".")) {
      // If a name has not been found among the superfeatures, look at the
      // individual components (exons, CDS, etc)
      outerloop:
      for (String txSuperType : FormatGTF.getTxSubFeatures()) {
        for (IntervalFeature x : txFeatures) {
          if (x.getFeature().toLowerCase().equals(txSuperType)) {
            txName = x.getName();
          }
          if (txName != null && !txName.isEmpty() && !txName.equals(".")) {
            break outerloop; // A name found, break
          }
        }
      }
    }
    transcript.setName(txName);
    transcript.setIdeogram(ideogram, false);
    return transcript;
  }

  protected List<IntervalFeature> getIntervalFeatureList() {
    return intervalFeatureList;
  }

  protected void setIntervalFeatureList(List<IntervalFeature> intervalFeatureList) {
    this.intervalFeatureList = intervalFeatureList;
  }

  @Override
  public ArrayList<String> getChromosomeNames() {
    ArrayList<String> x = new ArrayList<String>(this.getReader().getChromosomes());
    Collections.sort(x);
    return x;
  }

  private TabixBigBedReader getReader() {

    if (this.bigBedReader != null && this.tabixReader != null) {
      System.err.println("You cannot have both tabix and bigBed readers set!");
      throw new RuntimeException();
    }

    if (this.bigBedReader != null) {
      return new TabixBigBedReader(this.bigBedReader);
    } else if (this.tabixReader != null) {
      return new TabixBigBedReader(this.tabixReader);
    } else {
      System.err.println("Tabix and bigBed reader both null.");
      throw new RuntimeException();
    }
  }

  private TabixReader getTabixReader(String tabixFile) throws IOException {
    return new TabixReader(new File(tabixFile).getAbsolutePath());
  }

  /** This setter is for TrackBookmark to work. */
  protected void setTabixReader(TabixReader tabixReader) {
    this.tabixReader = tabixReader;
  }

  protected TabixReader getTabixReader() {
    return this.tabixReader;
  }

  @Override
  protected List<String> getRecordsAsStrings() {

    List<String> featureList = new ArrayList<String>();
    for (IntervalFeature ift : intervalFeatureList) {
      if (this.getTrackFormat().equals(TrackFormat.VCF) && this.getPrintNormalizedVcf()) {
        List<String> line =
            this.normalizeVcfRecordBySample(
                this.getVcfHeader().getSampleNamesInOrder(), ift.getRaw());
        featureList.addAll(line);
      } else {
        featureList.add(ift.getRaw());
      }
    }
    return featureList;
  }

  private List<String> normalizeVcfRecordBySample(List<String> sampleNames, String rawVcfLine) {
    List<String> tsv = new ArrayList<String>();
    if (sampleNames.size() == 0) {
      tsv.add(rawVcfLine);
      return tsv;
    }
    List<String> vcfList = Splitter.on("\t").splitToList(rawVcfLine);
    for (int i = 0; i < sampleNames.size(); i++) {
      List<String> samples = vcfList.subList(9, vcfList.size());
      List<String> fixed = vcfList.subList(0, 8);
      String fmtTags = vcfList.get(8);
      String tabLine =
          Joiner.on("\t").join(fixed)
              + "\t"
              + sampleNames.get(i)
              + "\t"
              + fmtTags
              + "\t"
              + samples.get(i);
      tsv.add(tabLine);
    }
    return tsv;
  }

  @Override
  protected void setColorForRegex(List<Argument> xcolorForRegex) {
    if (xcolorForRegex == null) {
      this.colorForRegex = null;
      return;
    } else {
      if (this.colorForRegex == null) {
        this.colorForRegex = new ArrayList<Argument>();
      }
      for (Argument p : xcolorForRegex) {
        this.colorForRegex.add(p);
      }
    }
  }

  private List<Argument> getColorForRegex() {
    return this.colorForRegex;
  }

  @Override
  protected void changeFeatureColor(List<Argument> list)
      throws InvalidColourException, IOException {
    if (list == null) {
      return;
    }

    for (Argument arg : list) {
      String regex = arg.getKey();
      String color = arg.getArg();
      String[] rawrecs = new String[this.getIntervalFeatureList().size()];
      for (int i = 0; i < this.getIntervalFeatureList().size(); i++) {
        rawrecs[i] = this.getIntervalFeatureList().get(i).getRaw();
      }
      boolean[] matched = Utils.matchByAwkOrRegex(rawrecs, regex);
      for (int i = 0; i < matched.length; i++) {
        boolean m = matched[i];
        if (arg.isInvert()) {
          m = !m;
        }
        if (m) {
          for (FeatureChar f : this.getIntervalFeatureList().get(i).getIdeogram(false, false)) {
            f.setBgColor(color);
            f.setFgColor(Xterm256.getContrastColor(color));
          }
        }
      }
    }
  }

  @Override
  public void reload()
      throws InvalidGenomicCoordsException,
          IOException,
          ClassNotFoundException,
          InvalidRecordException,
          SQLException {
    if (!Files.isSameFile(Paths.get(this.getWorkFilename()), Paths.get(this.getFilename()))) {
      TrackIntervalFeature tr = new TrackIntervalFeature(this.getFilename(), this.getGc());
      String fname = this.getWorkFilename();
      Files.move(
          Paths.get(tr.getWorkFilename()),
          Paths.get(fname),
          java.nio.file.StandardCopyOption.REPLACE_EXISTING);
      Files.move(
          Paths.get(tr.getWorkFilename() + FileExtensions.TABIX_INDEX),
          Paths.get(fname + FileExtensions.TABIX_INDEX),
          java.nio.file.StandardCopyOption.REPLACE_EXISTING);
    }
    this.tabixReader = this.getTabixReader(this.getWorkFilename());
    this.update();
  }

  @Override
  public void setFeatureName(String nameFieldOrAttribute) {
    if (this.getTrackFormat().equals(TrackFormat.GFF)
        || this.getTrackFormat().equals(TrackFormat.GTF)) {
      this.gtfAttributeForName = nameFieldOrAttribute;
    } else if (this.getTrackFormat().equals(TrackFormat.BED)
        || this.getTrackFormat().equals(TrackFormat.BIGBED)) {
      if (nameFieldOrAttribute.equals("-na")) {
        this.bedFieldForName = -1;
      } else {
        try {
          this.bedFieldForName =
              Integer.valueOf(nameFieldOrAttribute)
                  - 1; // User's input is 1-based, convert ot 0-based
        } catch (NumberFormatException e) {
          System.err.println("Cannot convert " + nameFieldOrAttribute + " to integer");
          throw e;
        }
      }
    }
  }

  protected int getScoreColIdx() {
    return scoreColIdx;
  }

  protected void setScoreColIdx(int scoreColIdx)
      throws ClassNotFoundException,
          IOException,
          InvalidGenomicCoordsException,
          InvalidRecordException,
          SQLException {
    this.scoreColIdx = scoreColIdx;
    this.update();
  }
}
