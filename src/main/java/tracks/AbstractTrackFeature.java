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
import java.util.*;
import java.util.regex.Pattern;
import org.apache.commons.lang3.StringUtils;
import org.broad.igv.bbfile.BBFileReader;
import samTextViewer.GenomicCoords;
import samTextViewer.Utils;

public abstract class AbstractTrackFeature<T extends IntervalFeature>
    extends AbstractTrack {
  /** For GTF/GFF data: Use this attribute to get the feature names */
  protected TabixReader tabixReader; // Leave *protected* for TrackBookmark to work

  protected BBFileReader bigBedReader;
  protected List<T> intervalFeatureList = new ArrayList<>();
  private List<Argument> colourForRegex = null;

  protected abstract T createFeature(String line) throws InvalidGenomicCoordsException;
  protected abstract Map<String, List<T>> groupByGFFAttribute();
  protected abstract Map<String, List<T>> groupByGTFAttribute();
  protected abstract T collapseGFFTranscript(List<T> features, List<Double> mapToScreen)
      throws InvalidGenomicCoordsException, InvalidColourException;

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
    for (T ift : this.intervalFeatureList) {
      ift.mapToScreen(this.getGc().getMapping());
    }
  }

  protected List<T> getFeaturesInInterval(String chrom, int from, int to)
      throws IOException, InvalidGenomicCoordsException {
    if (from < 1) {
      System.err.println("from < 1: " + from + "; resetting to 1.");
      from = 1;
    }
    if (from > to) {
      throw new InvalidGenomicCoordsException("Invalid coordinates: from: "
          + from
          + "; to: "
          + to
          + "; Resetting to initial 1-"
          + Integer.MAX_VALUE);
    }
    List<T> xFeatures = new ArrayList<>();
    this.removeInvisibleFeatures(xFeatures);
    return xFeatures;
  }

  protected void removeInvisibleFeatures(List<T> iftList)
      throws InvalidGenomicCoordsException, IOException {

    for (int i = 0; i < iftList.size(); i++) {

      T tr = iftList.get(i);
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

    if (!this.getAwk().isEmpty()) {
      this.getAwk();
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

    T nextFeature;
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

  protected T getNextFeature(String startChrom, int from)
      throws IOException, InvalidGenomicCoordsException {

    T next = this.getNextFeatureOnChrom(startChrom, from);
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
  protected T getNextFeatureOnChrom(String chrom, int from)
      throws IOException, InvalidGenomicCoordsException {

    int qend = Math.max((from - 1), 0);

    TabixBigBedIterator iter = this.getReader().query(chrom, qend, Integer.MAX_VALUE);

    // We accumulate here a number of feature to test for visibility and we return the first one
    // that passes. Testing one by one would be slow.
    List<T> buffer = new ArrayList<>();

    while (true) {
      String line = iter.next();
      if (buffer.size() < 1000 && line != null) {
        T x = this.createFeature(line);
        buffer.add(x);
        continue;
      }
      if (!buffer.isEmpty()) {
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
    }
  }

  private T getPreviousFeature(String startChrom, int pos)
      throws IOException, InvalidGenomicCoordsException {

    T prev = getPreviousFeatureOnChrom(startChrom, pos);
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
  private T getPreviousFeatureOnChrom(String chrom, int pos)
      throws IOException, InvalidGenomicCoordsException {

    int chunkFrom = Math.max(0, pos - 1000000);
    int chunkTo = pos - 1; // -1 because we don't include the current position
    T last = null;

    while (chunkTo > 0) {

      TabixBigBedIterator iter = this.getReader().query(chrom, chunkFrom, chunkTo);

      List<T> buffer = new ArrayList<T>();

      while (true) {
        // Find the last feature in this chunk where end coordinate is less then pos
        String line = iter.next();
        if (buffer.size() < 1000 && line != null) {
          T x = createFeature(line);
          buffer.add(x);
          continue;
        }
        if (!buffer.isEmpty()) {
          this.removeInvisibleFeatures(buffer);
          for (T intervalFeature : buffer) {
            if (intervalFeature.getTo() < pos) {
              last = intervalFeature;
            }
          }
          buffer.clear();
        }
        if (line == null) {
          break;
        }
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
    List<T> xLast = new ArrayList<>();
    xLast.add(last);
    this.removeInvisibleFeatures(xLast);
    if (last.getTo() > pos || xLast.isEmpty()) {
      return null;
    }
    return last;
  }

  protected GenomicCoords startEndOfNextFeature(GenomicCoords currentGc, boolean getPrevious)
      throws InvalidGenomicCoordsException, IOException {
    T nextFeature;
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

    T nextFeature = findNextRegexInGenome(pattern, currentGc.getChrom(), currentGc.getTo());
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

    List<T> matchedFeatures = this.findAllChromMatchInGenome(pattern, currentGc);

    if (matchedFeatures.isEmpty()) {
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
  private List<T> findAllChromMatchInGenome(Pattern pattern, GenomicCoords currentGc)
      throws IOException, InvalidGenomicCoordsException {

    // Accumulate features here
    List<T> matchedFeatures = new ArrayList<T>();

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
          T x = this.createFeature(line);
          matchedFeatures.add(x);
        }
      }
      // Now you have all the features matching regex on this chromosome. We need to exclude those
      // that are not visible
      this.removeInvisibleFeatures(matchedFeatures);

      if (!matchedFeatures.isEmpty()) {
        // At least one feature matching regex found on this chrom.
        // Check we are at the same position as the beginning. if so, continue to other chroms
        if (matchedFeatures.get(0).getChrom().equals(currentGc.getChrom())
            && matchedFeatures.get(0).getFrom() == currentGc.getFrom()
            && matchedFeatures.get(matchedFeatures.size() - 1).getTo() == currentGc.getTo()) {
          // Discard results and keep searching other chroms.
          matchedFeatures = new ArrayList<>();
        } else {
          break;
        }
      }
    } // Loop chrom
    return matchedFeatures;
  }

  protected List<List<T>> stackFeatures()
      throws InvalidGenomicCoordsException, IOException, InvalidColourException {

    List<T> intervals;
    List<T> flatListOfTx = this.flatListOfPrintableFeatures();
    if (this.getFeatureDisplayMode().equals(FeatureDisplayMode.COLLAPSED)) {
      intervals = Utils.mergeIntervalFeatures(flatListOfTx, false);
    } else if (this.getFeatureDisplayMode().equals(FeatureDisplayMode.ONELINE)) {
      intervals = Utils.mergeIntervalFeatures(flatListOfTx, true);
    } else {
      intervals = flatListOfTx;
    }

    // Make a copy of the T list. Items will be popped out as they are
    // added to individual lines.
    List<T> flatList = new ArrayList<>(intervals);

    List<List<T>> listOfLines = new ArrayList<>();
    if (flatList.isEmpty()) {
      return listOfLines;
    }
    List<T> line = new ArrayList<>();
    line.add(flatList.get(0));
    flatList.remove(0);
    listOfLines.add(line);

    while (true) {
      ArrayList<T> trToRemove = new ArrayList<T>();
      // Find a read in input whose start is greater then end of current
      for (int i = 0; i < flatList.size(); i++) {
        T intervalFeature = flatList.get(i);
        if (intervalFeature.getScreenFrom()
            > line.get(line.size() - 1).getScreenTo()
                + this.getGap()) { // +2 because we want some space between adjacent reads
          listOfLines.get(listOfLines.size() - 1).add(intervalFeature); // Append to the last line.
          trToRemove.add(intervalFeature);
        }
      } // At the end of the loop you have put in line as many reads as you can.
      for (T intervalFeature : trToRemove) {
        flatList.remove(intervalFeature);
      }
      // Create a new line, add the first T in list
      if (!flatList.isEmpty()) {
        line = new ArrayList<T>();
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
  protected String printToScreenOneLine(List<T> listToPrint) throws InvalidColourException {

    List<String> printable =
        new ArrayList<
            String>(); // Each item in this list occupies a character space in the terminal.
    // NB: Each item is String not char because it might contain the ansi formatting.
    for (int i = 0; i < this.getGc().getMapping().size(); i++) { // First create empty line
      printable.add(" ");
    }
    for (T intervalFeature : listToPrint) {
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
  private List<T> flatListOfPrintableFeatures()
      throws InvalidGenomicCoordsException, IOException, InvalidColourException {

    for (T x : this.getIntervalFeatureList()) {
      x.getIdeogram(true, false);
    }
    this.changeFeatureColour(this.getColourForRegex());

    List<T> flatList = new ArrayList<T>();

    List<Double> mapToScreen = this.getGc().getMapping();

    if (this.getTrackFormat().equals(TrackFormat.GFF)
        || this.getTrackFormat().equals(TrackFormat.GTF)) {

      Map<String, List<T>> tx;
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
          for (T x : tx.get(txId)) {
            x.getIdeogram(false, false);
            flatList.add(x);
          }
        } else {
          flatList.add(this.collapseGFFTranscript(tx.get(txId), mapToScreen));
        }
      }

    } else {
      for (T x : this.getIntervalFeatureList()) {
        flatList.add(x);
      }
    }

    Collections.sort(flatList);
    for (T x : flatList) {
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
    if (!title.isEmpty()) {
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
  protected T findNextRegexInGenome(Pattern pattern, String chrom, int from)
      throws IOException, InvalidGenomicCoordsException {

    int startingPoint = from - 1; // -1 because tabix.query from is 0 based (seems so at least)
    List<String> chromSearchOrder = this.getChromListStartingAt(chrom);
    chromSearchOrder.add(chrom);

    for (String curChrom : chromSearchOrder) {

      TabixBigBedIterator iter = this.getReader().query(curChrom, startingPoint, Integer.MAX_VALUE);
      List<T> buffer = new ArrayList<>();

      while (true) {
        String line = iter.next();
        if (line != null) {
          boolean matched = pattern.matcher(line).find();
          if (matched) {
            T x = this.createFeature(line);
            if (x.getFrom() > startingPoint) {
              buffer.add(x);
            }
          }
        }
        if (line == null || buffer.size() > 10) {
          this.removeInvisibleFeatures(buffer);
          if (!buffer.isEmpty()) {
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

  abstract protected List<T> getIntervalFeatureList();

  protected void setIntervalFeatureList(List<T> intervalFeatureList) {
    this.intervalFeatureList = intervalFeatureList;
  }

  @Override
  public ArrayList<String> getChromosomeNames() {
    ArrayList<String> x = new ArrayList<String>(this.getReader().getChromosomes());
    Collections.sort(x);
    return x;
  }

  TabixBigBedReader getReader() {

    if (this.bigBedReader != null && this.tabixReader != null) {
      throw new RuntimeException("You cannot have both tabix and bigBed readers set!");
    }

    if (this.bigBedReader != null) {
      return new TabixBigBedReader(this.bigBedReader);
    } else if (this.tabixReader != null) {
      return new TabixBigBedReader(this.tabixReader);
    } else {
      throw new RuntimeException("Tabix and bigBed reader both null.");
    }
  }

  protected TabixReader getTabixReader(String tabixFile) throws IOException {
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

    List<String> featureList = new ArrayList<>();
    for (T ift : intervalFeatureList) {
      featureList.add(ift.getRaw());
    }
    return featureList;
  }

  @Override
  protected void setColourForRegex(List<Argument> xcolourForRegex) {
    if (xcolourForRegex == null) {
      this.colourForRegex = null;
      return;
    } else {
      if (this.colourForRegex == null) {
        this.colourForRegex = new ArrayList<Argument>();
      }
      for (Argument p : xcolourForRegex) {
        this.colourForRegex.add(p);
      }
    }
  }

  private List<Argument> getColourForRegex() {
    return this.colourForRegex;
  }

  @Override
  protected void changeFeatureColour(List<Argument> list)
      throws InvalidColourException, IOException {
    if (list == null) {
      return;
    }

    for (Argument arg : list) {
      String regex = arg.getKey();
      String colour = arg.getArg();
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
            f.setBgColour(colour);
            f.setFgColour(Xterm256.getContrastColour(colour));
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

  private VCFCodec getVCFCodec() {
    return null;
  }

  @Override
  public void close() {
    if (this.tabixReader != null) {
      this.tabixReader.close();
    }
    if (this.bigBedReader != null) {
      this.bigBedReader.close();
    }
  }
}
