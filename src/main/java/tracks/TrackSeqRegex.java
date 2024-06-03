package tracks;

import exceptions.InvalidColourException;
import exceptions.InvalidGenomicCoordsException;
import exceptions.InvalidRecordException;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.tribble.index.tabix.TabixFormat;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import org.apache.commons.io.FileUtils;
import samTextViewer.GenomicCoords;
import samTextViewer.Utils;
import sortBgzipIndex.MakeTabixIndex;

public class TrackSeqRegex extends TrackIntervalFeature {

  private String seqRegex = "a^"; // Match nothing
  private final String noRe = "a^";
  private boolean isCaseSensitive = false;
  private boolean isIupac = false;

  public TrackSeqRegex() {}

  public TrackSeqRegex(GenomicCoords gc)
      throws ClassNotFoundException,
          IOException,
          InvalidGenomicCoordsException,
          InvalidRecordException,
          SQLException {
    super(gc);

    this.setGc(gc);
    this.setFilename(new File(gc.getFastaFile()).getAbsolutePath());
    this.setWorkFilename(new File(gc.getFastaFile()).getAbsolutePath());
    this.setTrackTag(new File(gc.getFastaFile()).getName());
    this.setHideTrack(true);
    this.setTrackFormat(TrackFormat.BED);
  }

  @Override
  /** Find regex matches and update the screen map. See also parent method. */
  public void update()
      throws ClassNotFoundException,
          IOException,
          InvalidGenomicCoordsException,
          InvalidRecordException,
          SQLException {

    this.findRegex();
    for (IntervalFeature ift : this.getIntervalFeatureList()) {
      ift.mapToScreen(this.getGc().getMapping());
    }
  }

  private void writeRegionsToFile(BufferedWriter wr, Set<String> regionListPos, String strand)
      throws InvalidGenomicCoordsException, IOException {

    List<IntervalFeature> iftList = new ArrayList<IntervalFeature>();

    for (String reg : regionListPos) {
      reg = reg + "\t.\t" + strand;
      IntervalFeature x =
          new IntervalFeature(reg, this.getTrackFormat(), null, this.getScoreColIdx());
      iftList.add(x);
    }
    this.removeInvisibleFeatures(iftList);
    for (int i = 0; i < iftList.size(); i++) {
      wr.write(iftList.get(i).getRaw() + "\n");
    }
  }

  /** Find regex matches in current genomic interval and update the IntervalFeature set and list. */
  private void findRegex()
      throws IOException,
          InvalidGenomicCoordsException,
          ClassNotFoundException,
          InvalidRecordException,
          SQLException {

    // Find matches
    // ============
    Pattern pattern = Pattern.compile(this.seqRegex, Pattern.CASE_INSENSITIVE);
    if (this.isCaseSensitive) {
      pattern = Pattern.compile(this.seqRegex);
    }

    byte[] seq = this.getGc().getSequenceFromFasta();
    Matcher matcher = pattern.matcher(new String(seq));

    // One list for matches on forward, one for reverse, one for palindromic
    Set<String> regionListPos = new HashSet<String>();
    Set<String> regionListNeg = new HashSet<String>();
    Set<String> regionListPalind = new HashSet<String>();

    // Forward match
    while (matcher.find()) {
      int matchStart = this.getGc().getFrom() + matcher.start() - 1;
      int matchEnd = this.getGc().getFrom() + matcher.end() - 1;
      String reg =
          this.getGc().getChrom()
              + "\t"
              + matchStart
              + "\t"
              + matchEnd
              + "\t"
              + this.trimMatch(matcher.group(), 100);
      regionListPos.add(reg);
    }
    // Reverse comp match
    SequenceUtil.reverseComplement(seq);

    matcher = pattern.matcher(new String(seq));

    while (matcher.find()) {
      int matchStart = this.getGc().getTo() - matcher.end();
      int matchEnd = this.getGc().getTo() - matcher.start();
      String reg =
          this.getGc().getChrom()
              + "\t"
              + matchStart
              + "\t"
              + matchEnd
              + "\t"
              + this.trimMatch(matcher.group(), 100);
      if (regionListPos.contains(reg)) {
        regionListPos.remove(reg);
        regionListPalind.add(reg);
      } else {
        regionListNeg.add(reg);
      }
    }

    // Prepare tmp file
    // ================
    String tmpname =
        Utils.createTempFile(".asciigenome.regexMatchTrack", ".tmp.bed", true).getAbsolutePath();
    // This is a hack: Copy the newly created tmp file to another file. This overcomes some
    // permission (?) problems later with the tabix indexing.
    String regexMatchFile = tmpname.replaceAll("\\.tmp\\.bed$", ".bed");
    FileUtils.copyFile(new File(tmpname), new File(regexMatchFile));
    new File(tmpname).delete();
    new File(regexMatchFile).deleteOnExit();
    BufferedWriter wr = new BufferedWriter(new FileWriter(new File(regexMatchFile)));

    // Write sets of matches to file
    // =============================

    this.writeRegionsToFile(wr, regionListPos, "+");
    this.writeRegionsToFile(wr, regionListNeg, "-");
    this.writeRegionsToFile(wr, regionListPalind, ".");
    wr.close();

    //		for(String reg : regionListPos){
    //			reg += "\t.\t+\n";
    //			if(this.featureIsVisible(reg)){
    //				wr.write(reg);
    //			}
    //		}
    //		for(String reg : regionListNeg){
    //			reg += "\t.\t-\n";
    //			if(this.featureIsVisible(reg)){
    //				wr.write(reg);
    //			}
    //		}
    //		for(String reg : regionListPalind){
    //			reg += "\t.\t.\n";
    //			if(this.featureIsVisible(reg)){
    //				wr.write(reg);
    //			}
    //		}

    // Compress, index, read back as list of IntervalFeatures
    // ======================================================
    File regexMatchBgzip = new File(regexMatchFile + ".gz");
    File regexMatchIndex = new File(regexMatchFile + ".gz.tbi");
    regexMatchBgzip.deleteOnExit();
    regexMatchIndex.deleteOnExit();

    new MakeTabixIndex(regexMatchFile, regexMatchBgzip, TabixFormat.BED);
    new File(regexMatchFile).delete();

    TrackIntervalFeature regexMatchTrack =
        new TrackIntervalFeature(regexMatchBgzip.getAbsolutePath(), this.getGc());
    regexMatchBgzip.delete();
    regexMatchIndex.delete();

    this.intervalFeatureList = regexMatchTrack.getIntervalFeatureList();
  }

  /** Trim the String x if longer then x and return it with the trimmed part annotated */
  private String trimMatch(String x, int maxLen) {
    if (x.length() > maxLen) {
      x = x.substring(0, maxLen) + "[" + (x.length() - maxLen) + "]";
    }
    return x;
  }

  /**
   * Convert iupac chars in regex string to corresponding regular expression
   *
   * @throws SQLException
   * @throws InvalidRecordException
   * @throws InvalidGenomicCoordsException
   * @throws IOException
   * @throws ClassNotFoundException
   */
  private String iupac2Regex(String rawRegex)
      throws ClassNotFoundException,
          IOException,
          InvalidGenomicCoordsException,
          InvalidRecordException,
          SQLException {

    // Prepare map
    // ---------------------------
    HashMap<Character, String> iupacToRe = new HashMap<Character, String>();
    iupacToRe.put('R', "[AG]");
    iupacToRe.put('Y', "[CT]");
    iupacToRe.put('S', "[GC]");
    iupacToRe.put('W', "[AT]");
    iupacToRe.put('K', "[GT]");
    iupacToRe.put('M', "[AC]");
    iupacToRe.put('B', "[CGT]");
    iupacToRe.put('D', "[AGT]");
    iupacToRe.put('H', "[ACT]");
    iupacToRe.put('V', "[ACG]");
    iupacToRe.put('N', "[ACGT]");

    // To preserve case: Same keys in lowercase.
    Set<Character> upper = new HashSet<Character>();
    for (char X : iupacToRe.keySet()) {
      upper.add(X);
    }
    for (char X : upper) {
      iupacToRe.put(Character.toLowerCase(X), iupacToRe.get(X).toLowerCase());
    }
    // ---------------------------

    String iupacRe = "";
    for (char c : rawRegex.toCharArray()) {
      if (iupacToRe.containsKey(c)) {
        iupacRe += iupacToRe.get(c);
      } else {
        iupacRe += c;
      }
    }
    return iupacRe;
  }

  /**
   * @param append Append to existing file.
   * @throws IOException
   */
  protected void saveIntervalsToFile(String filename, boolean append) throws IOException {

    BufferedWriter wr = new BufferedWriter(new FileWriter(filename, append));

    for (IntervalFeature feature : this.getIntervalFeatureList()) {
      wr.write(feature.getRaw() + "\n");
    }
    wr.close();
  }

  @Override
  public String getSeqRegex() {
    return seqRegex;
  }

  @Override
  public void setSeqRegex(String regex)
      throws ClassNotFoundException,
          IOException,
          InvalidGenomicCoordsException,
          InvalidRecordException,
          SQLException {
    if (regex.isEmpty()) {
      this.setHideTrack(true);
      regex = noRe; // Match nothing
    } else {
      this.setHideTrack(false);
    }

    // Convert as required
    if (this.isIupac) {
      regex = this.iupac2Regex(regex);
    }

    this.seqRegex = regex;
    this.update();
  }

  @Override
  public String getTitle()
      throws InvalidColourException, InvalidGenomicCoordsException, IOException {
    String title = this.getUnformattedTitle() + "; re: " + this.getSeqRegex();
    // title= Utils.padEndMultiLine(title, this.getGc().getUserWindowSize());
    return this.formatTitle(title) + "\n";
  }

  protected boolean isCaseSensitive() {
    return isCaseSensitive;
  }

  protected void setCaseSensitive(boolean isCaseSensitive) {
    this.isCaseSensitive = isCaseSensitive;
  }

  public boolean isIupac() {
    return isIupac;
  }

  /**
   * NB: setIupac() must be called BEFORE setSeqRegex() in order to have effect! See unit test.
   * Thius should be changed!!.
   */
  public void setIupac(boolean isIupac) {
    this.isIupac = isIupac;
  }
}
