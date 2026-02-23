package tracks;

import colouring.Config;
import colouring.ConfigKey;
import exceptions.InvalidColourException;
import exceptions.InvalidGenomicCoordsException;
import exceptions.InvalidRecordException;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.regex.Pattern;
import java.util.stream.Stream;
import org.apache.commons.lang3.NotImplementedException;
import org.apache.commons.lang3.StringUtils;
import samTextViewer.GenomicCoords;
import samTextViewer.Utils;

public class TrackReads extends AbstractTrack {

  private List<List<SamSequenceFragment>> readStack;
  private List<SAMRecord> readsPassingFilter;
  // private boolean withReadName= false;
  private long nRecsInWindow = -1;
  private int userWindowSize;
  private List<Argument> colourForRegex = null;
  private long alnRecCnt = -1;

  /* C o n s t r u c t o r s */

  public TrackReads(String bam, GenomicCoords gc)
      throws IOException,
          InvalidGenomicCoordsException,
          ClassNotFoundException,
          InvalidRecordException,
          SQLException {

    this.setTrackFormat(TrackFormat.BAM);

    if (!Utils.bamHasIndex(bam)) {
      File temp = Utils.createTempFile(".asciigenome.", ".bam", true);
      String fasta = null;
      if (gc != null) {
        fasta = gc.getFastaFile();
      }
      Utils.sortAndIndexSamOrBam(bam, temp.getAbsolutePath(), true, fasta);
      this.setWorkFilename(temp.getAbsolutePath());
    } else {
      this.setWorkFilename(bam);
    }
    this.setFilename(bam);
    this.setGc(gc);
  }

  /* M e t h o d s */

  @Override
  public void close() {}

  public void update()
      throws InvalidGenomicCoordsException,
          IOException,
          SQLException,
          InvalidRecordException,
          ClassNotFoundException {

    if (this.getyMaxLines() == 0) {
      return;
    }

    this.userWindowSize = this.getGc().getUserWindowSize();

    if (this.getGc().getGenomicWindowSize() < this.MAX_REGION_SIZE) {
      try (SamReader samReader =
              Utils.getSamReader(this.getWorkFilename(), this.getGc().getFastaFile());
          Stream<SAMRecord> passFilter =
              this.filterReads(
                  samReader,
                  this.getGc().getChrom(),
                  this.getGc().getFrom(),
                  this.getGc().getTo())) {
        List<TextRead> textReads = this.reservoirSampling(passFilter);
        this.readStack = stackReads(textReads);
      }
    } else {
      this.readStack = new ArrayList<>();
      this.nRecsInWindow = -1;
    }
  }

  private List<TextRead> reservoirSampling(Stream<SAMRecord> passFilter)
      throws InvalidGenomicCoordsException, IOException {
    float max_reads = Float.parseFloat(Config.get(ConfigKey.max_reads_in_stack));
    IndexedLinkedHashMap<String, List<TextRead>> selectedReads = new IndexedLinkedHashMap<>();

    Random rnd = new Random();
    Iterator<SAMRecord> iter = passFilter.iterator();
    this.nRecsInWindow = 0;
    int nTemplates = 0;
    while (iter.hasNext()) {
      this.nRecsInWindow += 1;
      SAMRecord rec = iter.next();
      String name = Utils.templateNameFromSamReadName(rec.getReadName());
      if (selectedReads.containsKey(name)) {
        // Already selected, just add the read
        selectedReads.get(name).add(new TextRead(rec, getGc(), getShowSoftClip()));
        continue;
      }
      nTemplates += 1;
      if (selectedReads.size() < max_reads) {
        // Fill up reservoir
        selectedReads.put(name, new ArrayList<>());
        selectedReads.get(name).add(new TextRead(rec, getGc(), getShowSoftClip()));
      } else {
        int j = rnd.nextInt(nTemplates + 1);
        if (j < max_reads) {
          String templateToRemove = selectedReads.getKeyAtIndex(j);
          selectedReads.remove(templateToRemove);
          selectedReads.put(name, new ArrayList<>());
          selectedReads.get(name).add(new TextRead(rec, getGc(), getShowSoftClip()));
        }
      }
    }
    List<TextRead> textReads = new ArrayList<>();
    for (List<TextRead> x : selectedReads.values()) {
      textReads.addAll(x);
    }
    return textReads;
  }

  @Override
  public String printToScreen() throws InvalidGenomicCoordsException, InvalidColourException {

    int yMaxLines = (this.getyMaxLines() < 0) ? Integer.MAX_VALUE : this.getyMaxLines();

    // If there are more lines (inner lists) than desired lines of output (yMaxLines), get a
    // representative sample
    List<Double> keep;
    if (this.readStack.isEmpty()) {
      return "";
    } else if (this.readStack.size() > yMaxLines) {
      keep = Utils.seqFromToLenOut(0, yMaxLines, yMaxLines);
    } else {
      keep = Utils.seqFromToLenOut(0, this.readStack.size() - 1, this.readStack.size());
    }
    StringBuilder printable = new StringBuilder();
    // this.changeFeatureColour(null);
    for (Double idx : keep) {
      List<SamSequenceFragment> line = this.readStack.get((int) Math.rint(idx));
      try {
        printable.append(linePrinter(line, this.bisulf, this.isNoFormat()));
        printable.append("\n");
      } catch (IOException e) {
        e.printStackTrace();
      }
    }
    return printable.toString().replaceAll("\n$", "");
  }

  /**
   * Put in the same list reads that will go in the same line of text Example Input, a list of
   * TextRead's: AAAAAAAAAAAA CCCCCCCCCCCC TTTTTTTTTTT GGGGGGGGGGG AAAAAAAA
   *
   * <p>Output, each line is a list of TextRead: [AAAAAAAAAAAA TTTTTTTTTTT GGGGGGGGGGG] [
   * CCCCCCCCCCCC AAAAAAAA]
   *
   * @throws IOException
   * @throws InvalidGenomicCoordsException
   */
  private List<List<SamSequenceFragment>> stackReads(List<TextRead> textReads)
      throws InvalidGenomicCoordsException, IOException {

    List<List<SamSequenceFragment>> listOfLines = new ArrayList<List<SamSequenceFragment>>();
    if (textReads.size() == 0) {
      return listOfLines;
    }

    List<SamSequenceFragment> fragments = this.makeFragments(textReads, this.getReadsAsPairs());
    List<SamSequenceFragment> line =
        new ArrayList<SamSequenceFragment>(); // All fragments going in the same line
    line.add(fragments.get(0));
    fragments.remove(0);
    listOfLines.add(line);
    int gap =
        (this.getGc().isSingleBaseResolution())
            ? 1
            : 0; // If reads are very compressed, do not add space between adjacent ones.
    while (true) {
      ArrayList<SamSequenceFragment> fragToRemove = new ArrayList<SamSequenceFragment>();
      // Find a fragment in input whose start is greater then end of current
      for (int i = 0; i < fragments.size(); i++) {
        SamSequenceFragment frag = fragments.get(i);
        if (frag.getTextStart()
            > line.get(line.size() - 1).getTextEnd()
                + gap) { // +2 because we want some space between adjacent reads
          listOfLines.get(listOfLines.size() - 1).add(frag); // Append to the last line.
          fragToRemove.add(frag);
        }
      } // At the end of the loop you have put in line as many reads as you can.
      for (SamSequenceFragment tr : fragToRemove) {
        fragments.remove(fragments.indexOf(tr));
      }
      // Create a new line, add the first textRead in list
      if (fragments.size() > 0) {
        line = new ArrayList<SamSequenceFragment>();
        line.add(fragments.get(0));
        listOfLines.add(line);
        fragments.remove(0);
      } else {
        break;
      }
    }
    return listOfLines;
  }

  /**
   * Match reads in textReads list to return a list fragments. Fragments are returned sorted by
   * start position.
   */
  private List<SamSequenceFragment> makeFragments(List<TextRead> textReads, boolean asPair) {

    List<SamSequenceFragment> fragments = new ArrayList<SamSequenceFragment>();

    while (!textReads.isEmpty()) {
      // Keep going until all reads have been moved to the list of fragments.
      TextRead tr = textReads.get(0);
      textReads.remove(0);
      if (!asPair
          || !(tr.getSamRecord().getReadPairedFlag() && tr.getSamRecord().getProperPairFlag())) {
        SamSequenceFragment frag = new SamSequenceFragment(tr);
        if (!asPair) {
          frag.setSingleton(true);
        }
        fragments.add(frag);
      } else {
        // Find the mate of this read, if present.
        TextRead mate = null;
        for (TextRead candidateMate : textReads) {
          if ((candidateMate.getSamRecord().getReadPairedFlag()
                  && candidateMate.getSamRecord().getProperPairFlag())
              && Utils.equalReadNames(
                  tr.getSamRecord().getReadName(), candidateMate.getSamRecord().getReadName())
              && tr.getSamRecord().getAlignmentStart()
                  == candidateMate.getSamRecord().getMateAlignmentStart()) {
            mate = candidateMate;
            break;
          }
        } // After this loop either we have found a mate or not. Either way, create a fragment from
        // a singleton or a pair.
        if (mate == null) {
          fragments.add(new SamSequenceFragment(tr));
        } else {
          fragments.add(new SamSequenceFragment(tr, mate));
          textReads.remove(mate);
        }
      }
    }
    sortFragmentsByStartPosition(fragments); // This may be redundant.
    return fragments;
  }

  private void sortFragmentsByStartPosition(List<SamSequenceFragment> fragments) {
    Collections.sort(
        fragments,
        new Comparator<SamSequenceFragment>() {
          @Override
          public int compare(final SamSequenceFragment frag1, final SamSequenceFragment frag2) {
            return Integer.compare(
                frag1.getLeftRead().getAlignmentStart(), frag2.getLeftRead().getAlignmentStart());
          }
        });
  }

  /** Prepare a printable string of each output line. */
  private String linePrinter(List<SamSequenceFragment> fragments, boolean bs, boolean noFormat)
      throws IOException, InvalidGenomicCoordsException, InvalidColourException {
    StringBuilder sb = new StringBuilder();

    int curPos = 0; // Position on the line, needed to pad with blanks btw reads.
    for (SamSequenceFragment frag : fragments) {
      String line =
          StringUtils.repeat(" ", (frag.getTextStart() - 1) - curPos); // Left pad with spaces
      sb.append(line);
      String printableRead = frag.getPrintableFragment(bs, noFormat);
      sb.append(printableRead);
      curPos = frag.getTextEnd();
    }
    int nchars = Utils.stripAnsiCodes(sb.toString()).length();
    if (nchars < this.userWindowSize) {
      sb.append(StringUtils.repeat(' ', this.userWindowSize - nchars));
    }
    return this.highlightMidCharacter(sb.toString());
  }

  /** Find the mid character and add some formatting to highlight it. */
  private String highlightMidCharacter(String fmtLine) {
    if (this.isNoFormat() || !Utils.asBoolean((Config.get(ConfigKey.highlight_mid_char)))) {
      return fmtLine;
    }
    List<Integer> idx = Utils.indexOfCharsOnFormattedLine(fmtLine);
    if (idx.size() <= 7) {
      return fmtLine;
    }
    int mid = idx.get(this.userWindowSize / 2);
    char midChar = fmtLine.charAt(mid);
    if (midChar != ' ') {
      String hLine =
          fmtLine.substring(0, mid)
              + "\033[1;7m"
              + midChar
              + "\033[21;27m"
              + fmtLine.substring(mid + 1);
      return hLine;
    }
    return fmtLine;
  }

  @Override
  public String getTitle()
      throws InvalidColourException, InvalidGenomicCoordsException, IOException {

    if (this.isHideTitle()) {
      return "";
    }

    String libsize = "";
    if (this.alnRecCnt != -1) {
      libsize = "/" + this.alnRecCnt;
    }
    String xtitle =
        this.getTrackTag()
            + "; Reads: "
            + this.nRecsInWindow
            + libsize
            + this.getTitleForActiveFilters();
    return this.formatTitle(xtitle) + "\n";
  }

  @Override
  protected String getTitleForActiveFilters() {
    List<String> title = new ArrayList<String>();
    if (!this.getAwk().equals(Filter.DEFAULT_AWK.getValue())) {
      title.add("awk");
    }
    if (!this.getSystemCommand().equals(Filter.DEFAULT_SYSTEM_COMMAND.getValue())) {
      title.add("stream");
    }
    if (!this.getShowRegex().pattern().equals(Filter.DEFAULT_SHOW_REGEX.getValue())
        || !this.getHideRegex().pattern().equals(Filter.DEFAULT_HIDE_REGEX.getValue())) {
      title.add("grep");
    }
    if (this.get_f_flag() != Integer.parseInt(Filter.DEFAULT_f_FLAG.getValue())
        || this.get_F_flag() != Integer.parseInt(Filter.DEFAULT_F_FLAG.getValue())) {
      title.add("bit-flag");
    }
    if (this.getMapq() != Integer.parseInt(Filter.DEFAULT_MAPQ.getValue())) {
      title.add("mapq");
    }
    if (!this.getFeatureFilter()
        .getVariantChrom()
        .equals(Filter.DEFAULT_VARIANT_CHROM.getValue())) {
      title.add("var-read");
    }
    if (!title.isEmpty()) {
      return "; filters: " + title;
    } else {
      return "";
    }
  }

  @Override
  protected List<String> getRecordsAsStrings() {
    List<String> featureList = new ArrayList<String>();

    for (List<SamSequenceFragment> x : this.readStack) {
      for (SamSequenceFragment frag : x) {
        featureList.add(frag.getLeftRead().getSamRecord().getSAMString());
        if (frag.getRightRead() != null) {
          featureList.add(frag.getRightRead().getSamRecord().getSAMString());
        }
      }
    }
    return featureList;
  }

  /* S e t t e r s   and   G e t t e r s */

  @Override
  public void setReadsAsPairs(boolean readsAsPairs)
      throws InvalidGenomicCoordsException,
          IOException,
          SQLException,
          InvalidRecordException,
          ClassNotFoundException {
    this.readsAsPairs = readsAsPairs;
    this.update();
  }

  @Override
  public void setShowSoftClip(boolean showSoftClip)
      throws InvalidGenomicCoordsException,
          IOException,
          SQLException,
          InvalidRecordException,
          ClassNotFoundException {
    this.showSoftClip = showSoftClip;
    this.update();
  }

  @Override
  public void changeFeatureColour(List<Argument> args) {
    List<List<SamSequenceFragment>> stack = this.readStack;
    for (List<SamSequenceFragment> frags : stack) {
      for (SamSequenceFragment frag : frags) {
        List<TextRead> reads = new ArrayList<TextRead>();
        reads.add(frag.getLeftRead());
        if (frag.getRightRead() != null) {
          reads.add(frag.getRightRead());
        }
        for (TextRead tr : reads) {
          for (Argument arg : args) {
            String regex = arg.getKey();
            // String colour= arg.getArg();
            boolean matched =
                Pattern.compile(regex).matcher(tr.getSamRecord().getSAMString()).find();
            if (arg.isInvert()) {
              matched = !matched;
            }
            if (matched) {
              // tr.getTextReadAsFeatureChars(this.isBisulf());
              // f.setFgColour(colour);
            }
          }
        }
      }
    }
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

  //	private List<Argument> getColourForRegex() {
  //		return this.colourForRegex;
  //	}

  @Override
  public void reload()
      throws InvalidGenomicCoordsException,
          IOException,
          ClassNotFoundException,
          InvalidRecordException,
          SQLException {
    if (!Files.isSameFile(Paths.get(this.getWorkFilename()), Paths.get(this.getFilename()))) {
      TrackReads tr = new TrackReads(this.getFilename(), this.getGc());
      String fname = this.getWorkFilename();
      Files.move(
          Paths.get(tr.getWorkFilename()),
          Paths.get(fname),
          java.nio.file.StandardCopyOption.REPLACE_EXISTING);
      Files.move(
          Paths.get(tr.getWorkFilename().replaceAll("\\.bam$", ".bai")),
          Paths.get(fname.replaceAll("\\.bam$", ".bai")),
          java.nio.file.StandardCopyOption.REPLACE_EXISTING);
    }
    this.update();
  }

  @Override
  public void setFeatureName(String gtfAttributeForName) {}

  @Override
  public ArrayList<String> getChromosomeNames() {
    List<SAMSequenceRecord> samSeqRecs = this.getGc().getSamSeqDict().getSequences();
    ArrayList<String> chromosomeNames = new ArrayList<String>();
    for (SAMSequenceRecord x : samSeqRecs) {
      chromosomeNames.add(x.getContig());
    }
    return chromosomeNames;
  }

  /**
   * A LinkedHashMap that can be queried by the index position of a key credit:
   * https://stackoverflow.com/a/57550038/1114453
   */
  class IndexedLinkedHashMap<K, V> extends LinkedHashMap<K, V> {

    private static final long serialVersionUID = 1L;

    ArrayList<K> al_Index = new ArrayList<K>();

    @Override
    public V put(K key, V val) {
      if (!super.containsKey(key)) al_Index.add(key);
      V returnValue = super.put(key, val);
      return returnValue;
    }

    public V getValueAtIndex(int i) {
      return (V) super.get(al_Index.get(i));
    }

    public K getKeyAtIndex(int i) {
      return (K) al_Index.get(i);
    }

    public int getIndexOf(K key) {
      return al_Index.indexOf(key);
    }

    @Override
    public V remove(Object key) {
      V v = super.remove(key);
      al_Index.remove(key);
      return v;
    }

    @Override
    public void clear() {
      super.clear();
      al_Index.clear();
    }

    @Override
    public void putAll(Map<? extends K, ? extends V> m) {
      throw new NotImplementedException("");
    }
  }
}
