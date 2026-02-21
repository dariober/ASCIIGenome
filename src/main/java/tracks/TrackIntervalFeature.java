package tracks;

import exceptions.InvalidColourException;
import exceptions.InvalidGenomicCoordsException;
import exceptions.InvalidRecordException;
import htsjdk.samtools.util.FileExtensions;
import htsjdk.tribble.readers.TabixReader;
import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.*;
import java.util.stream.Stream;

import org.apache.commons.lang3.StringUtils;
import org.broad.igv.bbfile.BBFileReader;
import samTextViewer.GenomicCoords;
import samTextViewer.Utils;
import sortBgzipIndex.MakeTabixIndex;

public class TrackIntervalFeature extends AbstractTrackFeature<IntervalFeature> {

  protected int scoreColIdx = -1;

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

  protected TrackIntervalFeature() {}

  /* M e t h o d s */
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
              Integer.parseInt(nameFieldOrAttribute)
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

  @Override
  protected IntervalFeature createFeature(String line) throws InvalidGenomicCoordsException {
    return new IntervalFeature(line, this.getTrackFormat(), this.getScoreColIdx());
  }

  @Override
  public String printToScreen() throws InvalidGenomicCoordsException, InvalidColourException {
    for (IntervalFeature x : this.getFeatureList()) {
      if (this.getTrackFormat().equals(TrackFormat.GFF)
          || this.getTrackFormat().equals(TrackFormat.GTF)) {
        x.setGtfAttributeForName(this.gtfAttributeForName);
      } else if (this.getTrackFormat().equals(TrackFormat.BED)
          || this.getTrackFormat().equals(TrackFormat.BIGBED)) {
        x.setBedFieldName(this.bedFieldForName);
      }
    }
    List<String> printable = new ArrayList<>();
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
    return StringUtils.join(printable, "\n").replaceAll("\n$", "");
  }

  @Override
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
          new IntervalFeature(q, this.getTrackFormat(), this.getScoreColIdx());
      xFeatures.add(intervalFeature);
    }
    this.removeInvisibleFeatures(xFeatures);
    return xFeatures;
  }

  /**
   * Group the features in this genomic window by GFF attribute (typically a transcripts). Features
   * that don't have the attribute make each a length=1 list.
   */
  @Override
  protected Map<String, List<IntervalFeature>> groupByGFFAttribute() {

    // * First collect the IDs of the transcripts

    // Key is transcript ID e.g. ENST00001234.
    // Values is all the IntervalFeatures captured by this ID and part of a transcript.
    // I.e. their are in txFeature set,
    Map<String, List<IntervalFeature>> txIds = new LinkedHashMap<>();

    // This key:value is for records which are not part transcripts. E.g. features like "chromosome"
    // or rRNA.
    txIds.put("_na_", new ArrayList<>());

    // Now populate the lists of values by assigning to each key the transcript records:
    for (IntervalFeature x : this.getFeatureList()) {

      if (FormatGTF.getTxSuperFeatures().contains(x.getFeature().toLowerCase())) {

        // Transcript feature. E.g.
        // 7 ensembl_havana mRNA 5527151 5530709 . - .
        // ID=transcript:ENST00000331789;Parent=gene:ENSG00000075624;Name=ACTB-...
        String txId = x.getGFFValueFromKey("ID");
        if (!txIds.containsKey(txId)) {
          txIds.put(txId, new ArrayList<>());
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
          txIds.put(txId, new ArrayList<>());
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
  protected Map<String, List<IntervalFeature>> groupByGTFAttribute() {

    // * First collect the IDs of the transcripts

    // Key is transcript ID e.g. ENST00001234.
    // Values is all the IntervalFeatures captured by this ID and part of a transcript.
    // I.e. their are in txFeature set,
    Map<String, List<IntervalFeature>> txIds = new LinkedHashMap<>();

    // This key:value is for records which are not part transcripts. E.g. features like "chromosome"
    // or rRNA.
    txIds.put("_na_", new ArrayList<>());

    // Now populate the lists of values by assigning to each key the transcript records:
    for (IntervalFeature x : this.getFeatureList()) {

      if (FormatGTF.getTxSuperFeatures().contains(x.getFeature().toLowerCase())
          || FormatGTF.getTxSubFeatures().contains(x.getFeature().toLowerCase())) {
        // Transcript feature. E.g.
        // chr7 hg19_wgEncodeGencodeBasicV19 exon       5566782 5567522 0.000000 - . gene_id
        // "ENST00000331789.5"; transcript_id "ENST00000331789.5";
        String txId = x.getGFFValueFromKey("transcript_id");
        if (!txIds.containsKey(txId)) {
          txIds.put(txId, new ArrayList<>());
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
   */
  protected IntervalFeature collapseGFFTranscript(
      List<IntervalFeature> txFeatures, List<Double> mapToScreen)
      throws InvalidGenomicCoordsException, InvalidColourException {

    if (txFeatures.isEmpty()) {
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

  protected List<IntervalFeature> getFeatureList() {
    return this.featureList;
  }

  @Override
  public void streamFeaturesThroughSystemCommand() throws IOException, InvalidGenomicCoordsException {

    Stream<String> textLines = this.getFeatureList()
            .stream()
            .map(IntervalFeature::getRaw);

    List<IntervalFeature> features = new ArrayList<>();

    try (Stream<String> out = Utils.streamLinesThroughSystemCommand(textLines, null, this.getSystemCommand())) {
      Iterator<String> iter = out.iterator();
      while (iter.hasNext()) {
        features.add(new IntervalFeature(iter.next(), this.getTrackFormat(), this.getScoreColIdx()));
      }
    }
    this.setFeatureList(features);
  };
}
