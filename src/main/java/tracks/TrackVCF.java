package tracks;

import com.google.common.base.Joiner;
import com.google.common.base.Splitter;
import exceptions.InvalidColourException;
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
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Pattern;

import org.apache.commons.lang3.StringUtils;
import samTextViewer.GenomicCoords;
import samTextViewer.Utils;
import sortBgzipIndex.MakeTabixIndex;

public class TrackVCF extends AbstractTrackIntervalFeature<VCFFeature> {

  protected List<VCFFeature> vcfFeatureList = new ArrayList<>();
  private VCFCodec vcfCodec;

  public TrackVCF(final String filename, GenomicCoords gc)
          throws IOException,
          InvalidGenomicCoordsException,
          ClassNotFoundException,
          InvalidRecordException,
          SQLException {

    this.setFilename(filename);
    this.setTrackFormat(TrackFormat.VCF);

    if (!Utils.hasTabixIndex(filename)) {
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

  @Override
  protected VCFFeature createFeature(String line)
          throws InvalidGenomicCoordsException {

    return new VCFFeature(
            line,
            getVCFCodec()
    );
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


  @Override
  protected VCFFeature mergeFeatures(
          VCFFeature a,
          VCFFeature b,
          boolean screenCoords)
          throws InvalidGenomicCoordsException, InvalidColourException {

    return new VCFFeature(a, b, screenCoords);
  }

  protected IntervalFeature findNextRegexInGenome(Pattern pattern, String chrom, int from)
          throws IOException, InvalidGenomicCoordsException {
    return super.findNextRegexInGenome(pattern, chrom, from);
  }

  @Override
  protected List<IntervalFeature> getFeaturesInInterval(String chrom, int from, int to)
      throws IOException, InvalidGenomicCoordsException {
    if (from < 1) {
      System.err.println("from < 1: " + from + "; resetting to 1.");
      from = 1;
    }
    if (from > to) {
      System.err.println(
          "Invalid coordinates: from: "
              + from
              + "; to: "
              + to
              + "; Resetting to initial 1-"
              + Integer.MAX_VALUE);
      throw new InvalidGenomicCoordsException();
    }
    List<IntervalFeature> xFeatures = this.getFeaturesInVCFInterval(chrom, from, to);
    this.removeInvisibleFeatures(xFeatures);
    return xFeatures;
  }

  private VCFCodec getVCFCodec() {
    if (this.vcfCodec == null) {
      VCFCodec vcfCodec = new VCFCodec();
      vcfCodec.setVCFHeader(this.getVcfHeader(), Utils.getVCFHeaderVersion(this.getVcfHeader()));
      this.vcfCodec = vcfCodec;
    }
    return this.vcfCodec;
  }



  private List<VCFFeature> getFeaturesInVCFInterval(String chrom, int from, int to)
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
    List<VCFFeature> xFeatures = new ArrayList<>();
    TabixBigBedIterator qry = this.getReader().query(chrom, from - 1, to);
    while (true) {
      String q = qry.next();
      if (q == null) {
        break;
      }
      VCFFeature intervalFeature =
          new VCFFeature(q, this.getVCFCodec(), this.getScoreColIdx());
      xFeatures.add(intervalFeature);
    }
    return xFeatures;
  }

  @Override
  public String printToScreen() throws InvalidGenomicCoordsException, InvalidColourException {
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

    // Genotype matrix
    try {
      String gtm =
              this.getGenotypeMatrix()
                      .printToScreen(
                              this.isNoFormat(),
                              this.vcfFeatureList,
                              this.getGc().getUserWindowSize(),
                              this.getVcfHeader());
      printable.add(gtm);
    } catch (InvalidColourException | IOException e) {
      e.printStackTrace();
    }
    return StringUtils.join(printable, "\n").replaceAll("\n$", "");
  }

  @Override
  protected List<String> getRecordsAsStrings() {
    List<String> featureList = new ArrayList<>();
    for (VCFFeature ift : this.intervalFeatureList) {
      if (this.getPrintNormalizedVcf()) {
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
    if (sampleNames.isEmpty()) {
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
}
