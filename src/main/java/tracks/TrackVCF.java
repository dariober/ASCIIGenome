package tracks;

import static htsjdk.variant.vcf.VCFFileReader.isBCF;

import com.github.lindenb.jvarkit.variant.bcf.BCFFileReader;
import com.google.common.base.Joiner;
import com.google.common.base.Splitter;
import exceptions.InvalidColourException;
import exceptions.InvalidGenomicCoordsException;
import exceptions.InvalidRecordException;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.FileExtensions;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.tribble.readers.TabixReader;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFReader;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Path;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.EnumSet;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import org.apache.commons.lang3.NotImplementedException;
import org.apache.commons.lang3.StringUtils;
import samTextViewer.GenomicCoords;
import samTextViewer.Utils;
import sortBgzipIndex.MakeTabixIndex;

public class TrackVCF extends AbstractTrackFeature<VCFFeature> {

  private List<VCFFeature> featureList = new ArrayList<>();
  private VCFCodec vcfCodec;
  private final VCFReader vcfReader;
  private final boolean isBCF;
  private GenotypeMatrix genotypeMatrix = new GenotypeMatrix();

  public TrackVCF(final String filename, GenomicCoords gc)
      throws IOException,
          InvalidGenomicCoordsException,
          ClassNotFoundException,
          InvalidRecordException,
          SQLException {

    this.isBCF = isBCF(Path.of(filename));

    this.setFilename(filename);
    this.setTrackFormat(TrackFormat.VCF);
    if (!Utils.hasTabixIndex(filename) && !isBCF(Path.of(filename))) {
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
    } else if (isBCF(Path.of(filename)) && !new File(Path.of(filename) + ".csi").exists()) {
      throw new NotImplementedException("Cannot sort and index bcf files.");
    } else { // This means the input is indexed.
      this.setWorkFilename(filename);
    }
    this.vcfReader = this.prepareVcfReader(this.getWorkFilename());
    if (!isBCF) {
      this.setTabixReader(new TabixReader(this.getWorkFilename()));
    } else {
      this.tabixReader = null;
    }
    this.setGc(gc);
  }

  protected GenotypeMatrix getGenotypeMatrix() {
    return genotypeMatrix;
  }

  @Override
  protected VCFFeature createFeature(String line) {
    return new VCFFeature(line, getVCFCodec());
  }

  @Override
  protected TrackFormat getTrackFormat() {
    return TrackFormat.VCF;
  }

  @Override
  public void setGc(GenomicCoords gc)
      throws ClassNotFoundException,
          IOException,
          InvalidGenomicCoordsException,
          InvalidRecordException,
          SQLException {
    this.gc = gc;
    this.update();
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
    this.featureList =
        this.getFeaturesInInterval(
            this.getGc().getChrom(), this.getGc().getFrom(), this.getGc().getTo());
    for (VCFFeature ift : this.getFeatureList()) {
      ift.mapToScreen(this.getGc().getMapping());
    }
  }

  private VCFCodec getVCFCodec() {
    if (this.vcfCodec == null) {
      VCFCodec vcfCodec = new VCFCodec();
      vcfCodec.setVCFHeader(this.getVcfHeader(), Utils.getVCFHeaderVersion(this.getVcfHeader()));
      this.vcfCodec = vcfCodec;
    }
    return this.vcfCodec;
  }

  @Override
  public VCFHeader getVcfHeader() {
    return this.vcfReader.getHeader();
  }

  @Override
  public void streamFeaturesThroughSystemCommand() {
    throw new NotImplementedException("TODO");
  }

  protected List<VCFFeature> getFeaturesInInterval(String chrom, int from, int to)
      throws IOException {

    // Collect feature
    List<VCFFeature> xFeatures = new ArrayList<>();
    List<VariantContext> ctxList = new ArrayList<>();

    try (CloseableIterator<VariantContext> iter = this.vcfReader.query(chrom, from, to)) {
      while (iter.hasNext()) {
        ctxList.add(iter.next());
      }
    }
    for (String vcfLine : variantsToVcfLines(ctxList, this.getVcfHeader())) {
      VCFFeature vcfFeature = new VCFFeature(vcfLine, this.getVCFCodec());
      xFeatures.add(vcfFeature);
    }
    this.filterFeaturesInPlace(xFeatures);
    return xFeatures;
  }

  private List<String> variantsToVcfLines(List<VariantContext> variants, VCFHeader header) {

    ByteArrayOutputStream baos = new ByteArrayOutputStream();

    VariantContextWriter writer =
        new VariantContextWriterBuilder()
            .setOutputVCFStream(baos)
            .setOptions(EnumSet.of(Options.ALLOW_MISSING_FIELDS_IN_HEADER))
            .build();
    writer.writeHeader(header);

    for (VariantContext vc : variants) {
      writer.add(vc);
    }

    writer.close();

    List<String> vcfLinesWithHeader =
        Splitter.on("\n").splitToList(baos.toString(StandardCharsets.UTF_8));
    List<String> vcfLines = new ArrayList<>();
    for (String line : vcfLinesWithHeader) {
      if (!line.startsWith("#") && !line.trim().isEmpty()) {
        vcfLines.add(line.trim());
      }
    }
    return vcfLines;
  }

  private VCFReader prepareVcfReader(String vcfFile) throws IOException {
    if (Utils.urlFileExists(vcfFile)) {
      AbstractFeatureReader<VariantContext, LineIterator> reader =
          AbstractFeatureReader.getFeatureReader(
              vcfFile, new VCFCodec(), true // requires index
              );
      return new VCFReaderAdapter(reader);
    }
    if (isBCF) {
      return new BCFFileReader(Path.of(vcfFile), true);
    } else {
      return new VCFFileReader(new File(vcfFile), true);
    }
  }

  @Override
  public String printToScreen() throws InvalidGenomicCoordsException, InvalidColourException {
    List<String> printable = new ArrayList<>();
    int nLines = 0;
    try {
      for (List<VCFFeature> listToPrint : this.stackFeatures()) {
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
                  this.getFeatureList(),
                  this.getGc().getUserWindowSize(),
                  this.getVcfHeader());
      printable.add(gtm);
    } catch (InvalidColourException | IOException e) {
      e.printStackTrace();
    }
    return StringUtils.join(printable, "\n").replaceAll("\n$", "");
  }

  @Override
  public void setFeatureName(String gtfAttributeForName) {}

  @Override
  protected List<String> getRecordsAsStrings() {
    List<String> featureList = new ArrayList<>();

    LinkedHashSet<String> samplesToKeep;
    try {
      samplesToKeep =
          this.getGenotypeMatrix().selectedSamples(this.getFeatureList(), this.getVcfHeader());
    } catch (IOException e) {
      throw new RuntimeException(e);
    }
    VCFHeader reHeader = this.copyHeaderWithSamples(samplesToKeep);

    for (VCFFeature vcfFeature : this.getFeatureList()) {
      VCFFeature reVcfFeature = (VCFFeature) vcfFeature.clone();

      reVcfFeature.stripVariantContext(Set.of(), Set.of(), new LinkedHashSet<>(samplesToKeep));
      reVcfFeature.setRaw(Utils.variantToVcfLine(reVcfFeature.getVariantContext(), reHeader));

      if (this.getPrintNormalizedVcf()) {
        List<String> line =
            this.normalizeVcfRecordBySample(
                reHeader.getSampleNamesInOrder(), reVcfFeature.getRaw());
        featureList.addAll(line);
      } else {
        featureList.add(reVcfFeature.getRaw());
      }
    }
    return featureList;
  }

  private VCFHeader copyHeaderWithSamples(LinkedHashSet<String> samplesToKeep) {
    Set<VCFHeaderLine> metaData = new HashSet<>(this.getVcfHeader().getMetaDataInInputOrder());
    Set<String> newSamples;

    if (samplesToKeep == null) {
      // keep all samples
      newSamples = new LinkedHashSet<>(this.getVcfHeader().getSampleNamesInOrder());
    } else {
      // keep only selected samples
      newSamples = new LinkedHashSet<>(samplesToKeep);
    }

    return new VCFHeader(metaData, newSamples);
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

  protected List<VCFFeature> getFeatureList() {
    return this.featureList;
  }

  @Override
  TabixBigBedReader getReader() {
    if (isBCF) {
      TabixBigBedReader tbb = new TabixBigBedReader(this.vcfReader);
      tbb.setVcfHeader(this.getVcfHeader());
      return tbb;
    }
    return new TabixBigBedReader(this.tabixReader);
  }

  // This is a bad but for now let's get on with it
  // --->
  @Override
  protected Map<String, List<VCFFeature>> groupByGFFAttribute() {
    return Map.of();
  }

  @Override
  protected Map<String, List<VCFFeature>> groupByGTFAttribute() {
    return Map.of();
  }

  @Override
  protected VCFFeature collapseGFFTranscript(List<VCFFeature> features, List<Double> mapToScreen) {
    return null;
  }

  // <----

  private static class VCFReaderAdapter implements VCFReader {
    private final AbstractFeatureReader<VariantContext, LineIterator> reader;

    public VCFReaderAdapter(AbstractFeatureReader<VariantContext, LineIterator> reader) {
      this.reader = reader;
    }

    @Override
    public VCFHeader getHeader() {
      return (VCFHeader) reader.getHeader();
    }

    @Override
    public CloseableIterator<VariantContext> query(String chrom, int start, int end) {
      try {
        return reader.query(chrom, start, end);
      } catch (IOException e) {
        throw new RuntimeException("Error querying VCF: " + e.getMessage(), e);
      }
    }

    @Override
    public boolean isQueryable() {
      return true;
    }

    @Override
    public CloseableIterator<VariantContext> iterator() {
      try {
        return reader.iterator();
      } catch (IOException e) {
        throw new RuntimeException(e);
      }
    }

    @Override
    public void close() throws IOException {
      reader.close();
    }
  }
}
