package tracks;

import com.google.common.base.Splitter;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.tribble.readers.TabixReader.Iterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderVersion;
import htsjdk.variant.vcf.VCFReader;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.EnumSet;
import java.util.List;
import java.util.NoSuchElementException;

import org.igv.feature.IGVFeature;
import org.igv.ucsc.bb.BBFeatureSource;
import org.igv.ucsc.bb.BBFile;
import utils.FlexibleTabixReader;

public class TabixBigBedIterator {

  Iterator tabixIterator;
  java.util.Iterator<IGVFeature> bigBedIterator;
  CloseableIterator<VariantContext> vcfIterator;
  private VCFHeader vcfHeader;
  private String chrom = null;

  protected TabixBigBedIterator(FlexibleTabixReader reader, String chrom, int start, int end) {
    this.chrom = chrom;
    this.tabixIterator = reader.query(chrom, start, end);
  }

  protected TabixBigBedIterator(BBFile reader, String chrom, int start, int end) throws IOException {
    BBFeatureSource bigbedFeatureSource = new BBFeatureSource(reader, null);
    this.bigBedIterator = bigbedFeatureSource.getFeatures(chrom, start, end);
    this.chrom = chrom;
  }

  protected TabixBigBedIterator(VCFReader reader, String chrom, int start, int end) {
    this.vcfIterator = reader.query(chrom, start, end);
    this.chrom = chrom;
  }

  private String variantToVcfLine(VariantContext ctx, VCFHeader header) {
    VCFHeader workHeader;
    if (header.getVCFHeaderVersion() != null
        && header.getVCFHeaderVersion().isAtLeastAsRecentAs(VCFHeaderVersion.VCF3_3)) {
      workHeader = new VCFHeader(header.getMetaDataInInputOrder(), header.getGenotypeSamples());
      workHeader.setVCFHeaderVersion(VCFHeaderVersion.VCF4_2);
    } else {
      workHeader = header;
    }
    ByteArrayOutputStream baos = new ByteArrayOutputStream();

    VariantContextWriter writer =
        new VariantContextWriterBuilder()
            .setOutputVCFStream(baos)
            .setOptions(EnumSet.of(Options.ALLOW_MISSING_FIELDS_IN_HEADER))
            .build();
    writer.writeHeader(workHeader);
    writer.add(ctx);
    writer.close();

    List<String> vcfLinesWithHeader =
        Splitter.on("\n").splitToList(baos.toString(StandardCharsets.UTF_8));
    List<String> vcfLines = new ArrayList<>();
    for (String line : vcfLinesWithHeader) {
      if (!line.startsWith("#") && !line.trim().isEmpty()) {
        vcfLines.add(line.trim());
      }
    }
    assert vcfLines.size() == 1;
    return vcfLines.get(0);
  }

  protected void setVcfHeader(VCFHeader vcfHeader) {
    this.vcfHeader = vcfHeader;
  }

  protected String next() throws IOException {

    if (this.tabixIterator != null) {
      return this.tabixIterator.next();

    } else if (this.bigBedIterator != null) {
      if (!this.bigBedIterator.hasNext()) {
        return null;
      }
      IGVFeature x = this.bigBedIterator.next();
      if (x == null) {
        return null;
      }
      StringBuilder sb = new StringBuilder();
      sb.append(this.chrom);
      sb.append("\t");
      sb.append(x.getStart());
      sb.append("\t");
      sb.append(x.getEnd());
      for (String key : x.getAttributeKeys()) {
        sb.append("\t");
        sb.append(x.getAttribute(key));
      }
      return sb.toString();
    } else if (this.vcfIterator != null) {
      VariantContext ctx;
      try {
        ctx = vcfIterator.next();
      } catch (NoSuchElementException e) {
        this.vcfIterator.close();
        return null;
      }
      return this.variantToVcfLine(ctx, this.vcfHeader);
    } else {
      throw new RuntimeException();
    }
  }

  protected void close() {
    if (this.vcfIterator != null) {
      this.vcfIterator.close();
    }
  }
}
