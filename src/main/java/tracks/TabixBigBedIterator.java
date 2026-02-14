package tracks;

import com.google.common.base.Splitter;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.tribble.readers.TabixReader;
import htsjdk.tribble.readers.TabixReader.Iterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFReader;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.EnumSet;
import java.util.List;
import java.util.NoSuchElementException;
import org.broad.igv.bbfile.BBFileReader;
import org.broad.igv.bbfile.BedFeature;
import org.broad.igv.bbfile.BigBedIterator;

public class TabixBigBedIterator {

  Iterator tabixIterator;
  BigBedIterator bigBedIterator;
  CloseableIterator<VariantContext> vcfIterator;
  private VCFHeader vcfHeader;

  protected TabixBigBedIterator(TabixReader reader, String chrom, int start, int end) {
    this.tabixIterator = reader.query(chrom, start, end);
  }

  protected TabixBigBedIterator(BBFileReader reader, String chrom, int start, int end) {
    this.bigBedIterator = reader.getBigBedIterator(chrom, start, chrom, end, false);
  }

  protected TabixBigBedIterator(VCFReader reader, String chrom, int start, int end) {
    this.vcfIterator = reader.query(chrom, start, end);
  }

  private String variantToVcfLine(VariantContext ctx, VCFHeader header) {
    ByteArrayOutputStream baos = new ByteArrayOutputStream();

    VariantContextWriter writer =
        new VariantContextWriterBuilder()
            .setOutputVCFStream(baos)
            .setOptions(EnumSet.of(Options.ALLOW_MISSING_FIELDS_IN_HEADER))
            .build();
    writer.writeHeader(header);
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
      BedFeature x = this.bigBedIterator.next();
      if (x == null) {
        return null;
      }
      StringBuilder sb = new StringBuilder();
      sb.append(x.getChromosome());
      sb.append("\t");
      sb.append(x.getStartBase());
      sb.append("\t");
      sb.append(x.getEndBase());
      for (String field : x.getRestOfFields()) {
        sb.append("\t");
        sb.append(field);
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
