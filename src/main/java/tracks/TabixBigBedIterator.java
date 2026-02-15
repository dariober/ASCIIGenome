package tracks;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.tribble.readers.TabixReader;
import htsjdk.tribble.readers.TabixReader.Iterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFReader;
import java.io.IOException;
import java.util.NoSuchElementException;
import org.broad.igv.bbfile.BBFileReader;
import org.broad.igv.bbfile.BedFeature;
import org.broad.igv.bbfile.BigBedIterator;
import samTextViewer.Utils;

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
      return Utils.variantToVcfLine(ctx, this.vcfHeader);
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
