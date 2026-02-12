package tracks;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.tribble.readers.TabixReader;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import htsjdk.variant.vcf.VCFReader;
import org.broad.igv.bbfile.BBFileReader;

/** Adapter to make tabixReader and bigBed behave in the same way. */
public class TabixBigBedReader {

  private TabixReader tabixReader;
  private BBFileReader bigBedReader;
  private VCFReader vcfReader;

  protected TabixBigBedReader(TabixReader tabixReader) {
    this.tabixReader = tabixReader;
  }

  protected TabixBigBedReader(BBFileReader bigBedReader) {
    this.bigBedReader = bigBedReader;
  }

  protected TabixBigBedReader(VCFReader vcfReader) {
    this.vcfReader = vcfReader;
  }
  ;
  protected TabixBigBedIterator query(String chrom, int start, int end) {

    if (this.tabixReader != null) {
      return new TabixBigBedIterator(this.tabixReader, chrom, start, end);
    } else if (this.bigBedReader != null) {
      return new TabixBigBedIterator(this.bigBedReader, chrom, start, end);
    } else if (this.vcfReader != null) {
      return new TabixBigBedIterator(this.vcfReader, chrom, start, end);
    } else {
      throw new RuntimeException();
    }
  }

  public Set<String> getChromosomes() {
    if (this.tabixReader != null) {
      return this.tabixReader.getChromosomes();
    } else if (this.bigBedReader != null) {
      return new HashSet<>(this.bigBedReader.getChromosomeNames());
    } else if (this.vcfReader != null) {
      System.err.println(this.vcfReader.getHeader());
      List<SAMSequenceRecord> seqs = this.vcfReader.getHeader().getSequenceDictionary().getSequences();
      HashSet<String> chroms = new HashSet<>();
      seqs.forEach(x -> chroms.add(x.getSequenceName()));
      return chroms;
    } else {
      throw new RuntimeException("There is no reader available to get chromosome names");
    }
  }
}
