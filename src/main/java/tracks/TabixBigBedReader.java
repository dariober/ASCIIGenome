package tracks;

import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFReader;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.apache.commons.lang3.NotImplementedException;
import org.igv.ucsc.bb.BBFeatureSource;
import org.igv.ucsc.bb.BBFile;
import org.igv.ucsc.bb.ChromTree;
import utils.FlexibleTabixReader;

/** Adapter to make tabixReader and bigBed behave in the same way. */
public class TabixBigBedReader {

  private FlexibleTabixReader tabixReader;
  private BBFile bigBedReader;
  private VCFReader vcfReader;
  private VCFHeader vcfHeader;
  private String fileName;

  protected TabixBigBedReader(FlexibleTabixReader tabixReader) {
    this.tabixReader = tabixReader;
  }

  protected TabixBigBedReader(BBFile bigBedReader, String fileName) {
    this.fileName = fileName;
    this.bigBedReader = bigBedReader;
  }

  protected TabixBigBedReader(VCFReader vcfReader) {
    this.vcfReader = vcfReader;
  }

  protected TabixBigBedIterator query(String chrom, int start, int end) throws IOException {

    if (this.tabixReader != null) {
      return new TabixBigBedIterator(this.tabixReader, chrom, start, end);
    } else if (this.bigBedReader != null) {
      return new TabixBigBedIterator(this.bigBedReader, chrom, start, end);
    } else if (this.vcfReader != null) {
      TabixBigBedIterator tbx = new TabixBigBedIterator(this.vcfReader, chrom, start, end);
      tbx.setVcfHeader(this.vcfHeader);
      return tbx;
    } else {
      throw new RuntimeException();
    }
  }

  protected void setVcfHeader(VCFHeader vcfHeader) {
    this.vcfHeader = vcfHeader;
  }

  private ArrayList<String> getChromosomeNamesFromBigBed() {
    ChromTree chromTree;
    try {
      chromTree = new ChromTree(this.fileName, this.bigBedReader.getHeader().chromTreeOffset);
    } catch (IOException e) {
      throw new RuntimeException(e);
    }
    ArrayList<String> chromNames = new ArrayList<>();
    int i = 0;
    while (true) {
      String name = chromTree.getNameForId(i);
      if (name == null) {
        break;
      }
      chromNames.add(name);
      i++;
    }
    return chromNames;
  }
  public Set<String> getChromosomes() throws IOException {
    if (this.tabixReader != null) {
      // This is also for vcf files
      return this.tabixReader.getChromosomes();
    } else if (this.bigBedReader != null) {
      return new HashSet<>(this.getChromosomeNamesFromBigBed());
    } else if (this.vcfReader != null) {
      // This is for BCF files which don't have tabix index but have contigs in the header
      List<SAMSequenceRecord> seqs =
          this.vcfReader.getHeader().getSequenceDictionary().getSequences();
      HashSet<String> chroms = new HashSet<>();
      seqs.forEach(x -> chroms.add(x.getSequenceName()));
      return chroms;
    } else {
      throw new RuntimeException("There is no reader available to get chromosome names");
    }
  }
}
