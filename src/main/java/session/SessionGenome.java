package session;

import exceptions.InvalidGenomicCoordsException;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import java.io.File;
import java.io.IOException;
import samTextViewer.GenomicCoords;
import samTextViewer.Utils;

public class SessionGenome {

  public String samSeqDictSource;
  public String chrom;
  public Integer from;
  public Integer to;
  public String fastaFile;

  public SessionGenome() {}

  /* Handle GenomicCoords object to make suitable for serialization **/
  public SessionGenome(GenomicCoords gc) {
    this.chrom = gc.getChrom();
    this.from = gc.getFrom();
    this.to = gc.getTo();
    this.fastaFile =
        gc.getOriginalFastaFile() == null ? null : new File(gc.getOriginalFastaFile()).getAbsolutePath();
    this.samSeqDictSource = gc.getSamSeqDictSource();
  }

  public GenomicCoords toGenomicCoords() throws InvalidGenomicCoordsException, IOException {
    String region = this.chrom + ':' + this.from + '-' + this.to;
    String ff = null;
    if (this.fastaFile != null && new File(this.fastaFile).exists()) {
      ff = this.fastaFile;
    }
    SAMSequenceDictionary samSeqDict = null;
    if (this.samSeqDictSource != null && new File(this.samSeqDictSource).exists()) {
      SamReaderFactory srf = SamReaderFactory.make();
      try (SamReader samReader = srf.open(new File(this.samSeqDictSource))) {
        samSeqDict = samReader.getFileHeader().getSequenceDictionary();
      }
    }
    return new GenomicCoords(region, Utils.getTerminalWidth(), samSeqDict, ff);
  }

  @Override
  public String toString() {
    return "SessionGenome{"
        + "chrom='"
        + chrom
        + '\''
        + ", from="
        + from
        + ", to="
        + to
        + ", fastaFile='"
        + fastaFile
        + '\''
        + '}';
  }
}
