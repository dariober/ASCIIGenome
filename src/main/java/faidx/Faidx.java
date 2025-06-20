package faidx;

import com.google.common.base.CharMatcher;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.nio.Buffer;
import java.nio.ByteBuffer;
import java.nio.channels.FileChannel;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.zip.GZIPInputStream;
import java.util.zip.ZipException;

public class Faidx {

  public Faidx(File fasta) throws UnindexableFastaFileException, IOException {
    FastaFormat fastaFormat = FastaFormat.FASTA;
    if (fasta.getName().toLowerCase().endsWith(".gff")
        || fasta.getName().toLowerCase().endsWith(".gff3")
        || fasta.getName().toLowerCase().endsWith(".gtf")) {
      fastaFormat = FastaFormat.GFX;
    }
    new Faidx(fasta, fastaFormat);
  }

  /**
   * Read input fasta and write out the corresponding .fai index. See tests/faidx for examples.
   * Basically just do:
   *
   * <p>new Faidx(new File("genome.fa"));
   *
   * <p>Index will be genome.fa.fai
   */
  public Faidx(File fasta, FastaFormat fastaFormat)
      throws IOException, UnindexableFastaFileException {

    if (this.isCompressed(fasta)) {
      // System.err.println(fasta.getAbsolutePath() + " is gzip compressed. Indexing of gzip file is
      // not supported.");
      throw new UnindexableFastaFileException(
          fasta.getAbsolutePath()
              + " is gzip compressed. Compressed files must be bgzip'd and have fai and gzi"
              + " indexes.");
    }

    FileChannel fileChannel = FileChannel.open(Paths.get(fasta.getAbsolutePath()));
    int noOfBytesRead = 0;
    StringBuilder sb = new StringBuilder();

    boolean isFirstSeqLine = false;
    long currOffset = 0;
    long prevOffset = 0;
    boolean isLast = false; // True when line is expected to be the last one of sequence
    boolean startFastaSection = fastaFormat.equals(FastaFormat.FASTA) ? true : false;

    Set<String> seqNames = new HashSet<String>();
    List<FaidxRecord> records = new ArrayList<FaidxRecord>();
    FaidxRecord faidxRecord = null;

    while (noOfBytesRead != -1) {
      ByteBuffer buffer = ByteBuffer.allocate(100000);
      noOfBytesRead = fileChannel.read(buffer);
      ((Buffer) buffer)
          .flip(); // Cast to Buffer because of https://jira.mongodb.org/browse/JAVA-2559

      while (buffer.hasRemaining()) {
        char x = (char) buffer.get();

        if (startFastaSection && !CharMatcher.ascii().matches(x)) {
          throw new UnindexableFastaFileException(
              "Non ascii characters found in " + fasta.getAbsoluteFile());
        }
        currOffset++;
        sb.append(x);
        if (x == '\n') { // One full line read.
          String line = sb.toString();
          sb.setLength(0);
          if (line.trim().isEmpty()) {
            isLast = true;
            continue;
          }
          if (fastaFormat.equals(FastaFormat.GFX) && line.trim().equals("##FASTA")) {
            startFastaSection = true;
            continue;
          }
          if (fastaFormat.equals(FastaFormat.GFX) && !startFastaSection) {
            continue;
          }
          if (line.startsWith(">")) {
            isLast = false;
            if (faidxRecord != null) {
              records.add(faidxRecord);
            }
            faidxRecord = new FaidxRecord();
            faidxRecord.makeSeqNameFromRawLine(line);

            if (seqNames.contains(faidxRecord.getSeqName())) {
              throw new UnindexableFastaFileException(
                  fasta.getAbsolutePath()
                      + ": Duplicate sequence name found for "
                      + faidxRecord.getSeqName());
            } else {
              seqNames.add(faidxRecord.getSeqName());
            }
            faidxRecord.byteOffset = currOffset;
            isFirstSeqLine = true;
          } else {
            if (isLast) {
              throw new UnindexableFastaFileException(
                  fasta.getAbsolutePath()
                      + ": Different line length in "
                      + faidxRecord.getSeqName());
            }
            int seqLen = line.replaceAll("\\s", "").length();
            faidxRecord.seqLength += seqLen;
            if (isFirstSeqLine) {
              faidxRecord.lineLength = seqLen;
              faidxRecord.lineFullLength = (int) (currOffset - prevOffset);
              isFirstSeqLine = false;
            } else if (faidxRecord.lineLength != seqLen) {
              isLast = true;
            }
          }
          prevOffset = currOffset;
        } // End of processing one full line
      } // End reading chunk of bytes
    } // End of reading channel.

    records.add(faidxRecord); // Add last record

    if (fastaFormat.equals(FastaFormat.GFX) && !startFastaSection) {
      throw new UnindexableFastaFileException(
          "No fasta sequence found in " + fasta.getAbsolutePath());
    }

    // Write out index
    BufferedWriter wr =
        new BufferedWriter(new FileWriter(new File(fasta.getAbsolutePath() + ".fai")));
    for (FaidxRecord rec : records) {
      wr.write(rec.toString() + "\n");
    }
    wr.close();
  }

  private boolean isCompressed(File fasta) throws IOException {
    try {
      InputStream fileStream = new FileInputStream(fasta);
      InputStream gzipStream = new GZIPInputStream(fileStream);
      gzipStream.close();
    } catch (ZipException e) {
      return false;
    }
    return true;
  }
}
