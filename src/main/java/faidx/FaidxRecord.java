package faidx;

/**
 * This class represents an entry of the fasta index (i.e. a line of the .fai file)
 *
 * <p>Example of .fai file: regular 122 9 32 33 sequenc_in_one_line 64 156 64 65 with 122 241 32 33
 *
 * <p>Columns are:
 *
 * <p>1. Sequence name 2. Sequence length 3. Byte pos of sequence start 4. Line length (only chars)
 * 5. Line length (including line terminators)
 */
class FaidxRecord {
  private String seqName = null;
  protected int seqLength = 0;
  protected long byteOffset = 0; // Byte position where the sequence starts.
  protected int lineLength = 0; // Only nucleotides
  protected int lineFullLength = 0; // Including line terminators.

  protected void makeSeqNameFromRawLine(String line) throws UnindexableFastaFileException {
    if (!line.startsWith(">")) {
      throw new UnindexableFastaFileException("Invalid name: Does not start with '>'");
    }
    String name = line.substring(1).trim().replaceAll("\\s.*", "");
    this.seqName = name;
  }

  protected String getSeqName() {
    return seqName;
  }

  @Override
  public String toString() {
    return seqName
        + "\t"
        + seqLength
        + "\t"
        + byteOffset
        + "\t"
        + lineLength
        + "\t"
        + lineFullLength;
  }
}
