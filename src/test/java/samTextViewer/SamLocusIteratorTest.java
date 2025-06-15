package samTextViewer;

import static org.junit.Assert.assertEquals;

import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.SamLocusIterator;
import java.io.ByteArrayInputStream;
import org.junit.Test;

public class SamLocusIteratorTest {
  private SamReader createSamFileReader(final String samExample) {
    final ByteArrayInputStream inputStream = new ByteArrayInputStream(samExample.getBytes());
    return SamReaderFactory.makeDefault().open(SamInputResource.of(inputStream));
  }

  private SamLocusIterator createSamLocusIterator(final SamReader samReader) {
    final SamLocusIterator ret = new SamLocusIterator(samReader);
    ret.setEmitUncoveredLoci(false);
    return ret;
  }

  @Test
  public void testIteratorWithIntron() {
    // See https://github.com/samtools/htsjdk/issues/838
    // Prepare a sam header
    String sqHeader = "@HD\tSO:coordinate\tVN:1.0\n" + "@SQ\tSN:chr1\tAS:HG18\tLN:10000000\n";

    // Prepare one read with a 500,000 bases skipped
    String cigar = "18M500000N18M";
    String s1 =
        "read1\t0\tchr1\t1\t255\t" + cigar + "\t*\t0\t0\tACCTACGTTCAATATTACAGGCGAACATACTTACTA\t*\n";

    // Prepare sam input and samReader
    String exampleSam = sqHeader + s1 + s1;
    ByteArrayInputStream inputStream = new ByteArrayInputStream(exampleSam.getBytes());
    SamReader samReader = SamReaderFactory.makeDefault().open(SamInputResource.of(inputStream));

    // A small interval to iterate over:
    IntervalList il = new IntervalList(samReader.getFileHeader());
    il.add(new Interval("chr1", 1, 100));

    SamLocusIterator sli = new SamLocusIterator(samReader, il, true);

    // Iterate
    long t0 = System.currentTimeMillis();
    int n = 0;
    for (SamLocusIterator.LocusInfo li : sli) {
      n++;
    }
    long t1 = System.currentTimeMillis();
    System.err.println("Time to iterate " + n + " loci: " + (t1 - t0) / 1000.0 + " sec");
    sli.close();
  }

  @Test
  public void testBasicIterator() {

    final String sqHeader = "@HD\tSO:coordinate\tVN:1.0\n@SQ\tSN:chrM\tAS:HG18\tLN:100000\n";
    final String seq1 = "ACCTACGTTCAATATTACAGGCGAACATACTTACTA";
    final String qual1 = "*"; // phred 10
    final String s1 = "3851612\t16\tchrM\t165\t255\t36M\t*\t0\t0\t" + seq1 + "\t" + qual1 + "\n";
    final String exampleSam = sqHeader + s1 + s1;

    final SamReader samReader = createSamFileReader(exampleSam);
    final SamLocusIterator sli = createSamLocusIterator(samReader);

    // make sure we accumulated depth of 2 for each position
    int pos = 165;
    for (final SamLocusIterator.LocusInfo li : sli) {
      assertEquals(pos++, li.getPosition());
      assertEquals(2, li.getRecordAndOffsets().size());
    }
  }
}
