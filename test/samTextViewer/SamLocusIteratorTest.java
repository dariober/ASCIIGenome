package samTextViewer;

import static org.junit.Assert.*;

import java.io.ByteArrayInputStream;

import org.junit.Test;

import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

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
            assertEquals(2, li.getRecordAndPositions().size());
        }

    }
}