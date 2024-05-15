package session;

import coloring.Config;
import exceptions.InvalidConfigException;
import exceptions.InvalidGenomicCoordsException;
import exceptions.SessionException;
import org.junit.Test;
import samTextViewer.GenomicCoords;

import java.io.IOException;

import static org.junit.Assert.assertEquals;

public class SessionTest {
    /** TESTS:
        * genome.file not found
        * genome.file not valid
        * genome.file empty or not present (session doesn't have a genome)
        * genome.region Invalid (malformed, extending beyond contig, contig not found in genome.file)
        * genome.region empty or not present
     */
    @Test
    public void canReadValidGenome() throws IOException, InvalidConfigException, SessionException, InvalidGenomicCoordsException {
        new Config(null);
        Session session = new Session("test_data/session.yaml", "S2");
        GenomicCoords gc = session.getGenome();
        assertEquals((long)gc.getFrom(), 1);
        assertEquals((long)gc.getTo(), 10000);
        assertEquals(gc.getSamSeqDict().getSequence("chr7").getSequenceName(), "chr7");
    }
}
