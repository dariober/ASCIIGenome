package tracks;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import org.junit.Before;
import org.junit.Test;

import coloring.Config;
import coloring.Xterm256;
import exceptions.BamIndexNotFoundException;
import exceptions.InvalidConfigException;
import exceptions.InvalidGenomicCoordsException;
import exceptions.InvalidRecordException;
import samTextViewer.GenomicCoords;

/**Test failing on Travis for unknown reasons
 * */
public class TravisExcludedTest {

    @Before
    public void config() throws IOException, InvalidConfigException{
        new Config(null);
        new Xterm256();
    }
    
    @Test
    public void canGetListOfOpenedFiles() throws ClassNotFoundException, IOException, BamIndexNotFoundException, InvalidGenomicCoordsException, InvalidRecordException, SQLException{

        GenomicCoords gc= new GenomicCoords("chr7:5565052-5571960", 80, null, null);
        
        TrackSet trackSet= new TrackSet(new ArrayList<String>(), gc);
        trackSet.addTrackFromSource("test_data/refSeq.hg19.bed.gz", gc, null);
        trackSet.addTrackFromSource("test_data/hg19_genes.gtf.gz", gc, null);
        trackSet.addTrackFromSource("test_data/refSeq.hg19.bed.gz", gc, null);
        trackSet.addTrackFromSource("ftp://ftp.ensembl.org/pub/release-86/gff3/homo_sapiens/Homo_sapiens.GRCh38.86.chromosome.18.gff3.gz", gc, null);
        
        assertTrue(new File(trackSet.getOpenedFiles().iterator().next()).isAbsolute());
        assertEquals(3, trackSet.getOpenedFiles().size());
        assertTrue(new File(trackSet.getOpenedFiles().iterator().next()).isAbsolute());
    }
    
    @Test
    public void canConstructTrackSetFromURL() throws InvalidGenomicCoordsException, IOException, ClassNotFoundException, BamIndexNotFoundException, InvalidRecordException, SQLException{

        GenomicCoords gc= new GenomicCoords("chr7:5565052-5571960", 80, null, null);
        
        TrackSet trackSet= new TrackSet(new ArrayList<String>(), gc);
        trackSet.addTrackFromSource("http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeCaltechRnaSeq/wgEncodeCaltechRnaSeqGm12878R2x75Il400GeneGencV3cRep2V3.gtf.gz", gc, null);
        trackSet.addTrackFromSource("ftp://ftp.ensembl.org/pub/release-86/gff3/homo_sapiens/Homo_sapiens.GRCh38.86.chromosome.18.gff3.gz", gc, null);
    }
  
}
