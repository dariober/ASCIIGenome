package samTextViewer;

import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.List;

import org.junit.Test;

import com.itextpdf.text.DocumentException;

import coloring.Config;
import exceptions.BamIndexNotFoundException;
import exceptions.InvalidColourException;
import exceptions.InvalidCommandLineException;
import exceptions.InvalidConfigException;
import exceptions.InvalidGenomicCoordsException;
import exceptions.InvalidRecordException;
import tracks.TrackSet;

public class TrackProcessorTest {

    @Test
    public void canProcessTracks() throws InvalidGenomicCoordsException, IOException, ClassNotFoundException, BamIndexNotFoundException, InvalidRecordException, SQLException, InvalidConfigException, InvalidCommandLineException, DocumentException, InvalidColourException{

        new Config(null);
        
        GenomicCoords gc= new GenomicCoords("chr7:1-100", 80, null, null);
        List<String> genome= new ArrayList<String>();
        genome.add("test_data/ear045.oxBS.actb.bam");
        gc.setGenome(genome, true);

        GenomicCoordsHistory gch= new GenomicCoordsHistory();
        gch.add(gc);
        
        TrackSet trackSet= new TrackSet(new ArrayList<String>(), gc);
        trackSet.addTrackFromSource("test_data/hg19_genes.gtf.gz", gc, null);
        trackSet.addTrackFromSource("test_data/ear045.oxBS.actb.bam", gc, null);
        trackSet.addTrackFromSource("test_data/ear045.oxBS.actb.tdf", gc, null);
        new TrackProcessor(trackSet, gch);
        // tp.iterateTracks();
    }
    
    @Test
    public void canProcessTracksNewline() throws InvalidGenomicCoordsException, IOException, ClassNotFoundException, BamIndexNotFoundException, InvalidRecordException, SQLException, InvalidConfigException, InvalidCommandLineException, DocumentException, InvalidColourException{

            new Config(null);
            
            GenomicCoords gc= new GenomicCoords("chr7:1-100", 80, null, null);
            List<String> genome= new ArrayList<String>();
            genome.add("test_data/ear045.oxBS.actb.bam");
            gc.setGenome(genome, true);

            GenomicCoordsHistory gch= new GenomicCoordsHistory();
            gch.add(gc);
            
            TrackSet trackSet= new TrackSet(new ArrayList<String>(), gc);
            trackSet.addTrackFromSource("test_data/hg19_genes.gtf.gz", gc, null);
            trackSet.editNamesForRegex(Utils.tokenize("editNames test_ foo\ntest_", " "));
            
            TrackProcessor tp = new TrackProcessor(trackSet, gch);
            tp.setNoFormat(true);
            tp.iterateTracks();
        }
}
