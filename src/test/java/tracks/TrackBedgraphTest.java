package tracks;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.IOException;
import java.security.NoSuchAlgorithmException;
import java.sql.SQLException;
import java.util.List;

import org.junit.Before;
import org.junit.Test;

import coloring.Config;
import exceptions.InvalidColourException;
import exceptions.InvalidCommandLineException;
import exceptions.InvalidConfigException;
import exceptions.InvalidGenomicCoordsException;
import exceptions.InvalidRecordException;
import samTextViewer.GenomicCoords;

public class TrackBedgraphTest {

    @Before
    public void prepareConfig() throws IOException, InvalidConfigException{
        new Config(null);
    }

    @Test
    public void canCloseReaders() throws ClassNotFoundException, IOException, InvalidRecordException, InvalidGenomicCoordsException, SQLException{
        GenomicCoords gc= new GenomicCoords("chr7:5540000-5570000", 80, null, null);
        TrackBedgraph tb= new TrackBedgraph("test_data/test.bedGraph", gc);
        tb.close();
    }
    
    @Test
    public void canPrintChromosomeNames() throws InvalidGenomicCoordsException, IOException, ClassNotFoundException, InvalidRecordException, SQLException{

        GenomicCoords gc= new GenomicCoords("chr7:5540000-5570000", 80, null, null);
        
        TrackBedgraph tb= new TrackBedgraph("test_data/test.bedGraph", gc);
        assertTrue(tb.getChromosomeNames().size() > 0);
        assertEquals("chr1", tb.getChromosomeNames().get(0));
    }
    
    @Test
    public void canGetDataColumnIndexForBedGraph() throws IOException, NoSuchAlgorithmException, InvalidGenomicCoordsException, InvalidRecordException, ClassNotFoundException, SQLException{
        String url= "test_data/test.bedGraph";
        GenomicCoords gc= new GenomicCoords("chr1:1-30", 80, null, null);
        TrackBedgraph tb= new TrackBedgraph(url, gc);
        tb.setScoreColIdx(5);
        assertEquals(0, tb.getScreenScores().get(0), 0.0001);
        
        tb= new TrackBedgraph(url, gc);
        assertEquals(1, tb.getScreenScores().get(0), 0.0001);

        tb.setScoreColIdx(20);
        assertEquals(Float.NaN, tb.getScreenScores().get(0), 0.0001);
    }

    @Test
    public void invalidBedgraphRecordAsNaN() throws IOException, NoSuchAlgorithmException, InvalidGenomicCoordsException, InvalidRecordException, ClassNotFoundException, SQLException{
        
        GenomicCoords gc= new GenomicCoords("chr1:1-10", 10, null, null);
        TrackBedgraph tb= new TrackBedgraph("test_data/invalid-1.bedgraph", gc);
        assertEquals(tb.getScreenScores().get(0), Float.NaN, 0.0001);
        assertEquals(tb.getScreenScores().get(1), 10, 0.0001);
        
    }
        
    @Test
    public void canParseNonBGZFFile() throws IOException, InvalidGenomicCoordsException, InvalidRecordException, ClassNotFoundException, SQLException{
        
        String url= "test_data/test2.bedGraph";
        GenomicCoords gc= new GenomicCoords("chr1:1-30", 80, null, null);
        new TrackBedgraph(url, gc);
                
    }
    
    @Test
    public void testYLimits() throws InvalidGenomicCoordsException, IOException, InvalidColourException, InvalidRecordException, ClassNotFoundException, SQLException{

        String url= "test_data/test.bedGraph.gz";
        GenomicCoords gc= new GenomicCoords("chr1:1-30", 80, null, null);
                
        TrackBedgraph tb= new TrackBedgraph(url, gc);
        tb.setYLimitMax((float)10.0);
        tb.setYLimitMin((float)-10.0);
        tb.setyMaxLines(10);
        String prof= tb.printToScreen();
        assertTrue(prof.contains(","));
        assertTrue(prof.contains(":"));
    }
    
    @Test
    public void testCloseToBorder() throws InvalidGenomicCoordsException, InvalidColourException, IOException, InvalidRecordException, ClassNotFoundException, SQLException{
        String url= "test_data/test.bedGraph.gz";
        int yMaxLines= 10;
        GenomicCoords gc= new GenomicCoords("chr1:1-800", 80, null, null);
        TrackBedgraph tb= new TrackBedgraph(url, gc);
        tb.setYLimitMax(Float.NaN);
        tb.setYLimitMin(Float.NaN);
        tb.setyMaxLines(yMaxLines);
        tb.printToScreen();
    } 
    
    @Test 
    public void canPrintLines() throws InvalidGenomicCoordsException, IOException, InvalidRecordException, ClassNotFoundException, SQLException, InvalidColourException, InvalidConfigException, InvalidCommandLineException{
        new Config(null);
            
        String url= "test_data/test.bedGraph.gz";
        GenomicCoords gc= new GenomicCoords("chr1:1-22", 80, null, null);
        TrackBedgraph tb= new TrackBedgraph(url, gc);
        tb.setPrintMode(PrintRawLine.CLIP);
        tb.setPrintRawLineCount(-1);
        tb.setNoFormat(true);
        assertTrue(tb.printLines().startsWith("chr1 "));
    }
    
    @Test
    public void canApplyAwk() throws IOException, InvalidGenomicCoordsException, ClassNotFoundException, InvalidRecordException, SQLException, InvalidColourException{
        
        String fileName= "test_data/test.bedGraph.gz";
        GenomicCoords gc= new GenomicCoords("chr1:1-22", 80, null, null);
        TrackBedgraph tb= new TrackBedgraph(fileName, gc);
        tb.setNoFormat(true);

        String before = tb.printToScreen();
        
        tb.setAwk("'$4 > 0");
        String after = tb.printToScreen();
        
        List<IntervalFeature> subset = tb.getFeaturesInInterval("chr1", 1, 500000000);
        for(IntervalFeature s : subset) {
            assertTrue(s.getScore() > 0);
        }

        assertTrue(!before.equals(after));
    }
    
    @Test
    public void canPrintTitle() throws ClassNotFoundException, IOException, InvalidGenomicCoordsException, InvalidRecordException, SQLException, InvalidColourException {
        String url= "test_data/test.bedGraph.gz";
        GenomicCoords gc= new GenomicCoords("chr1:1-22", 80, null, null);
        TrackBedgraph tb= new TrackBedgraph(url, gc);
        tb.setAwk("'$1 > 0");
        tb.setYLimitMax(-99);
        tb.setNoFormat(true);
        
        String title = tb.getTitle();

        assertTrue(title.contains("ylim[auto -99.0]"));
        assertTrue(title.contains("awk"));
    }
    
    @Test 
    public void canPrintBedGraph() throws InvalidGenomicCoordsException, IOException, InvalidRecordException, ClassNotFoundException, SQLException, InvalidColourException, InvalidConfigException{
        
        new Config(null);
        
        String url= "test_data/test.bedGraph.gz";
        int yMaxLines= 5;
        GenomicCoords gc= new GenomicCoords("chr1:1-22", 80, null, null);
        TrackBedgraph tb= new TrackBedgraph(url, gc);
        tb.setYLimitMax(Float.NaN);
        tb.setYLimitMin(Float.NaN);
        tb.setyMaxLines(yMaxLines);
        String prof= tb.printToScreen();
        System.out.println(prof);
        
        tb= new TrackBedgraph("test_data/positive.bedGraph.gz", gc);
        tb.setYLimitMax(Float.NaN);
        tb.setYLimitMin(Float.NaN);
        tb.setyMaxLines(5);
        prof= tb.printToScreen();
        System.out.println(prof);
        
        tb= new TrackBedgraph("test_data/negative.bedGraph.gz", gc);
        tb.setYLimitMax(Float.NaN);
        tb.setYLimitMin(Float.NaN);
        tb.setyMaxLines(5);
        
        System.out.println(prof);
        
        gc= new GenomicCoords("chr1:1-52", 80, null, null);
        tb= new TrackBedgraph("test_data/posNeg.bedGraph.gz", gc);
        tb.setYLimitMax(Float.NaN);
        tb.setYLimitMin(Float.NaN);
        tb.setyMaxLines(14);
        System.out.println(tb.printToScreen());
        
    }
}
