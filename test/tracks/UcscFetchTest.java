package tracks;

import static org.junit.Assert.*;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.Arrays;
import java.util.List;

import org.junit.Test;

import com.google.common.base.Splitter;

import exceptions.InvalidCommandLineException;
import exceptions.InvalidGenomicCoordsException;
import exceptions.InvalidRecordException;
import exceptions.UnableToExecuteUtilException;
import samTextViewer.GenomicCoords;

public class UcscFetchTest {

	@Test
	public void canGetGTFfromUCSC() throws InvalidGenomicCoordsException, IOException, InterruptedException, UnableToExecuteUtilException, InvalidCommandLineException, InvalidRecordException {
		
		System.out.print("This might take a while...");
		UcscFetch ucscTrack= new UcscFetch("dm6:refGene");
		System.out.println("Done");
		
		File gtf= ucscTrack.genePredToGtf();
		assertTrue(gtf.exists());
		assertTrue(gtf.length() > 1e6); // File exists and it's decent legnth. 
		assertTrue((new File(gtf.getAbsolutePath() + ".tbi")).length() > 10000); // Check indexing completed ok.
	}

	@Test
	public void canFailGracefullyOnNonExistantTable() throws InvalidGenomicCoordsException, IOException, InterruptedException, UnableToExecuteUtilException, InvalidCommandLineException, InvalidRecordException {
		
		boolean passed= false;
		try{
			UcscFetch ucscTrack= new UcscFetch("hg19:foobar");
		} catch(InvalidCommandLineException e){
			passed= true;
		}
		assertTrue(passed);
	}
	
	@Test
	public void canFailGracefullyOnWrongSyntax() throws InvalidGenomicCoordsException, IOException, InterruptedException, UnableToExecuteUtilException, InvalidCommandLineException, InvalidRecordException {
		
		boolean passed= false;
		try{
			UcscFetch ucscTrack= new UcscFetch("hg19-foobar");
		} catch(InvalidCommandLineException e){
			passed= true;
		}
		assertTrue(passed);
	}
	
	@Test
	public void canFailGracefullyOnNonGenePredTable() throws InvalidGenomicCoordsException, IOException, InterruptedException, UnableToExecuteUtilException, InvalidCommandLineException, InvalidRecordException {
		
		boolean passed= false;
		try{
			UcscFetch ucscTrack= new UcscFetch("hg19:chromInfo");
		} catch(InvalidCommandLineException e){
			passed= true;
		}
		assertTrue(passed);
	}
	
	@Test
	public void canFailGracefullyOnNonExistingDB() throws InvalidGenomicCoordsException, IOException, InterruptedException, UnableToExecuteUtilException, InvalidCommandLineException, InvalidRecordException {
		
		boolean passed= false;
		try{
			UcscFetch ucscTrack= new UcscFetch("hg999:knownGene");
		} catch(InvalidCommandLineException e){
			passed= true;
		}
		assertTrue(passed);
	}
}
