package coloring;

import static org.junit.Assert.assertTrue;

import java.io.IOException;
import java.sql.SQLException;

import org.junit.Test;

import exceptions.InvalidColourException;
import exceptions.InvalidConfigException;
import exceptions.InvalidGenomicCoordsException;
import exceptions.InvalidRecordException;
import samTextViewer.GenomicCoords;
import tracks.TrackReads;

public class ConfigTest {

	@Test
	public void canPrintHelp() throws IOException, InvalidConfigException{
		new Config(null);
		assertTrue(Config.help().length() > 100);
	}
	
	@Test
	public void canSetConfig() throws ClassNotFoundException, IOException, InvalidGenomicCoordsException, InvalidRecordException, SQLException, InvalidColourException, InvalidConfigException{
		
		Config conf= new Config(null);
		conf.set(ConfigKey.seq_a, "grey");
		String bam= "test_data/adjacent.bam";
		GenomicCoords gc= new GenomicCoords("chr7:1-50", 80, null, null);
		TrackReads tr= new TrackReads(bam, gc);
		assertTrue(tr.printToScreen().contains("8m"));
		conf= new Config(null);
	}
	
	@Test
	public void canReadConfigFromResource() throws IOException, InvalidConfigException {
		new Config(null);
	}

	@Test
	public void canReadConfigFromFile() throws IOException, InvalidConfigException {
		new Config("resources/config/black_on_white.conf");
	}

	@Test
	public void canReadConfigFromTag() throws IOException, InvalidConfigException {
		// If you add new themes add them here to test
		new Config("black_on_white");
		new Config("white_on_black");
		new Config("metal");
	}

	@Test
	public void failsOnInvalidSource() throws IOException, InvalidConfigException {
		boolean pass= false;
		try{
			new Config("foo");
		} catch(InvalidConfigException e){
			pass= true;
		}
		assertTrue(pass);
	}
	
	@Test 
	public void canGetColorForConfigKey() throws InvalidColourException, IOException, InvalidConfigException{
		new Config(null);
		int i= 0;
		while(i < 100000){ // This loop should go reasonably fast.
			assertTrue(Config.get256Color(ConfigKey.background) > -1);
			i++;
		}
	}
	
}
