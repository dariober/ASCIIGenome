package coloring;

import static org.junit.Assert.assertTrue;

import java.io.IOException;

import org.junit.Test;

import exceptions.InvalidColourException;
import exceptions.InvalidConfigException;

public class ConfigTest {

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
			assertTrue(Config.getColor(ConfigKey.background) > -1);
			i++;
		}
	}
	
}
