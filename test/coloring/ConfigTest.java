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
		new Config("resources/config/black_on_white.conf");
	}

	@Test 
	public void canGetColorForConfigKey() throws InvalidColourException, IOException, InvalidConfigException{
		new Config(null);
		int i= 0;
		while(i < 100000){ // This loop should go reasonably fast.
			assertTrue(Config.getColor(ConfigKey.background) > -1);
			i++;
		}
		System.err.println(Config.getColor(ConfigKey.read_negative_strand));
	}
	
}
