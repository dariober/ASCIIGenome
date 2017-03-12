package coloring;

import java.io.IOException;
import java.io.InputStream;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.io.IOUtils;

import com.google.common.base.Splitter;

import exceptions.InvalidColourException;

/** Prepare and handle configuration. 
 * */
public class Config {
	
	// C O N S T R U C T O R 
	
	private static final Map<ConfigKey, String> config= new HashMap<ConfigKey, String>();

	public Config(String source) {
		if(source == null || source.isEmpty()){
			
			InputStream res= Config.class.getResourceAsStream("/config/black_on_white.conf");
			String xraw;
			try {
				xraw = IOUtils.toString(res, "UTF-8").replaceAll("\t", " ").toLowerCase();
			} catch (IOException e) {
				e.printStackTrace();
				throw new RuntimeException();
			}
			
			// This will give one string per line
			List<String> raw= Splitter.on("\n").omitEmptyStrings().trimResults().splitToList(xraw);

			//List<String> config= new ArrayList<String>();
			for(String x : raw){
				x= x.replaceAll("#.*", "").trim();
				if( ! x.isEmpty()){
					List<String> keyValuePair= Splitter.on(" ").omitEmptyStrings().trimResults().splitToList(x);
					if(ConfigKey.getValues().contains(keyValuePair.get(0))){
						config.put(ConfigKey.valueOf(keyValuePair.get(0)), keyValuePair.get(1));
					} else {
						System.err.println("Unrecognized configuration key: " + keyValuePair.get(0));
						throw new RuntimeException();
						// throw new InvalidConfigException();
					}
				}
			}
		}
		// Check all fields have been populated
		for(ConfigKey key : ConfigKey.values()){
			if( ! config.containsKey(key)){
				System.err.println("Missing configuration key: " + key);
				throw new RuntimeException();
				// throw new InvalidConfigException();
			}
		}
	}

	// M E T H O D S
	
	/** Get xterm256 color corresponding to this configuration key
	 * */
	public static int getColor(ConfigKey key) throws InvalidColourException{
		return Xterm256.colorNameToXterm256(config.get(key));
	}

	/** Get value associated to this configuration key 
	 * */
	public static String get(ConfigKey key) {
		return config.get(key);
	}		

	
}
