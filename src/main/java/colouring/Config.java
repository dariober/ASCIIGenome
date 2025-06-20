package colouring;

import com.google.common.base.Joiner;
import com.google.common.base.Splitter;
import exceptions.InvalidColourException;
import exceptions.InvalidConfigException;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import org.apache.commons.io.FileUtils;
import org.apache.commons.io.IOUtils;
import samTextViewer.Utils;

/** Prepare and handle configuration. */
public class Config {

  // C O N S T R U C T O R

  private static final Map<ConfigKey, String> config = new HashMap<ConfigKey, String>();

  public Config(String source) throws IOException, InvalidConfigException {

    new Xterm256();

    String RawConfigFile = Config.getConfigFileAsString(source).replaceAll("\t", " ").toLowerCase();

    // This will give one string per line
    List<String> raw =
        Splitter.on("\n").omitEmptyStrings().trimResults().splitToList(RawConfigFile);

    // List<String> config= new ArrayList<String>();
    for (String x : raw) {
      x = x.replaceAll("#.*", "").trim();
      if (!x.isEmpty()) {
        List<String> keyValuePair =
            Splitter.on(" ").omitEmptyStrings().trimResults().splitToList(x);
        if (ConfigKey.getValues().contains(keyValuePair.get(0))) {
          config.put(ConfigKey.valueOf(keyValuePair.get(0)), keyValuePair.get(1));
        } else {
          System.err.println("Unrecognized configuration key: " + keyValuePair.get(0));
          // throw new RuntimeException();
          throw new InvalidConfigException();
        }
      }
    }
    if (!config.containsKey(ConfigKey.low_mapq)) {
      // For backward compatibility: low_mapq introduced in v1.18 so previous
      // config fails will fail. This effectively is a default value if low_mapq is missing
      config.put(ConfigKey.low_mapq, "5");
    }

    // Check all fields have been populated
    for (ConfigKey key : ConfigKey.values()) {
      if (!config.containsKey(key)) {
        System.err.println("Missing configuration key: " + key);
        throw new InvalidConfigException();
      }
    }
    try {
      colourNameToInt();
    } catch (InvalidColourException e) {
      e.printStackTrace();
    }
  }

  // M E T H O D S
  public static String help() {
    List<String> help = new ArrayList<String>();
    for (ConfigKey key : ConfigKey.values()) {
      help.add(key + "\t" + Config.get(key) + "\t#\t" + key.getDescription());
    }
    List<String> table = Utils.tabulateList(help, -1, " ");
    return Joiner.on("\n").join(table);
  }

  private static String getConfigFileAsString(String source)
      throws IOException, InvalidConfigException {

    String rawConfigFile = "";

    // Is source null or empty string?
    if (source == null || source.isEmpty()) {
      // See if default config file exists
      File def = new File(System.getProperty("user.home"), ".asciigenome_config");
      if (def.isFile()) {
        rawConfigFile = FileUtils.readFileToString(def, "UTF-8");
      } else {
        // If not, read from resource
        InputStream res = Config.class.getResourceAsStream("/config/black_on_white.conf");
        rawConfigFile = IOUtils.toString(res, "UTF-8");
      }
      return rawConfigFile;
    }

    source = Utils.tildeToHomeDir(source);

    // Is source a local file? E.g. /path/to/my.conf
    if ((new File(source)).isFile()) {
      rawConfigFile = FileUtils.readFileToString(new File(source), "UTF-8");
      return rawConfigFile;
    }

    // Is source a tag matching a configuration file in resources? E.g. "black_on_white"
    try {
      InputStream res = Config.class.getResourceAsStream("/config/" + source + ".conf");
      rawConfigFile = IOUtils.toString(res, "UTF-8");
      return rawConfigFile;
    } catch (Exception e) {
      //
    }
    throw new InvalidConfigException();
  }

  /**
   * We convert the colour names to the corresponding int. This is because looking up by int is much
   * faster than by name.
   *
   * @throws InvalidColourException
   */
  private static void colourNameToInt() throws InvalidColourException {
    for (ConfigKey key : config.keySet()) {
      if (ConfigKey.colourKeys().contains(key)) {
        int colourInt = Xterm256.colourNameToXterm256(config.get(key));
        config.put(key, Integer.toString(colourInt));
      }
    }
  }

  /** Get xterm256 colour corresponding to this configuration key */
  public static int get256Colour(ConfigKey key) throws InvalidColourException {
    new Xterm256();
    return Xterm256.colourNameToXterm256(config.get(key));
  }

  /** Get value associated to this configuration key */
  public static String get(ConfigKey key) {
    return config.get(key);
  }

  public static void set(ConfigKey key, String value) throws InvalidColourException {
    if (ConfigKey.booleanKeys().contains(key)) {
      boolean bool = Utils.asBoolean(value);
      config.put(key, Boolean.toString(bool));
    } else if (ConfigKey.integerKeys().contains(key)) {
      Integer.valueOf(value);
      config.put(key, value);
    } else {
      if (key.equals(ConfigKey.shade_structural_variant) && value.equals("false")) {
        // This is very much a hack: If colour is false we reset to background
        value = config.get(ConfigKey.background);
      }
      config.put(key, value);
      colourNameToInt();
    }
  }

  //	/**Return the key-value map of configuration parameters
  //	 * */
  //	public static Map<ConfigKey, String> getConfigMap(){
  //		return config;
  //	}
}
