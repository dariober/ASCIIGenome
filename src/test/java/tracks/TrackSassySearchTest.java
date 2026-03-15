package tracks;

import static org.junit.Assert.assertTrue;

import colouring.Config;
import colouring.ConfigKey;
import exceptions.InvalidColourException;
import exceptions.InvalidCommandLineException;
import exceptions.InvalidConfigException;
import exceptions.InvalidGenomicCoordsException;
import exceptions.InvalidRecordException;
import java.io.IOException;
import java.sql.SQLException;
import java.util.Arrays;
import java.util.List;
import org.apache.commons.lang3.StringUtils;
import org.junit.Test;
import samTextViewer.GenomicCoords;

public class TrackSassySearchTest {

  @Test
  public void canInitializeTrack()
      throws InvalidGenomicCoordsException,
      IOException, InvalidColourException, InvalidConfigException {

    new Config(null);
    Config.set(ConfigKey.sassy, "test_data/sassy/sassy-x86_64-unknown-linux-gnu");
    GenomicCoords gc = new GenomicCoords("chr7:11001-11200", 200, null, "test_data/chr7.fa");
    List<String> sassyCliOpts = Arrays.stream("-p gggaggc -k 2".split(" ")).toList();
    TrackSassySearch ss = new TrackSassySearch(sassyCliOpts, gc);
    ss.setNoFormat(true);
    String out = ss.printToScreen();
    assertTrue(out.startsWith(".......   ..G...C "));
  }

  @Test
  public void testManyFeatures()
      throws InvalidGenomicCoordsException,
      IOException, InvalidColourException, InvalidConfigException {

    new Config(null);
    Config.set(ConfigKey.sassy, "test_data/sassy/sassy-x86_64-unknown-linux-gnu");
    GenomicCoords gc = new GenomicCoords("chr7:1-14557334", 200, null, "test_data/chr7.fa");
    List<String> sassyCliOpts = Arrays.stream("-p GCGGCCGC -k 1".split(" ")).toList();
    TrackSassySearch ss = new TrackSassySearch(sassyCliOpts, gc);
    ss.setNoFormat(true);
    String out = ss.printToScreen();
    int count = StringUtils.countMatches(out, ">");
    assertTrue(count > 50);
    count = StringUtils.countMatches(out, "<");
    assertTrue(count > 50);
  }
}
