package tracks;

import static org.junit.Assert.assertTrue;

import aligner.SassyTest;
import colouring.Config;
import colouring.ConfigKey;
import exceptions.InvalidColourException;
import exceptions.InvalidConfigException;
import exceptions.InvalidGenomicCoordsException;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import org.apache.commons.lang3.StringUtils;
import org.junit.Test;
import samTextViewer.GenomicCoords;

public class TrackSassySearchTest {

  @Test
  public void canInitializeTrack()
      throws InvalidGenomicCoordsException,
          IOException,
          InvalidColourException,
          InvalidConfigException {

    new Config(null);
    Config.set(ConfigKey.sassy, SassyTest.sassyExec().toString());
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
          IOException,
          InvalidColourException,
          InvalidConfigException {

    new Config(null);
    Config.set(ConfigKey.sassy, SassyTest.sassyExec().toString());
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
