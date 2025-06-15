package tracks;

import static org.junit.Assert.*;

import colouring.Config;
import colouring.Xterm256;
import exceptions.InvalidColourException;
import exceptions.InvalidConfigException;
import java.io.IOException;
import org.junit.Before;
import org.junit.Test;

public class TrackHeaderTest {

  @Before
  public void config() throws IOException, InvalidConfigException {
    new Config(null);
    new Xterm256();
  }

  @Test
  public void canInterpolateText() throws InvalidColourException {
    TrackHeader th = new TrackHeader("HEADER");
    assertTrue(th.format(true).contains("HEADER"));
    th.setHeaderText("foo {-} bar {-}");
    assertTrue(th.format(true).contains("foo HEADER bar HEADER"));
  }

  @Test
  public void canFormat() throws InvalidColourException {
    TrackHeader th = new TrackHeader("HEADER");
    th.setColour("yellow");
    assertTrue(th.format(false).endsWith("\033[0m"));
    assertTrue(th.format(false).contains("5;11m"));
    assertTrue(th.format(false).contains("HEADER"));
  }

  @Test
  public void canAlignRight() throws InvalidColourException {
    TrackHeader th = new TrackHeader("HEAD");
    th.setHeaderAlignmentPct(1);

    th.setTerminalWidth(4);
    assertEquals(th.format(true), "HEAD\n");

    th.setTerminalWidth(5);
    assertEquals(th.format(true), " HEAD\n");

    th.setTerminalWidth(6);
    assertEquals(th.format(true), "  HEAD\n");
  }

  @Test
  public void canAlignCenter() throws InvalidColourException {
    TrackHeader th = new TrackHeader("HEAD");
    th.setHeaderAlignmentPct(0.5);

    th.setTerminalWidth(4);
    assertEquals(th.format(true), "HEAD\n");

    th.setTerminalWidth(5);
    assertEquals(th.format(true), "HEAD\n");

    th.setTerminalWidth(6);
    assertEquals(th.format(true), " HEAD\n");

    th.setTerminalWidth(7);
    assertEquals(th.format(true), " HEAD\n");

    th.setTerminalWidth(8);
    assertEquals(th.format(true), "  HEAD\n");
  }
}
