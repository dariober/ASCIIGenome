package tracks;

import static org.junit.Assert.*;

import colouring.Config;
import exceptions.InvalidColourException;
import exceptions.InvalidConfigException;
import java.io.IOException;
import org.junit.Before;
import org.junit.Test;

public class FeatureCharTest {

  @Before
  public void config() throws IOException, InvalidConfigException {
    new Config(null);
  }

  @Test
  public void testFormatting() throws InvalidColourException {
    FeatureChar fc = new FeatureChar();
    // Formatting VCF
    fc.addFormatVCF('A');
    assertEquals("A", fc.format(true));
    assertTrue(fc.format(false).startsWith("\033["));

    fc.setBgColour("255");
    assertTrue(fc.format(false).contains("255"));
    fc.setFgColour("254");
    assertTrue(fc.format(false).contains("254"));

    // Formatting general text
    fc = new FeatureChar();
    fc.setText('X');
    fc.setFgColour("blue");
  }
}
