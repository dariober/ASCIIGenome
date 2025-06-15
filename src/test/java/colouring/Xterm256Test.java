package colouring;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import exceptions.InvalidColourException;
import java.awt.Color;
import org.junit.Test;

public class Xterm256Test {

  @Test
  public void canGetContrastColour() throws InvalidColourException {
    new Xterm256();
    assertEquals("grey85", Xterm256.getContrastColour("black").toLowerCase());
  }

  @Test
  public void testXterm256() throws InvalidColourException {
    Color c = Xterm256.xterm256ToColour(30);
    System.err.println(c.getRGB()); // Not sure how to interpret this number.
  }

  @Test
  public void testColourByInteger() throws InvalidColourException {
    new Xterm256();
    assertEquals(231, Xterm256.colourNameToXterm256("grey100"));
    assertEquals(231, Xterm256.colourNameToXterm256("231"));

    // Invalid colour
    boolean pass = false;
    try {
      Xterm256.colourNameToXterm256("foo");
    } catch (InvalidColourException e) {
      pass = true;
    }
    assertTrue(pass);

    // Invalid colour as int
    pass = false;
    try {
      Xterm256.colourNameToXterm256("256");
    } catch (InvalidColourException e) {
      pass = true;
    }
    assertTrue(pass);
  }

  @Test
  public void canGetColourByApproxMatching() throws InvalidColourException {
    new Xterm256();
    assertEquals(26, Xterm256.colourNameToXterm256("DodgerBl"));
  }

  @Test
  public void canShowColours() throws InvalidColourException {
    System.out.println(Xterm256.colourShowForTerminal());
  }
}
