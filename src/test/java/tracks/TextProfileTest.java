package tracks;

import static org.junit.Assert.*;

import java.util.ArrayList;
import java.util.List;
import org.junit.Test;

public class TextProfileTest {

  @Test
  public void testNaN() {
    System.out.println((int) Double.NaN);
  }

  @Test
  public void canSetMinMax() {
    List<Float> yValues = new ArrayList<Float>();
    yValues.add((float) 1);
    yValues.add((float) 2);
    yValues.add((float) 3);
    yValues.add((float) 4);

    TextProfile tp = new TextProfile(yValues, 10, Float.NaN, Float.NaN);
    assertEquals(1, tp.getYMinLimit(), 0.001);
    assertEquals(4, tp.getYMaxLimit(), 0.001);

    tp = new TextProfile(yValues, 10, (float) -1.0, Float.NaN);
    assertEquals(-1, tp.getYMinLimit(), 0.001);
    assertEquals(4, tp.getYMaxLimit(), 0.001);

    tp = new TextProfile(yValues, 10, Float.NaN, (float) 10.0);
    assertEquals(1, tp.getYMinLimit(), 0.001);
    assertEquals(10, tp.getYMaxLimit(), 0.001);
  }

  @Test
  public void canGraphProfile() {
    List<Float> yValues = new ArrayList<Float>();
    yValues.add((float) 1);
    yValues.add((float) 2);
    yValues.add((float) 3);
    yValues.add((float) 4);
    yValues.add((float) 5);
    yValues.add((float) 6);
    yValues.add((float) 7);
    yValues.add((float) 8);
    yValues.add((float) 9);
    yValues.add((float) 10);

    TextProfile tp = new TextProfile(yValues, 5, (float) 0.0, Float.NaN);
    for (List<String> x : tp.getProfile()) {
      System.out.println(x);
    }

    System.out.println(Double.NaN < 0);
  }

  @Test
  public void test() {

    List<Float> yValues = new ArrayList<Float>();
    yValues.add((float) 0);
    yValues.add((float) 0);
    yValues.add((float) 0);
    yValues.add((float) 0);

    TextProfile tp = new TextProfile(yValues, 10, Float.NaN, Float.NaN);
    // System.out.println(tp.getProfile().size());
    // System.out.println(tp.getProfile());
  }
}
