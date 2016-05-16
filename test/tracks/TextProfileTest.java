package tracks;

import static org.junit.Assert.*;

import java.util.ArrayList;
import java.util.List;

import org.junit.Test;

public class TextProfileTest {

	@Test
	public void testNaN(){
		System.out.println((int)Double.NaN);
	}
	
	@Test
	public void canSetMinMax(){
		List<Double> yValues= new ArrayList<Double>();
		yValues.add((double)1);
		yValues.add((double)2);
		yValues.add((double)3);
		yValues.add((double)4);
		
		TextProfile tp= new TextProfile(yValues, 10, Double.NaN, Double.NaN);
		assertEquals(1, tp.getYMinLimit(), 0.001);
		assertEquals(4, tp.getYMaxLimit(), 0.001);
		
		tp= new TextProfile(yValues, 10, -1.0, Double.NaN);
		assertEquals(-1, tp.getYMinLimit(), 0.001);
		assertEquals(4, tp.getYMaxLimit(), 0.001);
		
		tp= new TextProfile(yValues, 10, Double.NaN, 10.0);
		assertEquals(1, tp.getYMinLimit(), 0.001);
		assertEquals(10, tp.getYMaxLimit(), 0.001);
	}
	
	@Test
	public void canGraphProfile(){
		List<Double> yValues= new ArrayList<Double>();
		yValues.add((double)1);
		yValues.add((double)2);
		yValues.add((double)3);
		yValues.add((double)4);
		yValues.add((double)5);
		yValues.add((double)6);
		yValues.add((double)7);
		yValues.add((double)8);
		yValues.add((double)9);
		yValues.add((double)10);
		
		TextProfile tp= new TextProfile(yValues, 5, 0.0, Double.NaN);
		for(List<String> x : tp.getProfile()){
			System.out.println(x);
		}
		
		System.out.println(Double.NaN < 0);
		
	}
	
	@Test
	public void test() {
		
		List<Double> yValues= new ArrayList<Double>();
		yValues.add((double)0);
		yValues.add((double)0);
		yValues.add((double)0);
		yValues.add((double)0);
		
		TextProfile tp= new TextProfile(yValues, 10, Double.NaN, Double.NaN);
		//System.out.println(tp.getProfile().size());
		//System.out.println(tp.getProfile());
	}

}
