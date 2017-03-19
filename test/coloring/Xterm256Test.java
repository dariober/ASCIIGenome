package coloring;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.awt.Color;

import org.junit.Test;

import exceptions.InvalidColourException;

public class Xterm256Test {

	@Test
	public void testXterm256() throws InvalidColourException{
		
		Color c= Xterm256.xterm256ToColor(30);
		System.err.println(c.getRGB());
	}
	
	@Test
	public void testColorByInteger() throws InvalidColourException{
		
		assertEquals(231, Xterm256.colorNameToXterm256("grey100"));
		assertEquals(231, Xterm256.colorNameToXterm256("231"));
		
		// Invalid colour
		boolean pass= false;
		try{
			Xterm256.colorNameToXterm256("foo");
		} catch(InvalidColourException e){
			pass= true;
		}
		assertTrue(pass);
		
		// Invalid colour as int
		pass= false;
		try{
			Xterm256.colorNameToXterm256("256");
		} catch(InvalidColourException e){
			pass= true;
		}
		assertTrue(pass);
	}

	@Test
	public void canGetColorByApproxMatching() throws InvalidColourException{
		assertEquals(26, Xterm256.colorNameToXterm256("DodgerBl"));
	}
	
	@Test
	public void canShowColors() throws InvalidColourException{
		System.out.println(Xterm256.colorShowForTerminal());
	}
}
