package tracks;

import coloring.Xterm256;
import exceptions.InvalidColourException;

public class FormatVCF {

	public FormatVCF(){
				
	}

	/** Format alternative allele for text printing */
	public static char textForVariant(String refAllele, String altAllele){
		char text;
		if(refAllele.length() == 1 && altAllele.length() == 1){
			// SNP
			text= altAllele.charAt(0);
		} else if(refAllele.length() == 1 && altAllele.length() > 1){
			// Insertion into the reference
			text= 'I';
		} else if(refAllele.length() > 1 && altAllele.length() == 1){
			// Deletion into the reference
			text= 'D';
		} else {
			throw new RuntimeException();
		}
		return text; 
	}

	
	/** Format alternative allele for text printing 
	 * @throws InvalidColourException */
	public static String format(char textForVariant) throws InvalidColourException{

		// Format
		String formattedText= "";
		// For colour scheme see http://www.umass.edu/molvis/tutorials/dna/atgc.htm
		
		if(textForVariant == 'A' || textForVariant == 'a'){
			formattedText += "\033[38;5;231;48;5;" + Xterm256.colorNameToXterm256("blue") + "m" + textForVariant + "\033[0m";
		} else if(textForVariant == 'C' || textForVariant == 'c') {
			formattedText += "\033[38;5;231;48;5;" + Xterm256.colorNameToXterm256("red") + "m" + textForVariant + "\033[0m";
		} else if(textForVariant == 'G' || textForVariant == 'g') {
			formattedText += "\033[38;5;231;48;5;" + Xterm256.colorNameToXterm256("green") + "m" + textForVariant + "\033[0m";
		} else if(textForVariant == 'T' || textForVariant == 't') {
			formattedText += "\033[38;5;231;48;5;" + Xterm256.colorNameToXterm256("yellow") + "m" + textForVariant + "\033[0m";
		} else {
			formattedText += textForVariant;
		}
		return formattedText;
	}
	
}
