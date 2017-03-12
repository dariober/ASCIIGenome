package tracks;

import coloring.Config;
import coloring.ConfigKey;
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
		StringBuilder formattedText= new StringBuilder();
		formattedText.append("\033[7;48;5;"); // 7: Invert colour back/fore-ground
		formattedText.append(Config.getColor(ConfigKey.background));
		formattedText.append(";38;5;");
		if(textForVariant == 'A' || textForVariant == 'a'){
			formattedText.append(Config.getColor(ConfigKey.seq_a));
			// formattedText += "\033[38;5;231;48;5;" + Xterm256.colorNameToXterm256("blue") + "m" + textForVariant + "\033[0m";
		} else if(textForVariant == 'C' || textForVariant == 'c') {
			formattedText.append(Config.getColor(ConfigKey.seq_c));
			// formattedText += "\033[38;5;231;48;5;" + Xterm256.colorNameToXterm256("red") + "m" + textForVariant + "\033[0m";
		} else if(textForVariant == 'G' || textForVariant == 'g') {
			formattedText.append(Config.getColor(ConfigKey.seq_g));
			// formattedText += "\033[38;5;231;48;5;" + Xterm256.colorNameToXterm256("green") + "m" + textForVariant + "\033[0m";
		} else if(textForVariant == 'T' || textForVariant == 't') {
			formattedText.append(Config.getColor(ConfigKey.seq_t));
			// formattedText += "\033[38;5;231;48;5;" + Xterm256.colorNameToXterm256("yellow") + "m" + textForVariant + "\033[0m";
		} else {
			formattedText.append(Config.getColor(ConfigKey.seq_other));
			// formattedText += textForVariant;
		}
		formattedText.append("m");
		formattedText.append(textForVariant);
		formattedText.append("\033[0m\033[38;5;");
		formattedText.append(Config.getColor(ConfigKey.foreground));
		formattedText.append(";48;5;");
		formattedText.append(Config.getColor(ConfigKey.background));
		formattedText.append("m");
		return formattedText.toString();
	}
	
}
