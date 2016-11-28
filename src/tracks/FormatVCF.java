package tracks;

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

	
	/** Format alternative allele for text printing */
	public static String format(char textForVariant){

		// Format
		String formattedText= "";
		// For colour scheme see http://www.umass.edu/molvis/tutorials/dna/atgc.htm
		
		if(textForVariant == 'A' || textForVariant == 'a'){
			formattedText += "\033[107;34m" + textForVariant + "\033[48;5;231m";
		} else if(textForVariant == 'C' || textForVariant == 'c') {
			formattedText += "\033[107;31m" + textForVariant + "\033[48;5;231m";
		} else if(textForVariant == 'G' || textForVariant == 'g') {
			formattedText += "\033[107;32m" + textForVariant + "\033[48;5;231m";
		} else if(textForVariant == 'T' || textForVariant == 't') {
			formattedText += "\033[107;33m" + textForVariant + "\033[48;5;231m";
		} else {
			formattedText += textForVariant;
		}
		return formattedText;
	}
	
}
