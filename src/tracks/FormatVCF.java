package tracks;

import org.apache.commons.lang3.StringUtils;

public class FormatVCF {

	public FormatVCF(){
				
	}
	
	/** Format alternative allele for text printing */
	public static String format(String refAllele, String altAllele, boolean noFormat){
		String text= "";
		if(refAllele.length() == 1 && altAllele.length() == 1){
			// SNP
			text= altAllele;
		} else if(refAllele.length() == 1 && altAllele.length() > 1){
			// Insertion into the reference
			text= "I";
		} else if(refAllele.length() > 1 && altAllele.length() == 1){
			// Deletion into the reference
			text= "D";
		}
		if(noFormat){
			return text; 
		}
		// Format
		String formattedText= "";
		for(int i= 0; i < text.length(); i++){
			// For colour scheme see http://www.umass.edu/molvis/tutorials/dna/atgc.htm
			char base= text.charAt(i);
			if(base == 'A' || base == 'a'){
				formattedText += "\033[107;34m" + base + "\033[0m";
			} else if(base == 'C' || base == 'c') {
				formattedText += "\033[107;31m" + base + "\033[0m";
			} else if(base == 'G' || base == 'g') {
				formattedText += "\033[107;32m" + base + "\033[0m";
			} else if(base == 'T' || base == 't') {
				formattedText += "\033[107;33m" + base + "\033[0m";
			} else {
				formattedText += base;
			} 
		}
		return formattedText;
	}
	
}
