package tracks;

import org.apache.commons.lang3.StringUtils;

import com.google.common.base.Joiner;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;

public class FormatVCF {

	public FormatVCF(VariantContext ctx){

	}

	/** Format alternative allele for text printing 
	 * */
	public static String textForVariant(VariantContext ctx){
		// Keep it simple, get the first alternate allele if more than one exist.
		if(ctx.isSNP()){
//			// SNV: Return the ALT base
			Allele alt = ctx.getAlleles().get(1);
			return alt.getDisplayString(); // ? or: alt.getBaseString();
		} 
		else if(ctx.isSimpleDeletion()) {
			Allele alt = ctx.getAlleles().get(0);
			return StringUtils.repeat("D", alt.length());
		} 
		else if(ctx.isSimpleInsertion()) {
			Allele alt = ctx.getAlleles().get(1);
			return StringUtils.repeat("I", alt.length());
		} 
		else {
			Allele alt = ctx.getAlleles().get(0);
			return StringUtils.repeat("|", alt.length());
		} 
//		// * DEL: Return the string corresponding to the deleted bases
//		
//		String altAllele= Joiner.on(",").join(ctx.getAlternateAlleles());
//		char text;
//		if(refAllele.length() == 1 && altAllele.length() == 1){
//			// SNP
//			text= altAllele.charAt(0);
//		} else if(altAllele.length() > refAllele.length()){
//			// Insertion into the reference
//			text= 'I';
//		} else if(refAllele.length() > altAllele.length()){
//			// Deletion into the reference
//			text= 'D';
//		} else {
//			throw new RuntimeException();
//		}
//		return text; 
	}

	
	/** Format alternative allele for text printing 
	 * @throws InvalidColourException */
//	@Deprecated
//	public static String format(char textForVariant) throws InvalidColourException{
//
//		// Format
//		StringBuilder formattedText= new StringBuilder();
//		formattedText.append("\033[7;48;5;"); // 7: Invert colour back/fore-ground
//		formattedText.append(Config.get256Color(ConfigKey.background));
//		formattedText.append(";38;5;");
//		if(textForVariant == 'A' || textForVariant == 'a'){
//			formattedText.append(Config.get256Color(ConfigKey.seq_a));
//			// formattedText += "\033[38;5;231;48;5;" + Xterm256.colorNameToXterm256("blue") + "m" + textForVariant + "\033[0m";
//		} else if(textForVariant == 'C' || textForVariant == 'c') {
//			formattedText.append(Config.get256Color(ConfigKey.seq_c));
//			// formattedText += "\033[38;5;231;48;5;" + Xterm256.colorNameToXterm256("red") + "m" + textForVariant + "\033[0m";
//		} else if(textForVariant == 'G' || textForVariant == 'g') {
//			formattedText.append(Config.get256Color(ConfigKey.seq_g));
//			// formattedText += "\033[38;5;231;48;5;" + Xterm256.colorNameToXterm256("green") + "m" + textForVariant + "\033[0m";
//		} else if(textForVariant == 'T' || textForVariant == 't') {
//			formattedText.append(Config.get256Color(ConfigKey.seq_t));
//			// formattedText += "\033[38;5;231;48;5;" + Xterm256.colorNameToXterm256("yellow") + "m" + textForVariant + "\033[0m";
//		} else {
//			formattedText.append(Config.get256Color(ConfigKey.seq_other));
//			// formattedText += textForVariant;
//		}
//		formattedText.append("m");
//		formattedText.append(textForVariant);
//		formattedText.append("\033[0m\033[38;5;");
//		formattedText.append(Config.get256Color(ConfigKey.foreground));
//		formattedText.append(";48;5;");
//		formattedText.append(Config.get256Color(ConfigKey.background));
//		formattedText.append("m");
//		return formattedText.toString();
//	}
	
}
