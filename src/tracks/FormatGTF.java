package tracks;

import java.util.HashMap;

/** Mapping of GTF/GFF features to text character to use to represent them.
 * */
public class FormatGTF {

	private static HashMap<Character, HashMap<String, Character>> featureToTextCharDict= 
			new HashMap<Character, HashMap<String, Character>>();
	
	/* C O N S T R U C T O R */
	
	public FormatGTF(){

	}

	/* G E T T E R S */
	
	public static HashMap<Character, HashMap<String, Character>> getFeatureToTextCharDict(){
		
		/* Map GTF features to characters. Forward capital LETTERS, reverse small letters  
		 * Feature names are case insensitive */
		HashMap<String, Character> fwdFeature= new HashMap<String, Character>();
		HashMap<String, Character> revFeature= new HashMap<String, Character>();
		HashMap<String, Character> unstrFeature= new HashMap<String, Character>();
		
		// MEMO: Makr these case insensitive by writing them lowercase.
		fwdFeature.put("exon",            'E'); revFeature.put("exon",            'e');
		fwdFeature.put("cds", 	          'C'); revFeature.put("cds",             'c');
		fwdFeature.put("start_codon",     'A'); revFeature.put("start_codon",     'a');
		fwdFeature.put("stop_codon",      'Z'); revFeature.put("stop_codon",      'z');
		fwdFeature.put("utr",             'U'); revFeature.put("utr",             'u');
		fwdFeature.put("3utr",            'U'); revFeature.put("3utr",            'u');
		fwdFeature.put("three_prime_utr", 'U'); revFeature.put("three_prime_utr", 'u');
		fwdFeature.put("5utr",           'W'); revFeature.put("5utr",            'w');
		fwdFeature.put("five_prime_utr", 'W'); revFeature.put("five_prime_utr",  'w');
		fwdFeature.put("gene",		     'G'); revFeature.put("gene",            'g');
		fwdFeature.put("transcript",     'T'); revFeature.put("transcript",      't');
		fwdFeature.put("mrna",           'M'); revFeature.put("mrna",            'm');
		fwdFeature.put("trna",           'X'); revFeature.put("trna",            'x');
		fwdFeature.put("rrna", 		     'R'); revFeature.put("rrna",            'r');
		fwdFeature.put("mirna",          'I'); revFeature.put("mirna",           'i');
		fwdFeature.put("ncrna",          'L'); revFeature.put("ncrna",           'l');
		fwdFeature.put("lncrna",         'L'); revFeature.put("lncrna",          'l');
		fwdFeature.put("sirna",          'S'); revFeature.put("sirna",           's');
		fwdFeature.put("pirna",          'P'); revFeature.put("pirna",           'p');
		fwdFeature.put("snorna",         'O'); revFeature.put("snorna",          'o');
		
		// For feature with strand not available, use forward encoding, unless feature unknown
		unstrFeature.putAll(fwdFeature);
		fwdFeature.put("other", '>'); 		
		revFeature.put("other", '<');
		unstrFeature.put("other", '|');
		
		featureToTextCharDict.put('+', fwdFeature);
		featureToTextCharDict.put('-', revFeature);
		featureToTextCharDict.put('.', unstrFeature);
		
		return featureToTextCharDict; 
	}
	
	/* M E T H O D S */
	
	/** Return text formatted according to strand. 
	 * */
	public static String format(char text, char strand){
		if(strand == '+') {
			return "\033[30;48;5;147m" + text + "\033[48;5;231m";
		} else if(strand == '-') {
			return "\033[30;48;5;225m" + text + "\033[48;5;231m";
		} else {
			return "\033[30;47m" + text + "\033[48;5;231m";
		}	
	}
}
