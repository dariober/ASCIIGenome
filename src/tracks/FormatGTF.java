package tracks;

import java.util.HashMap;
import java.util.LinkedHashSet;
import java.util.Set;

import coloring.Config;
import coloring.ConfigKey;
import exceptions.InvalidColourException;

/** Mapping of GTF/GFF features to text character to use to represent them.
 * */
public class FormatGTF {

	private static HashMap<Character, HashMap<String, Character>> featureToTextCharDict= 
			new HashMap<Character, HashMap<String, Character>>();

	private static final Set<String> txSuperFeatures= new LinkedHashSet<String>();
	private static final Set<String> txSubFeatures= new LinkedHashSet<String>();
	
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
	 * @throws InvalidColourException 
	 * */
	public static String format(char text, char strand) throws InvalidColourException{
		StringBuilder sb= new StringBuilder();
		sb.append("\033[48;5;");
		if(strand == '+') {
			sb.append(Config.get256Color(ConfigKey.feature_background_positive_strand));
			sb.append(";38;5;");
			sb.append(Config.get256Color(ConfigKey.foreground));
		} else if(strand == '-') {
			sb.append(Config.get256Color(ConfigKey.feature_background_negative_strand));
			sb.append(";38;5;");
			sb.append(Config.get256Color(ConfigKey.foreground));
		} else {
			sb.append(Config.get256Color(ConfigKey.feature_background_no_strand));
			sb.append(";38;5;");
			sb.append(Config.get256Color(ConfigKey.foreground));
		}	
		sb.append("m");
		sb.append(text);
		sb.append("\033[48;5;");
		sb.append(Config.get256Color(ConfigKey.background));
		sb.append("m");
		return sb.toString();
	}

	protected static Set<String> getTxSuperFeatures() {

		// Features that define a record as a transcript:
		// Manually extracted from ensembl Homo_sapiens.GRCh38.86.chromosome.7.gff3.gz  
		txSuperFeatures.add("mrna");
		txSuperFeatures.add("transcript");
		txSuperFeatures.add("processed_transcript");
		txSuperFeatures.add("aberrant_processed_transcript");
		txSuperFeatures.add("NMD_transcript_variant");
		txSuperFeatures.add("pseudogenic_transcript");
		txSuperFeatures.add("lincrna");

		return txSuperFeatures;
	}

	protected static Set<String> getTxSubFeatures() {
		
		// Features that make part of a transcript.
		// Order matters: Put first the features that should be overwritten on screen by later features. 
		txSubFeatures.add("intron");
		txSubFeatures.add("exon");
		txSubFeatures.add("utr");
		txSubFeatures.add("5utr");
		txSubFeatures.add("five_prime_utr");
		txSubFeatures.add("3utr");
		txSubFeatures.add("three_prime_utr");
		txSubFeatures.add("cds");
		txSubFeatures.add("start_codon");
		txSubFeatures.add("stop_codon");
		
		return txSubFeatures;
	}
}
