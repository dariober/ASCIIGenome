package tracks;

import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

import com.google.common.collect.HashMultiset;
import com.google.common.collect.Multiset;

import samTextViewer.Utils;

/** Class to store information about a single position.
 * */
class Locus {

	String chrom;
	int pos;
	// Map<String, Integer> counts= new HashMap<String, Integer>();
	Multiset<String> counts = HashMultiset.create();

//	private static final int MIN_BASE_QUAL= 20; // Min base quality for a read base to be counted
	private static final int MIN_DEPTH_ALT= 3; // Min read depth for alternative allele to be taken into account
	private static final double MIN_PCT_ALT= 0.01; // Min % of alternative allele to be taken into account  
	private static final double MIN_PCT_TOT= 0.98; // % (Reference + Alternative) must be above this otherwise set consensus to N.   
	
	/*   C O N S T R U C T O R   */
	
	protected Locus(String chrom, int pos) {
		this.chrom= chrom;
		this.pos= pos;
	}

	/*  M E T H O D S  */
	protected void add(char base, boolean isReverse, boolean isFirstOFPair) {
		base= Character.toUpperCase(base);
		String key= null;
		if(base == 'A'){
			if(isFirstOFPair){
				if(isReverse){
					key= "A_1R";
				} else {
					key= "A_1F";
				}
			} else {
				if(isReverse){
					key= "A_2R";
				} else {
					key= "A_2F";
				}				
			}
		}
		else if(base == 'C'){
			if(isFirstOFPair){
				if(isReverse){
					key= "C_1R";
				} else {
					key= "C_1F";
				}
			} else {
				if(isReverse){
					key= "C_2R";
				} else {
					key= "C_2F";
				}				
			}
		}
		else if(base == 'G'){
			if(isFirstOFPair){
				if(isReverse){
					key= "G_1R";
				} else {
					key= "G_1F";
				}
			} else {
				if(isReverse){
					key= "G_2R";
				} else {
					key= "G_2F";
				}				
			}
		}
		else if(base == 'T'){
			if(isFirstOFPair){
				if(isReverse){
					key= "T_1R";
				} else {
					key= "T_1F";
				}
			} else {
				if(isReverse){
					key= "T_2R";
				} else {
					key= "T_2F";
				}				
			}
		}
		else if(base == 'N'){
			if(isFirstOFPair){
				if(isReverse){
					key= "N_1R";
				} else {
					key= "N_1F";
				}
			} else {
				if(isReverse){
					key= "N_2R";
				} else {
					key= "N_2F";
				}				
			}
		} 
		else if(base == 'D'){
			if(isFirstOFPair){
				if(isReverse){
					key= "D_1R";
				} else {
					key= "D_1F";
				}
			} else {
				if(isReverse){
					key= "D_2R";
				} else {
					key= "D_2F";
				}				
			}
		} else {
			throw new RuntimeException();
		}
		// int count= this.counts.get(key);	
		this.counts.add(key);
	}

	protected int getDepth(){
		int depth= 0;
		for(String key : this.counts.elementSet()){
			depth += this.counts.count(key);
		}
		return depth; 
	}

	private Map<Character, Integer> getGroupedCounts(){
		Map<Character, Integer> grpCounts= new HashMap<Character, Integer>();
		grpCounts.put('A', this.counts.count("A_1F") + this.counts.count("A_1R") + this.counts.count("A_2F") + this.counts.count("A_2R"));
		grpCounts.put('C', this.counts.count("C_1F") + this.counts.count("C_1R") + this.counts.count("C_2F") + this.counts.count("C_2R"));
		grpCounts.put('G', this.counts.count("G_1F") + this.counts.count("G_1R") + this.counts.count("G_2F") + this.counts.count("G_2R"));
		grpCounts.put('T', this.counts.count("T_1F") + this.counts.count("T_1R") + this.counts.count("T_2F") + this.counts.count("T_2R"));
		grpCounts.put('N', this.counts.count("N_1F") + this.counts.count("N_1R") + this.counts.count("N_2F") + this.counts.count("N_2R"));
		grpCounts.put('D', this.counts.count("D_1F") + this.counts.count("D_1R") + this.counts.count("D_2F") + this.counts.count("D_2R"));
		return grpCounts;
	}
	
	/** Call consensus base based on calls
	 * */
	protected char getConsensus(){
		
		Map<Character, Integer> grpCounts = this.getGroupedCounts();
		
		// * Sort by count 
		Iterator<Character> iter = Utils.sortByValue(grpCounts).keySet().iterator();
		char allele1= iter.next();
		char allele2= iter.next();

		// Is allele2 supported by at least n calls? 
		// Is allele2 making up more than x % of the total?
		char consensus;
		if(grpCounts.get(allele1) == 0){
			consensus= ' ';
		}
		else if((float)(grpCounts.get(allele1) + grpCounts.get(allele2))/this.getDepth() < MIN_PCT_TOT){
			consensus= 'N';
		} else if(grpCounts.get(allele2) >= MIN_DEPTH_ALT 
				&& (float)grpCounts.get(allele2)/this.getDepth() >= MIN_PCT_ALT ){
			consensus= this.iupacAmbiguity(allele1, allele2);
		} else {
			consensus= allele1;
		}
		return consensus;
	}
	
	private char iupacAmbiguity(char x, char y){
		
		if((x == 'A' && y == 'G') || (x == 'G' && y == 'A')){ return 'R'; }
		if((x == 'C' && y == 'T') || (x == 'T' && y == 'C')){ return 'Y'; }
		if((x == 'G' && y == 'C') || (x == 'C' && y == 'G')){ return 'S'; }
		if((x == 'A' && y == 'T') || (x == 'T' && y == 'A')){ return 'W'; }
		if((x == 'G' && y == 'T') || (x == 'T' && y == 'G')){ return 'K'; }
		if((x == 'A' && y == 'C') || (x == 'C' && y == 'A')){ return 'M'; }
		return 'N';
	}

}
