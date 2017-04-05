package tracks;

import java.util.HashMap;
import java.util.Map;

/** Class to store information about a single position.
 * */
class Locus {

	String chrom;
	int pos;
	Map<String, Integer> counts= new HashMap<String, Integer>();

	/*   C O N S T R U C T O R   */
	
	protected Locus(String chrom, int pos) {
		this.chrom= chrom;
		this.pos= pos;
		
//		// Counts for first in pair
		this.counts.put("A_1F", 0);
		this.counts.put("C_1F", 0);
		this.counts.put("G_1F", 0);
		this.counts.put("T_1F", 0);
		this.counts.put("N_1F", 0);
		this.counts.put("D_1F", 0);
		
		this.counts.put("A_1R", 0);
		this.counts.put("C_1R", 0);
		this.counts.put("G_1R", 0);
		this.counts.put("T_1R", 0);
		this.counts.put("N_1R", 0);
		this.counts.put("D_1R", 0);
		
//		// Counts for second in pair
		this.counts.put("A_2F", 0);
		this.counts.put("C_2F", 0);
		this.counts.put("G_2F", 0);
		this.counts.put("T_2F", 0);
		this.counts.put("N_2F", 0);
		this.counts.put("D_2F", 0);
		
		this.counts.put("A_2R", 0);
		this.counts.put("C_2R", 0);
		this.counts.put("G_2R", 0);
		this.counts.put("T_2R", 0);
		this.counts.put("N_2R", 0);
		this.counts.put("D_2R", 0);
		
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
		if(base == 'C'){
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
		if(base == 'G'){
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
		if(base == 'T'){
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
		if(base == 'N'){
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
		if(base == 'D'){
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
		} 
		int count= this.counts.get(key);
		this.counts.put(key, count+1);
	}

	protected int getDepth(){
		int depth= 0;
		for(int x : this.counts.values()){
			depth += x;
		}
		return depth;
	}
	
	protected Map<String, Integer> getCounts(){
		return this.counts;
	}
}
