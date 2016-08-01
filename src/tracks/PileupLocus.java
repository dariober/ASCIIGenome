package tracks;

import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.Map.Entry;

import samTextViewer.SamLocusIterator.LocusInfo;
import samTextViewer.SamLocusIterator.RecordAndOffset;
import samTextViewer.Utils;

public class PileupLocus {

	private static final int MIN_BASE_QUAL= 20; // Min base quality for a read base to be counted
	private static final int MIN_DEPTH_ALT= 3; // Min read depth for alternative allele to be taken into account
	private static final double MIN_PCT_ALT= 0.01; // Min % of alternative allele to be taken into account  
	private static final double MIN_PCT_TOT= 0.98; // % (Reference + Alternative) must be above this otherwise set consensus to N.   
	
	String chrom;
	int pos= -1;
	char ref= '.';
	private LinkedHashMap<Character, Integer> baseCount;
	
	/* C o n s t r u c t o r  */
	
	/** Populate the ACTG counts for the given locus info. ref is the reference base if available.  
	 * The LocusInfo object is expected to be already filtered for desired flags. Note threshold 
	 * to mapping quality to count bases. */
	PileupLocus(LocusInfo locus, char ref){
		
		this.chrom= locus.getSequenceName();
		this.pos= locus.getPosition();
		this.ref= ref;

		this.baseCount= new LinkedHashMap<Character, Integer>();
		this.baseCount.put('A', 0);
		this.baseCount.put('C', 0);
		this.baseCount.put('G', 0);
		this.baseCount.put('T', 0);
		this.baseCount.put('N', 0);

		
		for(RecordAndOffset recOff : locus.getRecordAndPositions()){
			if(((int)recOff.getBaseQuality()) > PileupLocus.MIN_BASE_QUAL){
				char base = Character.toUpperCase((char) recOff.getReadBase());
				int count= this.baseCount.get(base) + 1;
				this.baseCount.put(base, count);
			}
		}
	}
	
	/* M e t h o d s */
	
	/** Call consensus base based on calls
	 * */
	protected char getConsensus(){
		// * Sort by count 
		
		Iterator<Character> iter = Utils.sortByValue(this.baseCount).keySet().iterator();
		char allele1= iter.next();
		char allele2= iter.next();

		// Is allele2 supported by at least n calls? 
		// Is allele2 making up more than x % of the total?
		char consensus;
		if(this.baseCount.get(allele1) == 0){
			consensus= ' ';
		}
		else if((float)(this.baseCount.get(allele1) + this.baseCount.get(allele2))/this.depth() < PileupLocus.MIN_PCT_TOT){
			consensus= 'N';
		} else if(this.baseCount.get(allele2) >= PileupLocus.MIN_DEPTH_ALT 
				&& (float)this.baseCount.get(allele2)/this.depth() >= PileupLocus.MIN_PCT_ALT ){
			consensus= this.iupacAmbiguity(allele1, allele2);
		} else {
			consensus= allele1;
		}
		// Is the consensus the same as reference?
		if(this.ref == consensus){
			return '=';
		} else {
			return consensus; 
		}
	}
	
	private char iupacAmbiguity(char x, char y){
		
		if((x == 'A' && y == 'G') || (x == 'G' && y == 'A')){ return 'R'; }
		if((x == 'C' && y == 'T') || (x == 'T' && y == 'C')){ return 'Y'; }
		if((x == 'G' && y == 'C') || (x == 'C' && y == 'G')){ return 'S'; }
		if((x == 'A' && y == 'T') || (x == 'T' && y == 'A')){ return 'W'; }
		if((x == 'G' && y == 'T') || (x == 'T' && y == 'G')){ return 'K'; }
		if((x == 'A' && y == 'C') || (x == 'C' && y == 'A')){ return 'M'; }
		throw new RuntimeException();
	}
	
	private int depth(){
		int depth= 0;
		
		Iterator<Entry<Character, Integer>> iter = this.baseCount.entrySet().iterator();
		while(iter.hasNext()){
			depth += iter.next().getValue();
		}
		return depth; 
	}
	
	public String toString(){
		StringBuilder sb= new StringBuilder();
		sb.append(this.chrom); sb.append("\t");
		sb.append(this.pos); sb.append("\t");
		sb.append(this.ref); sb.append("\t");
		sb.append(this.baseCount.get('A')); sb.append("\t");
		sb.append(this.baseCount.get('C')); sb.append("\t");
		sb.append(this.baseCount.get('G')); sb.append("\t");
		sb.append(this.baseCount.get('T')); sb.append("\t");
		sb.append(this.baseCount.get('N'));
		return sb.toString();
	}	
}
