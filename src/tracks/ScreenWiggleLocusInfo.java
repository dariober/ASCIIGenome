package tracks;

/** Info about a screen position for wiggle-like data. Somewhat counterpart to ScreenLocusInfo */
public class ScreenWiggleLocusInfo {
	
	private int cntGenomicLoci= 0; // Count of genomic loci mapped to this screen position
	private float sumScore= 0; // Sum of scores accumulated from wiggle sites mapped to this screen locus
	
	/* C o n s t r u c t o r */
	public ScreenWiggleLocusInfo(){ }
	
	/* M e t h o d s */
	/** Increment attributes by given score */
	public void increment(float score){
		cntGenomicLoci++;
		sumScore += score;
	}

	public String toString(){
		String str= "cntGenomicLoci: " + this.cntGenomicLoci + "; sumScores: " + this.sumScore;
		return str;
	}
	
	/*   G e t t e r s   */
	protected float getMeanScore(){
		return (float)this.sumScore / this.cntGenomicLoci; 
	}
}
