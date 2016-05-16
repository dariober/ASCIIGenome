package tracks;

import filter.ReadFromTopStrandFilter;
import samTextViewer.SamLocusIterator.LocusInfo;
import samTextViewer.SamLocusIterator.RecordAndOffset;

/**
 * Store info about a locus on screen. Sort of counterpart of LocusInfo which refers to genomic
 * position, here we refer to sreen postion */
class ScreenLocusInfo {
	
	// private int screenPosition= -1; 
	private int cntGenomicLoci= 0; // Count of genomic loci mapped to this screen position
	private int depth= 0; // Sum of depths of all the loci mapped at this screen pos.
	private int cntM= 0; // For BS-Seq: Count of unconverted cytosines (methylated)
	private int cntU= 0; // Same as cntM, but for converted cyt
	
	/* Constructor */
	protected ScreenLocusInfo(){ }
	
	/* Methods */
	/** Increment attributes using locusInfo and refSeq. refSeq can be null which means
	 * some attributes are not updated */
	protected void increment(LocusInfo locusInfo, byte refBase, boolean bs){

		refBase= (byte) Character.toUpperCase(refBase);
		
		cntGenomicLoci++;
		depth += locusInfo.getRecordAndPositions().size();
		if(! bs){ 
			return;
		}
		
		// Code for BS-Seq data:
		if(refBase == '\0' || (refBase != 'C' && refBase != 'G')){
			return;
		}
		for(RecordAndOffset recOff : locusInfo.getRecordAndPositions()){

			boolean readIsTopStrand= !(new ReadFromTopStrandFilter(true)).filterOut(recOff.getRecord());
			char readBase= Character.toUpperCase((char)recOff.getReadBase());
			
			if(refBase == 'C'){

				if( readIsTopStrand	){
					if(readBase == 'C'){
						cntM++;
					} else if(readBase == 'T'){
						cntU++;
					} 
				}  					
			} else if (refBase == 'G'){

				if(	!readIsTopStrand ){
					if(readBase == 'G'){
						cntM++;
					} else if(readBase == 'A'){
						cntU++;
					} 
				}  
			} else {
				System.err.println("Unexpected strand or base!");
				throw new RuntimeException();
			}
		} //end iterating locusInfo
	}
	
	public String toString(){
		String str= "cntGenomicLoci: " + this.cntGenomicLoci + "; depth: " + 
				this.depth + "; cntM: " + this.cntM + "; cntU: " + this.cntU;
		return str;
	}

	/*   G e t t e r s   */
	protected Double getMeanDepth(){
		return (double)this.depth / this.cntGenomicLoci; 
	}
	protected Double getMeanCntM(){
		return (double)this.cntM / this.cntGenomicLoci; 
	}
	protected Double getMeanCntU(){
		return (double)this.cntU / this.cntGenomicLoci; 
	}
	protected int getCntM(){
		return cntM; 
	}
	protected int getCntU(){
		return cntU; 
	}
	
}
