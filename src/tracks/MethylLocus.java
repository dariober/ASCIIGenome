package tracks;

import filter.ReadFromTopStrandFilter;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
// import htsjdk.samtools.util.SamLocusIterator.LocusInfo;
import samTextViewer.SamLocusIterator.LocusInfo;
// import htsjdk.samtools.util.SamLocusIterator.RecordAndOffset;
import samTextViewer.SamLocusIterator.RecordAndOffset;

/**
 * Info at a cytosine site necessary for methylation calling
 */
class MethylLocus {

	private LocusInfo locus;
	// private boolean isCorGLocus= false;
	private char refBase; // This going to be alwasy upper case
	private Integer cntM= 0;
	private Integer cntU= 0;
//	private int cntN;
	
	             /* C o n s t r u c t o r s */
	
	public MethylLocus(LocusInfo locus, IndexedFastaSequenceFile faSeqFile){
		
		String chrom= locus.getSequenceName();
		this.refBase= Character.toUpperCase((char)faSeqFile.getSubsequenceAt(chrom, locus.getPosition(), locus.getPosition()).getBases()[0]);
		
		makeMethylLocus(locus, refBase);
	}
	
	/**
	 * Constructor based on locusInfo and reference sequence as byte[].
	 * Offset is the genomic position of the first base in subRefSeq. I.e. offset=1 if
	 * subRefSeq starts from the start of the chrom. offset=n if byte[] is a substring
	 * of reference chromosome starting at n.
	 * */
	public MethylLocus(LocusInfo locus, byte[] subRefSeq, int offset){
		
		if(locus.getPosition() < offset || locus.getPosition() > offset + subRefSeq.length){
			System.err.println("Required position is not in the supplied reference sequence.");
		}
		int posOnSubRefSeq= (locus.getPosition() - offset);	
		this.refBase= Character.toUpperCase((char)subRefSeq[posOnSubRefSeq]);
		makeMethylLocus(locus, this.refBase);	
	}
	
	                   /* M e t h o d s */
	/**
	 * Fill up fields cntM, cntU etc
	 */
	private void populateMethylLocus(){
		
		if(this.refBase != 'C' && this.refBase != 'G'){
			cntM= null;
			cntU= null;
			return;
		} 
		
		for(RecordAndOffset recOff : this.locus.getRecordAndPositions()){

			boolean readIsTopStrand= !(new ReadFromTopStrandFilter(true)).filterOut(recOff.getRecord());
			char readBase= Character.toUpperCase((char)recOff.getReadBase());
			
			if(this.refBase == 'C'){

				if( readIsTopStrand	){
					if(readBase == 'C'){
						cntM++;
					} else if(readBase == 'T'){
						cntU++;
					} 
				}  					
			} else if (this.refBase == 'G'){

				if(	!readIsTopStrand ){
					if(readBase == 'G'){
						cntM++;
					} else if(readBase == 'A'){
						cntU++;
					} 
				}  
			} else {
				System.err.println("Unexpected strand or base!");
				System.exit(1);
			}
		}		
	}

	private void makeMethylLocus(LocusInfo locus, char refBase){

		this.locus= locus;
		this.populateMethylLocus();
	}

	public String toString(){
		StringBuilder sb= new StringBuilder();
		sb.append("Locus: " + this.locus + "\n");
		sb.append("Reference base: " + this.refBase + "\n");
		sb.append("Count methylated (C/T): " + this.cntM + "\n");
		sb.append("Count unmethylated (G/A): " + this.cntU + "\n");
		return sb.toString();
	}
	
	        /* S e t t e r s   and   G e t t e r s */
	
	public LocusInfo getLocus() {
		return locus;
	}

	public void setLocus(LocusInfo locus) {
		this.locus = locus;
	}

	public Integer getCntM() {
		return cntM;
	}

	public Integer getCntU() {
		return cntU;
	}

	public char getRefBase() {
		return refBase;
	}
}
