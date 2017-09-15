package tracks;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import coloring.Config;
import coloring.ConfigKey;
import exceptions.InvalidColourException;
import exceptions.InvalidGenomicCoordsException;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import samTextViewer.GenomicCoords;
import samTextViewer.Utils;

/** Class holding features of a read necessary to represent the read in text
 * format.
 * 
 * Nomenclature: dna*:  Bases as ACTG
 *               cons*: Bases to highlight consensus to reference, e.g. '.', ','
 *               text*: Bases as they will be printed on screen.
 * 
 * Missing/Todo: Insertions to the reference are not visible.
 * @author berald01
          10        20        30        40        
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
AACCTTGGCC---------------------------------------
-AACCTTGGCC--------------------------------------
--AACCTTGGCC-------------------------------------
-------------AACCTTGGCC--------------------------
--------------AA----CCTT-------------------------
 */
class TextRead extends IntervalFeature{
	
	// Characters for methylation coding. NB: M and U are valid DNA chars in IUPAC!
	private final static char charM= 'M';
	private final static char charU= 'U';
	private final static char charm= 'm';
	private final static char charu= 'u';
	private final static char charFwd= '>';
	private final static char charRev= '<';
	private final static int  SHADE_MAPQ= 5;
	// private final static int  SHADE_BASEQ= 13;
	
	/** Char to represent deletions from the reference. I.e. gaps in the read */
	final private char DEL= '-';
	/** Char to represent region skip. I.e. gaps in the read */
	final private char SKIP= '_';
	
	/** Start position of the read in window coordinates. 1-based. If in genomic
	 * coords read starts aligning at pos 100 and window span is 100:150, then textStart= 1.*/
	private int textStart;

	/** End position of the read in window coordinates. 1-based. If in genomic
	 * coords read ends pos 150 and window span is 100:150, then textEnd= 1. */
	private int textEnd;
	
	private SAMRecord samRecord;
	private GenomicCoords gc;
	private List<int[]> textPositionsOfSkippedBases= new ArrayList<int[]>();
	
	/*    C o n s t r u c t o r s    */
		
	protected TextRead(SAMRecord rec, GenomicCoords gc) throws InvalidGenomicCoordsException, IOException{
		// At least part of the read must be in the window
		//            |  window  |
		//                         |------| read
		if(rec.getAlignmentStart() > gc.getTo()){
			System.err.println("Alignment starts beyond text window!");
			System.err.println(rec.getSAMString());
			System.err.println("Aln starts: " +  rec.getAlignmentStart());
			System.err.println("Aln ends: " +  rec.getAlignmentEnd());
			System.err.println("Window: " + gc.toString());
			throw new RuntimeException();
		}
		//            |  window  |
		//|---------| read
		if(rec.getAlignmentEnd() < gc.getFrom()){
			System.err.println("Alignment ends before text window!");
			System.err.println(rec.getSAMString());
			System.err.println("Aln starts: " +  rec.getAlignmentStart());
			System.err.println("Aln ends: " +  rec.getAlignmentEnd());
			System.err.println("Window: " + gc.toString());
			throw new RuntimeException();
		}		
		this.gc= gc;
		this.samRecord= rec;
		this.setTextStart();
		this.setTextEnd();
		this.setTextPositionsOfSkippedBases();
	}
	
	/*       M e t h o d s       */

	/**
	 * Return read ready to be printed on track. 
	 * @param refSeq Reference sequence. Can be null in which case bases are displayed as they are.
	 * @param bs Should the read be converted to BS-Seq mode?
	 * @param noFormat Do not apply string formatting (colours, etc.)
	 * @param withReadName Print the read name instead of the bases.
	 * @return
	 * @throws IOException 
	 * @throws InvalidGenomicCoordsException 
	 * @throws InvalidColourException 
	 */
	public String getPrintableTextRead(boolean bs, boolean noFormat, boolean withReadName) throws IOException, InvalidGenomicCoordsException, InvalidColourException{
		List<FeatureChar> fmt = this.getTextReadAsFeatureChars(bs);
		StringBuilder sb= new StringBuilder();
		for(FeatureChar x : fmt){
			sb.append(x.format(noFormat));
		}
		return sb.toString();
//		if(withReadName){ // Replace bases with read name. As long as it fits
//			int upTo= (unformatted.size() < this.samRecord.getReadName().length()) 
//					? unformatted.size() 
//					: this.samRecord.getReadName().length();  
//			for(int i= 0; i< upTo; i++){
//				char nchar= (char) this.samRecord.getReadName().getBytes()[i];
//				unformatted.set(i, nchar);
//			}
//			if(unformatted.size() > upTo){ // Add a separator btw read name and sequence
//				unformatted.set(upTo, '/');
//			}
//			return Joiner.on("").join(unformatted);
//		} else {
//			return readFormatter(unformatted, noFormat, bs);
//		}
	}
	
	/**
	 * Obtain the start position of the read on screen. Screen positions are 1-based. 
	 * @return
	 * @throws IOException 
	 * @throws InvalidGenomicCoordsException 
	 */
	private void setTextStart() throws InvalidGenomicCoordsException, IOException{		
		if(samRecord.getAlignmentStart() <= gc.getFrom()){ // Read starts right at the window start or even earlier
			this.textStart= 1;
			return;
		}		
		this.textStart= Utils.getIndexOfclosestValue(samRecord.getAlignmentStart(), gc.getMapping()) + 1;
		return;
	}
	
	private void setTextEnd() throws InvalidGenomicCoordsException, IOException{
		this.textEnd= Utils.getIndexOfclosestValue(samRecord.getAlignmentEnd(), gc.getMapping()) + 1;
		return;
	}
	
	/**List of positions on screen where the skipped bases (cigar op: N) start
	 * @throws IOException 
	 * @throws InvalidGenomicCoordsException 
	 * */
	private void setTextPositionsOfSkippedBases() throws InvalidGenomicCoordsException, IOException{
		int genomicPosition= this.getSamRecord().getAlignmentStart();
		List<CigarElement> cigar = this.getSamRecord().getCigar().getCigarElements();
		for(CigarElement el : cigar){
			if(el.getOperator().equals(CigarOperator.SKIPPED_REGION)){
				int[] textPositions= new int[2];
				// +1 because textPosition is 1-based
				textPositions[0]= Utils.getIndexOfclosestValue(genomicPosition, this.gc.getMapping()) + 1; 
				textPositions[1]= Utils.getIndexOfclosestValue(genomicPosition + el.getLength(), this.gc.getMapping()) + 1;
				this.textPositionsOfSkippedBases.add(textPositions);
			};
			if(el.getOperator().consumesReferenceBases()){
				genomicPosition += el.getLength();
			}
		}
	}

	/** If the windowSize and genomic span are not mapped 1:1, i.e. 1 bp : 1 char, then
	 * represent reads as simplified bases */
	private List<FeatureChar> getSquashedRead(){
		
		ArrayList<FeatureChar> squashedRead= new ArrayList<FeatureChar>();
		char xc;
		if(this.samRecord.getReadNegativeStrandFlag()){
			xc= charRev;
		} else {
			xc= charFwd;
		}
		for(int i= this.textStart; i <= this.textEnd; i++){
			FeatureChar sq= new FeatureChar();
			// Set char to print
			if(this.textPositionIsSkipped(i)){
				sq.setText(this.SKIP);
			} else {
				sq.setText(xc);
			}
			// Set formatting
			if(!this.samRecord.getReadNegativeStrandFlag()){
				sq.setFgColor(Config.get(ConfigKey.feature_background_positive_strand));
			} 
			else if(this.samRecord.getReadNegativeStrandFlag()){
				sq.setFgColor(Config.get(ConfigKey.feature_background_negative_strand));
			} 

			if(this.samRecord.getMappingQuality() < SHADE_MAPQ){
				sq.setBgColor(Config.get(ConfigKey.shade_low_mapq));
			}
			squashedRead.add(sq);
		}
		return squashedRead;
	}
	
	private boolean textPositionIsSkipped(int textPos){
		for(int[] skippedRegion : this.textPositionsOfSkippedBases){
			if(textPos >= skippedRegion[0] && textPos <= skippedRegion[1]){
				return true; 
			}
		}
		return false;
	}
	
	/** Get a representation of the read as it appears aligned to the reference. 
	 * I.e. clipped ends omitted and deletions appearing as gaps (empty byte).
	 * Only the portion contained between the genomic coords from:to is returned.
	 * @return
	 * @throws IOException 
	 * @throws InvalidGenomicCoordsException 
	 */
	protected List<FeatureChar> getTextReadAsFeatureChars(boolean bs) throws InvalidGenomicCoordsException, IOException {

		if( ! this.gc.isSingleBaseResolution){
			return this.getSquashedRead();
		}
		
		// Accumulate here the read bases inside the window 
		ArrayList<FeatureChar> dnaRead= new ArrayList<FeatureChar>();
		byte[] readBases= samRecord.getReadBases();
		byte[] baseQual= this.samRecord.getBaseQualities();
		
		byte[] ref= null;
		if(this.gc.getRefSeq() != null){
			ref= Arrays.copyOfRange(this.gc.getRefSeq(), this.textStart-1, this.textEnd);
		}
		int SHADE_BASEQ= Integer.parseInt(Config.get(ConfigKey.shade_baseq));
		
		// Walk along the aligned read and append bases to textRead as long as
		// the genomic position of the base is inside the genomic coords of the window
		int curBaseGenomicPos= samRecord.getAlignmentStart();
		int curBaseReadPos= 0; // Position on read. Start from zero walk along the read
		List<CigarElement> cigarEls= samRecord.getCigar().getCigarElements();
		for(CigarElement el : cigarEls){
			if(el.getOperator() == CigarOperator.M || el.getOperator() == CigarOperator.EQ || el.getOperator() == CigarOperator.X){
				// Add nucleotide chars to growing read
				for(int i= 0; i < el.getLength(); i++){
					FeatureChar xc= new FeatureChar();
					if(curBaseGenomicPos >= gc.getFrom() && curBaseGenomicPos <= gc.getTo()){
						// If base is inside window:
						if(readBases.length > 0){
							xc.setText((char)readBases[curBaseReadPos]);
						} else { // If sam record has no read seq stored put N
							xc.setText('N');
						}
						
						if(ref != null){
							char refBase= Character.toUpperCase((char) ref[dnaRead.size()]);
							if(bs){
								this.convertDnaBaseToTextBS(xc, refBase);
							}
							if(Character.toUpperCase(xc.getText()) == refBase){
								if(this.samRecord.getReadNegativeStrandFlag()){
									xc.setText(',');
								} else {
									xc.setText('.');
								}
							}
						}
						
						// Add formatting as appropriate
						if(this.samRecord.getMappingQuality() < SHADE_MAPQ){
							xc.setBgColor(Config.get(ConfigKey.shade_low_mapq));
							xc.setFgColor(Config.get(ConfigKey.foreground));
						}
						else if(bs && Character.toUpperCase(xc.getText()) == charM){
							xc.setBgColor(Config.get(ConfigKey.methylated_background));
							xc.setFgColor(Config.get(ConfigKey.methylated_foreground));
						} 
						else if(bs && Character.toUpperCase(xc.getText()) == charU){
							xc.setBgColor(Config.get(ConfigKey.unmethylated_background));
							xc.setFgColor(Config.get(ConfigKey.unmethylated_foreground));
						} 
						else if(Character.toUpperCase(xc.getText()) == 'A'){
							xc.setFgColor(Config.get(ConfigKey.seq_a));
						} 
						else if(Character.toUpperCase(xc.getText()) == 'C') {
							xc.setFgColor(Config.get(ConfigKey.seq_c));
						} 
						else if(Character.toUpperCase(xc.getText()) == 'G') {
							xc.setFgColor(Config.get(ConfigKey.seq_g));
						} 
						else if(Character.toUpperCase(xc.getText()) == 'T') {
							xc.setFgColor(Config.get(ConfigKey.seq_t));
						} 
						else if(!bs && !this.samRecord.getReadNegativeStrandFlag()){
								xc.setFgColor(Config.get(ConfigKey.feature_background_positive_strand));
						} 
						else if(!bs && this.samRecord.getReadNegativeStrandFlag()){
								xc.setFgColor(Config.get(ConfigKey.feature_background_negative_strand));
						}
						
						if(baseQual.length > 0){
							int bq= (int) baseQual[i];
							if(bq < SHADE_BASEQ){
								xc.setBgColor(Config.get(ConfigKey.shade_low_mapq));
							}
						}
						dnaRead.add(xc);
					}
					curBaseGenomicPos++; // M consumes read and ref bases. So increment them
					curBaseReadPos++;
				}
			} else if(el.getOperator() == CigarOperator.D || el.getOperator() == CigarOperator.N){
				// Add gap chars to growing read
				for(int i= 0; i < el.getLength(); i++){
					FeatureChar xc= new FeatureChar();
					if(curBaseGenomicPos >= gc.getFrom() && curBaseGenomicPos <= gc.getTo()){
						if(el.getOperator() == CigarOperator.D){
							xc.setText(this.DEL);
							xc.setInvertFgBgColor(true);
						} else if(el.getOperator() == CigarOperator.N){ 
							xc.setText(this.SKIP);
						} else {
							System.err.println("Unexpected operator");
							throw new RuntimeException();
						}
						dnaRead.add(xc);
					}
					curBaseGenomicPos++;
				}
			} else if(el.getOperator() == CigarOperator.I) {
				if(dnaRead.size() > 0){ // If the insertion is outside the terminal window, there is no base to mark
					dnaRead.get(dnaRead.size()-1).setInvertFgBgColor(true);
				}
				curBaseReadPos += el.getLength();
			} else if(el.getOperator() == CigarOperator.S){
				curBaseReadPos += el.getLength();
			} else if(el.getOperator() == CigarOperator.H){
				// Nothing to do
			} else if(el.getOperator() == CigarOperator.P){
				// Nothing to do: NOT SURE is is correct to just ignore padding!
			} else {
				System.err.println("Unexpected operator in cigar string for record\n" + samRecord.getSAMString()); 
				throw new RuntimeException();
			}
		}
		for(FeatureChar x : dnaRead){
			if(this.samRecord.getReadNegativeStrandFlag()){
				x.setText(Character.toLowerCase(x.getText()));
			} else {
				x.setText(Character.toUpperCase(x.getText()));
			}
		}
		return dnaRead;
	}
	
	/** Convert textRead, the actual bases found in sam, to represent match, mismatch and strandness.
	 * @param refSeq The reference sequence spanning the window.
	 * @return
	 * @throws IOException 
	 * @throws InvalidGenomicCoordsException 
	 */
//	private List<FeatureChar> getConsRead(boolean bs) throws IOException, InvalidGenomicCoordsException {
//		
//		List<FeatureChar> dnaRead= this.getDnaRead(bs);
//		if(this.gc.getRefSeq() == null){
//			return dnaRead;
//		}
//		List<FeatureChar> consRead= new ArrayList<FeatureChar>();
//		int posOnRead= 0;
//		for(int i= this.textStart - 1; i < this.textEnd; i++){
//			char base= Character.toUpperCase( dnaRead.get(posOnRead).getText() );
//			char ref= (char) Character.toUpperCase(this.gc.getRefSeq()[i]);
//			if( base == ref){
//				FeatureChar xc= new FeatureChar();
//				if(this.samRecord.getReadNegativeStrandFlag()){
//					xc.setText(',');
//				} else {
//					xc.setText('.');
//				}
//				consRead.add(xc);
//			} else {
//				FeatureChar xc= new FeatureChar();
//				if(this.samRecord.getReadNegativeStrandFlag()){
//					xc.setText(Character.toLowerCase(base));
//				} else {
//					xc.setText(Character.toUpperCase(base));
//				}	
//				consRead.add(xc);
//			}
//			posOnRead++;
//		}
//		return consRead;
//	}

	/** Memo: You compare the reads with these bases on the reference:
	   1st |  2nd
	-------+------
	+ve  C |  G
	-------+------
 	-ve  g |  c
 	
	 * @param refSeq
	 * @return
	 * @throws IOException 
	 * @throws InvalidGenomicCoordsException 
	 */
	private void convertDnaBaseToTextBS(FeatureChar dnaBase, char refBase) throws IOException, InvalidGenomicCoordsException{
	
		if(this.gc.getRefSeq() == null){ // Effectively don't convert 
			return;
		}
		
		// For convenience extract flags from sam record
		boolean isSecondOfPair= false;
		if(samRecord.getReadPairedFlag() && samRecord.getSecondOfPairFlag()){
			isSecondOfPair= true;
		}
		boolean isForwardStrand= !this.samRecord.getReadNegativeStrandFlag();
		
		char dnaChar= Character.toUpperCase(dnaBase.getText());
		if( ( isForwardStrand && !isSecondOfPair ) || ( !isForwardStrand && isSecondOfPair )){
			// Look for C on the reference
			if(isForwardStrand){ // +strand, first in pair or unpaired
				if(refBase == 'C' && dnaChar == 'C'){
					dnaBase.setText(charM);
				} else if(refBase == 'C' && dnaChar == 'T'){
					dnaBase.setText(charU);
				} else {
					// Nothing to change
				}
			} else { // -ve strand, 2nd in pair
				// Look for c=',' -> m; 't' -> u
				if(refBase == 'C' && dnaChar == 'C'){
					dnaBase.setText(charm);
				} else if(refBase == 'C' && dnaBase.getText() == 'T'){
					dnaBase.setText(charu);
				} else {
					// Nothing to change
				}
			}
		} else if( ( !isForwardStrand && !isSecondOfPair ) || ( isForwardStrand && isSecondOfPair )){
			// Look for G on the reference
			if(!isForwardStrand){ // -ve strand; first in pair or unpaired
				if(refBase == 'G' && dnaChar == 'G'){
					dnaBase.setText('m');
				} else if(refBase == 'G' && dnaChar == 'A'){
					dnaBase.setText('u');
				} else {
					// Nothing to change
				}
			} else { // -ve strand, 2nd in pair
				if(refBase == 'G' && dnaChar == 'G'){
					dnaBase.setText('M');
				} else if(refBase == 'G' && dnaChar == 'A'){
					dnaBase.setText('U');
				} else {
					// Nothing to change
				}
			}
		}
	}
	
	public String toString(double bpPerScreenColumn){

		String txt= "";
		try {
			txt = this.getPrintableTextRead(false, true, false);
		} catch (IOException e) {
			e.printStackTrace();
		} catch (InvalidGenomicCoordsException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (InvalidColourException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		StringBuilder sb= new StringBuilder();
		sb.append("On-screen text start: " + this.textStart + "\n");
		sb.append("On-screen text end: " + this.textEnd + "\n");
		sb.append("Text read: [" + txt + "]\n");				
		return sb.toString();

	}
	
	/*      S e t t e r s   and   G e t t e r s     */
	public int getTextStart() {
		return textStart;
	}
	public int getTextEnd() {
		return textEnd;
	}

	protected SAMRecord getSamRecord(){
		return this.samRecord;
	}
	
}
