package tracks;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import com.google.common.base.Joiner;

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
class TextRead {
	
	// Characters for methylation coding. NB: M and U are valid DNA chars in IUPAC!
	private final static char charM= 'M';
	private final static char charU= 'U';
	private final static char charm= 'm';
	private final static char charu= 'u';
	private final static char charFwd= '>';
	private final static char charRev= '<';
	private final static int  SHADE_MAPQ= 5;
	
	/** Char to represent deletions from the reference. I.e. gaps in the read */
	final private char DEL= '-';
	/** Char to represent region skip. I.e. gaps in the read */
	final private char N= '_';
	
	/** Start position of the read in window coordinates. 1-based. If in genomic
	 * coords read starts aligning at pos 100 and window span is 100:150, then textStart= 1.*/
	private int textStart;

	/** End position of the read in window coordinates. 1-based. If in genomic
	 * coords read ends pos 150 and window span is 100:150, then textEnd= 1. */
	private int textEnd;
	
	private SAMRecord rec;
	private GenomicCoords gc;
	
	/*    C o n s t r u c t o r s    */
		
	public TextRead(SAMRecord rec, GenomicCoords gc) throws InvalidGenomicCoordsException, IOException{
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
		this.rec= rec;
		this.setTextStart();
		this.setTextEnd();
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
		List<Character> unformatted;
		if(!bs){
			unformatted= this.getConsRead();
		} else {
			unformatted= this.convertDnaReadToTextReadBS();
		}
		if(withReadName){ // Replace bases with read name. As long as it fits
			int upTo= (unformatted.size() < this.rec.getReadName().length()) 
					? unformatted.size() 
					: this.rec.getReadName().length();  
			for(int i= 0; i< upTo; i++){
				char nchar= (char) this.rec.getReadName().getBytes()[i];
				unformatted.set(i, nchar);
			}
			if(unformatted.size() > upTo){ // Add a separator btw read name and sequence
				unformatted.set(upTo, '/');
			}
			return Joiner.on("").join(unformatted);
		} else {
			return readFormatter(unformatted, noFormat, bs);
		}
	}
	
	/**
	 * Return the read chars in read as nicely formatted string. This method only adds non-ascii formatting
	 * to the individual bases (chars). The bases themselves and their capitalization are not touched.
	 * Formatting depends on the base (actg, m/u) and on the read it comes from (strand & mate).
	 * 
	 * For formatting see http://misc.flogisoft.com/bash/tip_colors_and_formatting
	 * and http://stackoverflow.com/questions/5947742/how-to-change-the-output-color-of-echo-in-linux
	 * 
	 * @param read
	 * @param noFormat
	 * @return
	 * @throws IOException 
	 * @throws InvalidGenomicCoordsException 
	 * @throws InvalidColourException 
	 */
	private String readFormatter(List<Character> read, boolean noFormat, boolean bs) throws InvalidGenomicCoordsException, IOException, InvalidColourException{

		if(noFormat){ // Essentially nothing to do in this case
			return Joiner.on("").join(read);
		}
		StringBuilder formatted= new StringBuilder();		
		for(char c : read){ // Each base is formatted independently from the others
			formatted.append("\033["); // Start format 
			if(this.rec.getReadPairedFlag() && this.rec.getSecondOfPairFlag()){
				formatted.append("4;"); // Underline 2nd in pair
			}					
			if(this.rec.getMappingQuality() < SHADE_MAPQ){ // Grey out low mapq
				formatted.append("48;5;"); 
				formatted.append(Config.getColor(ConfigKey.shade_low_mapq)); 
				formatted.append(";38;5;");
				formatted.append(Config.getColor(ConfigKey.foreground));
			} else if(Character.toUpperCase(c) == charM){
				formatted.append("48;5;"); 
				formatted.append(Config.getColor(ConfigKey.methylated_background)); 
				formatted.append(";38;5;");
				formatted.append(Config.getColor(ConfigKey.methylated_foreground));
			} else if(Character.toUpperCase(c) == charU){
				formatted.append("48;5;"); 
				formatted.append(Config.getColor(ConfigKey.unmethylated_background)); 
				formatted.append(";38;5;");
				formatted.append(Config.getColor(ConfigKey.unmethylated_foreground));
			} else if(Character.toUpperCase(c) == 'A'){
				formatted.append("1;38;5;"); 
				formatted.append(Config.getColor(ConfigKey.seq_a));
			} else if(Character.toUpperCase(c) == 'C') {
				formatted.append("1;38;5;"); 
				formatted.append(Config.getColor(ConfigKey.seq_c));
			} else if(Character.toUpperCase(c) == 'G') {
				formatted.append("1;38;5;"); 
				formatted.append(Config.getColor(ConfigKey.seq_g));
			} else if(Character.toUpperCase(c) == 'T') {
				formatted.append("1;38;5;"); 
				formatted.append(Config.getColor(ConfigKey.seq_t));
			} else if(!this.rec.getReadNegativeStrandFlag() && !(bs && !(gc.getBpPerScreenColumn() > 1))){
				formatted.append("48;5;"); 
				formatted.append(Config.getColor(ConfigKey.feature_background_positive_strand)); 
				formatted.append(";38;5;");
				formatted.append(Config.getColor(ConfigKey.foreground));
			} else if(this.rec.getReadNegativeStrandFlag() && !(bs && !(gc.getBpPerScreenColumn() > 1))){
				formatted.append("48;5;"); 
				formatted.append(Config.getColor(ConfigKey.feature_background_negative_strand)); 
				formatted.append(";38;5;");
				formatted.append(Config.getColor(ConfigKey.foreground));
			} else {
				formatted.append("48;5;"); 
				formatted.append(Config.getColor(ConfigKey.background)); 
				formatted.append(";38;5;");
				formatted.append(Config.getColor(ConfigKey.foreground));
			}
		}
		return formatted.toString();
	}
	
	/**
	 * Obtain the start position of the read on screen. Screen positions are 1-based. 
	 * @return
	 * @throws IOException 
	 * @throws InvalidGenomicCoordsException 
	 */
	private void setTextStart() throws InvalidGenomicCoordsException, IOException{		
		if(rec.getAlignmentStart() <= gc.getFrom()){ // Read starts right at the window start or even earlier
			this.textStart= 1;
			return;
		}		
		this.textStart= Utils.getIndexOfclosestValue(rec.getAlignmentStart(), gc.getMapping(gc.getUserWindowSize())) + 1;
		return;
	}
	
	private void setTextEnd() throws InvalidGenomicCoordsException, IOException{
		//if(rec.getAlignmentEnd() >= gc.getTo()){
		//	this.textEnd= gc.getUserWindowSize() < gc.getGenomicWindowSize() ?  
		//			gc.getUserWindowSize() : gc.getGenomicWindowSize();
		//	return;
		//}
		this.textEnd= Utils.getIndexOfclosestValue(rec.getAlignmentEnd(), gc.getMapping(gc.getUserWindowSize())) + 1;
		return;
	}
	
	/** If the windowSize and genomic span are not mapped 1:1, i.e. 1 bp : 1 char, then
	 * represent reads as simplified bases */
	private List<Character> getSquashedRead(){
				
		ArrayList<Character> squashedRead= new ArrayList<Character>();
		char xc;
		if(this.rec.getReadNegativeStrandFlag()){
			xc= charRev;
		} else {
			xc= charFwd;
		}
		for(int i= this.textStart; i <= this.textEnd; i++){
			squashedRead.add(xc);
		}
		return squashedRead;
	}
	
	/** Get a representation of the read as it appears aligned to the reference. 
	 * I.e. clipped ends omitted and deletions appearing as gaps (empty byte).
	 * Only the portion contained between the genomic coords from:to is returned.
	 * @return
	 * @throws IOException 
	 * @throws InvalidGenomicCoordsException 
	 */
	private List<Character> getDnaRead() throws InvalidGenomicCoordsException, IOException {

		if(this.gc.getBpPerScreenColumn() > 1){
			return this.getSquashedRead();
		}
		
		// Accumulate here the read bases inside the window 
		ArrayList<Character> dnaRead= new ArrayList<Character>();
		byte[] readBases= rec.getReadBases();
		// Walk along the aligned read and append bases to textRead as long as
		// the genomic position of the base is inside the genomic coords of the window
		int curBaseGenomicPos= rec.getAlignmentStart();
		int curBaseReadPos= 0; // Position on read. Start from zero walk along the read
		List<CigarElement> cigarEls= rec.getCigar().getCigarElements();
		for(CigarElement el : cigarEls){
			if(el.getOperator() == CigarOperator.M || el.getOperator() == CigarOperator.EQ || el.getOperator() == CigarOperator.X){
				for(int i= 0; i < el.getLength(); i++){
					if(curBaseGenomicPos >= gc.getFrom() && curBaseGenomicPos <= gc.getTo()){
						// If base is inside window:
						if(readBases.length > 0){
							dnaRead.add((char)readBases[curBaseReadPos]);
						} else { // If sam record has no read seq stored put N
							dnaRead.add('N');
						}
					}
					curBaseGenomicPos++; // M consumes read and ref bases. So increment them
					curBaseReadPos++;
				}
			} else if(el.getOperator() == CigarOperator.D || el.getOperator() == CigarOperator.N){
				for(int i= 0; i < el.getLength(); i++){
					if(curBaseGenomicPos >= gc.getFrom() && curBaseGenomicPos <= gc.getTo()){
						if(el.getOperator() == CigarOperator.D){
							dnaRead.add(DEL);
						} else if(el.getOperator() == CigarOperator.N){ 
							dnaRead.add(N);
						} else {
							System.err.println("Unexpected operator");
							throw new RuntimeException();
						}
					}
					curBaseGenomicPos++;
				}
			} else if(el.getOperator() == CigarOperator.I) {
				curBaseReadPos += el.getLength(); // Insertions in the reference are missed
			} else if(el.getOperator() == CigarOperator.S){
				curBaseReadPos += el.getLength();
			} else if(el.getOperator() == CigarOperator.H){
				// Nothing to do
			} else if(el.getOperator() == CigarOperator.P){
				// Nothing to do: NOT SURE is is correct to just ignore padding!
			} else {
				System.err.println("Unexpected operator in cigar string for record\n" + rec.getSAMString()); 
				throw new RuntimeException();
			}
		}
		for(int i= 0; i < dnaRead.size(); i++){
			if(this.rec.getReadNegativeStrandFlag()){
				dnaRead.set(i, Character.toLowerCase(dnaRead.get(i)));
			} else {
				dnaRead.set(i, Character.toUpperCase(dnaRead.get(i)));
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
	private List<Character> getConsRead() throws IOException, InvalidGenomicCoordsException {
		
		List<Character> dnaRead= this.getDnaRead();
		if(this.gc.getRefSeq() == null){
			return dnaRead;
		}
		List<Character> consRead= new ArrayList<Character>();
		int posOnRead= 0;
		for(int i= this.textStart - 1; i < this.textEnd; i++){
			char base= Character.toUpperCase( dnaRead.get(posOnRead) );
			char ref= (char) Character.toUpperCase(this.gc.getRefSeq()[i]);
			if( base == ref){
				if(this.rec.getReadNegativeStrandFlag()){
					consRead.add(',');
				} else {
					consRead.add('.');
				}
			} else {
				if(this.rec.getReadNegativeStrandFlag()){
					consRead.add(Character.toLowerCase(base));
				} else {
					consRead.add(base);
				}				
			}
			posOnRead++;
		}
		return consRead;
	}

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
	private List<Character> convertDnaReadToTextReadBS() throws IOException, InvalidGenomicCoordsException{
	
		if(this.gc.getRefSeq() == null){ // Effectively don't convert 
			return this.getConsRead();
		}
		
		// for(int i=0; i < this.gc.getRefSeq().length; i++){ // Make ref uppercase
		// 	this.gc.getRefSeq()[i]= (byte) Character.toUpperCase(this.gc.getRefSeq()[i]);
		// }
		
		// For convenience extract flags from sam record
		boolean isSecondOfPair= false;
		if(rec.getReadPairedFlag() && rec.getSecondOfPairFlag()){
			isSecondOfPair= true;
		}
		boolean isForwardStrand= !this.rec.getReadNegativeStrandFlag();
		
		List<Character> textReadBS= this.getConsRead(); // Iterate through each base to set methyl state
		for(int i= 0; i < textReadBS.size(); i++){
			char ref= (char) this.gc.getRefSeq()[i + this.textStart - 1];
			ref= Character.toUpperCase(ref);
			char read= textReadBS.get(i);
			if( ( isForwardStrand && !isSecondOfPair ) || ( !isForwardStrand && isSecondOfPair )){
				// Look for C on the reference
				if(isForwardStrand){ // +strand, first in pair or unpaired
					if(ref == 'C' && read == '.'){
						textReadBS.set(i, charM);
					} else if(ref == 'C' && read == 'T'){
						textReadBS.set(i, charU);
					} else {
						// Nothing to change
					}
				} else { // -ve strand, 2nd in pair
					// Look for c=',' -> m; 't' -> u
					if(ref == 'C' && read == ','){
						textReadBS.set(i, charm);
					} else if(ref == 'C' && read == 't'){
						textReadBS.set(i, charu);
					} else {
						// Nothing to change
					}
				}
			} else if( ( !isForwardStrand && !isSecondOfPair ) || ( isForwardStrand && isSecondOfPair )){
				// Look for G on the reference
				if(!isForwardStrand){ // -ve strand; first in pair or unpaired
					if(ref == 'G' && read == ','){
						textReadBS.set(i, 'm');
					} else if(ref == 'G' && read == 'a'){
						textReadBS.set(i, 'u');
					} else {
						// Nothing to change
					}
				} else { // -ve strand, 2nd in pair
					if(ref == 'G' && read == '.'){
						textReadBS.set(i, 'M');
					} else if(ref == 'G' && read == 'A'){
						textReadBS.set(i, 'U');
					} else {
						// Nothing to change
					}
				}
			}
		}
		return textReadBS;
	}
	
	public String toString(){

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

}
