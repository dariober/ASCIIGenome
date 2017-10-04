package samTextViewer;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import exceptions.InvalidCommandLineException;
import exceptions.InvalidGenomicCoordsException;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;

public class GenomicCoordsHistory {

	/** List of genomic positions in order as they have been visited.*/
	private List<GenomicCoords> currentSessionHistory= new ArrayList<GenomicCoords>();
	
	/** List index stating where we are in history */
	private int positionTracker= -1;

	/* Constructor */
	public GenomicCoordsHistory(){}
	
	/*  Methods */
	
	/** Add GenomicCoords obj to history provided this item is not equal in coordinates to 
	 * to the last one in history. 
	 * The position tracker is reset to the last when a new position is added */
	public void add(GenomicCoords gc) {
		if(this.currentSessionHistory.size() == 0 || !this.currentSessionHistory.get(this.currentSessionHistory.size() - 1).equalCoords(gc)){
			this.currentSessionHistory.add((GenomicCoords) gc.clone());
		}
		this.positionTracker= this.currentSessionHistory.size() - 1;
	}

	public GenomicCoords current() {
		if(positionTracker < 0){
			return null;
		}
		return this.currentSessionHistory.get(positionTracker);
	}

	public void previous() {

		int idx= this.positionTracker - 1;
		
		if(this.currentSessionHistory.size() == 0){
			this.positionTracker= -1;
		} else if(this.currentSessionHistory.size() == 1){
			this.positionTracker= 0;
		} else if (idx < 0){
			this.positionTracker= 0;
		} else {
			this.positionTracker= idx;
		}
		
	}

	public void next() {

		int idx= this.positionTracker + 1;
		
		if(this.currentSessionHistory.size() == 0){
			this.positionTracker= -1;
		} else if(this.currentSessionHistory.size() == 1){
			this.positionTracker= 0;
		} else if (idx >= this.currentSessionHistory.size()){
			this.positionTracker= this.currentSessionHistory.size() - 1;
		} else {
			this.positionTracker= idx;
		}
		
	}
	
	protected List<GenomicCoords> getCurrentSessionHistory() {
		return currentSessionHistory;
	}

	/** Set sequence dictionary and fasta ref, if available for all the genomic coords in this 
	 * history. 
	 * @throws InvalidGenomicCoordsException 
	 * */
	public void setGenome(List<String> tokens) throws InvalidCommandLineException, IOException, InvalidGenomicCoordsException {
		
		
		if(tokens.size() == 0){
			throw new InvalidCommandLineException();
		}

		List<GenomicCoords> toDelete= new ArrayList<GenomicCoords>();
		
		GenomicCoords now= this.current();
		for(GenomicCoords gc : this.getCurrentSessionHistory()){ // Set the genome for each position in history. Invalid positions are removed.
			gc.setGenome(tokens, true);
			if( ! this.isValidPosition(gc, gc.getSamSeqDict())){
				toDelete.add(gc);
			}
		}
		
		for(GenomicCoords gc : toDelete){
			this.getCurrentSessionHistory().remove(gc);
		}

		this.positionTracker= this.getCurrentSessionHistory().indexOf(now); // After having removed invalid positions, reset the position tracker to where we were.
		
		if(this.positionTracker < 0){ // The position from where the genome was set is not part of the dictionary 
			
			if(this.getCurrentSessionHistory().size() > 0){ // Try to move to the last valid position.
				this.positionTracker= this.getCurrentSessionHistory().size() - 1; 	
			
			} else { // There are no valid positions in the history so create a new one and move there

				int terminalWindowSize= Utils.getTerminalWidth();
				GenomicCoords defaultPos= new GenomicCoords("default", terminalWindowSize, null, null); 
				defaultPos.setGenome(tokens, true);
				
				LinkedHashMap<String, Integer> chromLen= new LinkedHashMap<String, Integer>();
				for(SAMSequenceRecord x : defaultPos.getSamSeqDict().getSequences()){
					chromLen.put(x.getSequenceName(), x.getSequenceLength());
				}
				String chrom= chromLen.keySet().iterator().next();
				GenomicCoords newGc= new GenomicCoords(chrom + ":" + 1, terminalWindowSize, defaultPos.getSamSeqDict(), defaultPos.getFastaFile());
				this.add(newGc);
				
			}
		}
	}

	/** Return true if the GenomicCoords object is valid given the SAMSequenceDictionary. I.e.
	 * gc contig and position is part of the dictionary.   
	 * */
	private boolean isValidPosition(GenomicCoords gc, SAMSequenceDictionary samSeqDict){

		if(samSeqDict == null || gc == null){
			return true;
		}
		
		HashMap<String, Integer> chromLen= new HashMap<String, Integer>();
		for(SAMSequenceRecord x : gc.getSamSeqDict().getSequences()){
			chromLen.put(x.getSequenceName(), x.getSequenceLength());
		}
		boolean isValid= false;
		if( chromLen.containsKey(gc.getChrom()) ){
			int len= chromLen.get(gc.getChrom());
			if(gc.getTo() <= len){
				isValid= true; 
			}
		} 
		return isValid;
		
	}

	/** Add to the GenomicCoordsHistory object the positions read from 
	 * historyFile.
	 * @throws IOException 
	 * */
	public void readHistory(String historyFileName, GenomicCoords checkGc) throws IOException {

		ASCIIGenomeHistory ag= new ASCIIGenomeHistory(historyFileName);
		int terminalWindowSize= Utils.getTerminalWidth();
		for(String reg : ag.getPositions()){
			try {
				// try to create genomicCoords object using checkGc as template
				// If success, add this position to history list.
				GenomicCoords gc= new GenomicCoords(reg, terminalWindowSize, checkGc.getSamSeqDict(), checkGc.getFastaFile(), false);
				this.add(gc);
			} catch (Exception e){
				//
			}
		}
	}

	/** Prepare a list of strings of positions ready to be written to 
	 * the ~/.asciigenome_history file.
	 * @throws IOException 
	 * */
	public List<String> prepareHistoryForHistoryFile(String historyFile, int maxPos) throws IOException {

		if(maxPos < 0){
			maxPos= 0;
		}

		List<String> positions= new ASCIIGenomeHistory(historyFile).getPositions(); 
		
		// Now add the positions visited in this session
		// We need to find the positions that were not inherited from
		// ASCIIGenomeHistory.
		List<String> newPositions= new ArrayList<String>();
		for(GenomicCoords gc : this.currentSessionHistory){
			newPositions.add(gc.toStringRegion());
		}
		
		for(String x : positions){
			// Remove from newPositions the positions taken from dot file
			if(newPositions.get(0).equals(x)){
				newPositions.remove(0);
			}
		}

		for(String x : newPositions){
			// Remove consecutive duplicates
			if( positions.size() == 0 || ! x.equals(positions.get(positions.size()-1))){
				positions.add(x);	
			}
		}

		// If necessary, trim the full list of positions to the max size
		if(positions.size() > maxPos){
			positions= positions.subList(positions.size() - maxPos, positions.size());
		}
		return positions;
	}
	
}