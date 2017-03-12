package samTextViewer;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;

import org.apache.commons.io.FileUtils;

import exceptions.InvalidCommandLineException;
import exceptions.InvalidGenomicCoordsException;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;

public class GenomicCoordsHistory {

	private final String MARKER_FOR_POS= "## pos ##";
	private int countHistoricPositions= 0; // Number of valid positions read from file ~/.asciigenome_history and
	                                       // put in list of historic positions.
	
	/** List of genomic positions in order as they have been visited.*/
	private List<GenomicCoords> history= new ArrayList<GenomicCoords>();
	
	/** List index stating where we are in history */
	private int positionTracker= -1;

	/* Constructor */
	public GenomicCoordsHistory(){}
	
	/*  Methods */
	
	/** Add GenomicCoords obj to history provided this item is not equal in coordinates to 
	 * to the last one in history. 
	 * The position tracker is reset to the last when a new position is added */
	public void add(GenomicCoords gc) {
		if(this.history.size() == 0 || !this.history.get(this.history.size() - 1).equalCoords(gc)){
			this.history.add((GenomicCoords) gc.clone());
		}
		this.positionTracker= this.history.size() - 1;
	}

	public GenomicCoords current() {
		if(positionTracker < 0){
			return null;
		}
		return this.history.get(positionTracker);
	}

	public void previous() {

		int idx= this.positionTracker - 1;
		
		if(this.history.size() == 0){
			this.positionTracker= -1;
		} else if(this.history.size() == 1){
			this.positionTracker= 0;
		} else if (idx < 0){
			this.positionTracker= 0;
		} else {
			this.positionTracker= idx;
		}
		
	}

	public void next() {

		int idx= this.positionTracker + 1;
		
		if(this.history.size() == 0){
			this.positionTracker= -1;
		} else if(this.history.size() == 1){
			this.positionTracker= 0;
		} else if (idx >= this.history.size()){
			this.positionTracker= this.history.size() - 1;
		} else {
			this.positionTracker= idx;
		}
		
	}
	
	protected List<GenomicCoords> getHistory() {
		return history;
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
		
		for(GenomicCoords gc : this.getHistory()){ // Set the genome for each position in history. Invalid positions are removed.
			gc.setGenome(tokens);
			if( ! this.isValidPosition(gc, gc.getSamSeqDict())){
				toDelete.add(gc);
			} 
		}

		for(GenomicCoords gc : toDelete){
			this.getHistory().remove(gc);
		}

		this.positionTracker= this.getHistory().indexOf(now); // After having removed invalid positions, reset the position tracker to where we were.
		
		if(this.positionTracker < 0){ // The position from where the genome was set is not part of the dictionary 
			
			if(this.getHistory().size() > 0){ // Try to move to the last valid position.
				this.positionTracker= this.getHistory().size() - 1; 	
			
			} else { // There are no valid positions in the history so create a new one and move there

				GenomicCoords defaultPos= new GenomicCoords("default", null, null); 
				defaultPos.setGenome(tokens);
				
				LinkedHashMap<String, Integer> chromLen= new LinkedHashMap<String, Integer>();
				for(SAMSequenceRecord x : defaultPos.getSamSeqDict().getSequences()){
					chromLen.put(x.getSequenceName(), x.getSequenceLength());
				}
				String chrom= chromLen.keySet().iterator().next();
				GenomicCoords newGc= new GenomicCoords(chrom + ":" + 1, defaultPos.getSamSeqDict(), defaultPos.getFastaFile());
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
	public void readHistory(File historyFile, GenomicCoords checkGc) throws IOException {
		
		String[] hist = null;
		try{
			hist= FileUtils.readFileToString(historyFile).split("\\n");
		} catch(Exception e){
			System.err.println("Cannot read positions from file " + historyFile);
			return;
		}
		for(String line : hist){
			if(line.startsWith(MARKER_FOR_POS)){
				// try to create genomicCoords object using checkGc as template
				// If success, add this position to history list.
				try {
					String reg= line.replaceFirst(MARKER_FOR_POS, "");
					GenomicCoords gc= new GenomicCoords(reg, checkGc.getSamSeqDict(), checkGc.getFastaFile(), false);
					this.add(gc);
					this.countHistoricPositions += 1;
				} catch (Exception e){
					//
				}
			}
		}
	}

	/** Prepare a list of strings of positions ready to be written to 
	 * the ~/.asciigenome_history file.
	 * @throws IOException 
	 * */
	public List<String> prepareHistoryForHistoryFile(File historyFile, int maxPos) throws IOException {

		if(maxPos < 0){
			maxPos= 0;
		}
		
		List<String> hist= new ArrayList<String>();
		
		// First read all the positions from file, valid and non
		try{
			String[] histFile= FileUtils.readFileToString(historyFile).split("\\n");
			for(String line : histFile){
				if(line.startsWith(MARKER_FOR_POS)){
					hist.add(line);
				}
			}
		} catch(Exception e){
			System.err.println("Note: cannot read history file " + historyFile);
		}
		// Now add the positions visited in this session
		for(int i= this.countHistoricPositions; i < this.getHistory().size(); i++){
			hist.add(MARKER_FOR_POS + this.getHistory().get(i).toStringRegion());
		}

		// If necessary, trim the full list of positions to the max size
		if(hist.size() > maxPos){
			hist= hist.subList(hist.size() - maxPos, hist.size());
		}
		return hist;
	}
	
	/** Reset window size according to current terminal screen. 
	 * If the user reshapes the terminal window size or the font size, 
	 * detect the new size and add it to the history. 
	 * */
//	public void resetWindowSize() throws InvalidGenomicCoordsException, IOException{
//
//		int currentWindowSize= this.current().getUserWindowSize(); 
//		
//		int newSize= jline.TerminalFactory.get().getWidth() - 1;
//		if(newSize != this.current().getUserWindowSize()){
//			// Replace the current genomicCoords obj with a new one having the same coordinates but different windowSize.
//			// NB: The current genomic obj might not be the last one in the history list.
//			currentWindowSize= newSize;
//			String newRegion= this.current().getChrom() + ":" + this.current().getFrom() + "-" + this.current().getTo(); 
//			this.getHistory().add(
//					this.getHistory().indexOf(this.current()), 
//					new GenomicCoords(newRegion, this.current().getSamSeqDict(), currentWindowSize, this.current().getFastaFile())
//			);
//		}
//	}	
}
