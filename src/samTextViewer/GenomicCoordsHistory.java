package samTextViewer;

import java.util.ArrayList;
import java.util.List;

public class GenomicCoordsHistory {

	/** List of genomic positions in order as they have been visited.*/
	private List<GenomicCoords> history= new ArrayList<GenomicCoords>();
	
	/** List index stating where we are in history */
	private int positionTracker= -1;

//	private String seqRegex= ""; 
	
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

//	public TrackIntervalFeature findRegex() throws ClassNotFoundException, IOException, InvalidGenomicCoordsException, InvalidRecordException, SQLException {
//
//		int prevHeight;
//		TrackIntervalFeature prevRegexMatchTrack= null;
//		if(this.getHistory().size() >= 1){
//			prevRegexMatchTrack = this.getHistory().get(this.getHistory().size()-1).getRegexMatchTrack();
//		} 
//		if(prevRegexMatchTrack == null){
//			 prevHeight= 10;
//		} else {
//			prevHeight= prevRegexMatchTrack.getyMaxLines();
//		}
//		TrackIntervalFeature seqRegexTrack = this.current().findRegex(this.seqRegex);
//		seqRegexTrack.setyMaxLines(prevHeight);
//
//		// Track is created, delete tmp files.
//		new File(seqRegexTrack.getFilename()).delete();
//		new File(seqRegexTrack.getFilename() + ".tbi").delete();
//		
//		return seqRegexTrack;
//	}

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
