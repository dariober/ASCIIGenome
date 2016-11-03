package samTextViewer;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.sql.SQLException;
import coloring.Png;
import exceptions.InvalidCommandLineException;
import exceptions.InvalidGenomicCoordsException;
import exceptions.InvalidRecordException;
import tracks.Track;
import tracks.TrackReads;
import tracks.TrackSet;

/** Process a TrackSet given the necessary elements
 * */
public class TrackProcessor {
	
	private TrackSet trackSet;
	private boolean noFormat= false;
	private GenomicCoordsHistory genomicCoordsHistory; 
	private String snapshotFile= null;
	// int windowSize= 160;
	
	/* C O N S T R U C T O R S */
	
	public TrackProcessor(TrackSet trackSet, GenomicCoordsHistory genomicCoordsHistory) throws IOException, InvalidRecordException {
		this.trackSet= trackSet;
		this.genomicCoordsHistory= genomicCoordsHistory;
		
		/* Initialize GC profile */
		//if(genomicCoordsHistory.current().getFastaFile() != null){
		//	TrackWiggles gcProfile = genomicCoordsHistory.current().getGCProfile();
		//	this.gcProfileHashCode= gcProfile.hashCode();
		//	this.getTrackSet().add(gcProfile, "gcProfile");
		//}		
	}
	
	/* M E T H O D S */

	public void iterateTracks() throws IOException, InvalidGenomicCoordsException, InvalidRecordException, ClassNotFoundException, SQLException, InvalidCommandLineException{
		
		GenomicCoords currentGC= this.genomicCoordsHistory.current();
		
		if(currentGC.getChromIdeogram(20, this.noFormat) != null){
			Utils.printer(currentGC.getChromIdeogram(20, this.noFormat) + "\n", this.snapshotFile);
		}			

		for(Track track : trackSet.getTrackList()){
			
			track.setNoFormat(this.noFormat);
			if( ! track.getGc().equalCoordsAndWindowSize(currentGC) ){
				track.setGc(currentGC);
			}
			if(track.getyMaxLines() > 0 && !track.isHideTrack()){
				Utils.printer(track.getTitle(), this.snapshotFile);
				Utils.printer(track.printToScreen() + "\n", this.snapshotFile);
				Utils.printer(track.getPrintableConsensusSequence(), this.snapshotFile);
				Utils.printer(track.printFeatures(currentGC.getUserWindowSize()), this.snapshotFile);
			}
		}

		// Ruler and sequence
		// ------------------
		Utils.printer(currentGC.printableRefSeq(noFormat), snapshotFile);
		String ruler= currentGC.printableRuler(10, noFormat);
		Utils.printer(ruler + "\n", snapshotFile);

		// Position, memory, etc
		// ---------------------
		String footer= currentGC.toString() + "; " + 
				Math.rint(currentGC.getBpPerScreenColumn() * 10d)/10d + " bp/char; " + this.getMemoryStat();
		if(!noFormat){
			Utils.printer("\033[48;5;231;34m" + footer + "\033[48;5;231;30m; \n", this.snapshotFile);
		} else {
			Utils.printer(footer + "\n", this.snapshotFile);
		}

		// Optionally convert to png
		// -------------------------
		if(this.snapshotFile != null && this.snapshotFile.endsWith("png")){
			(new Png(new File(this.snapshotFile))).convert(new File(this.snapshotFile));
		} 	
	}
	
	public void exportTrackSetSettings(String filename) throws IOException{
		
		BufferedWriter wr= new BufferedWriter(new FileWriter(new File(filename)));
		if(this.getTrackSet().getTrackList().size() == 0){
			System.err.println("No track found. Nothing exported.");
			wr.close();
			return;
		}
		wr.write("goto " + this.getGenomicCoordsHistory().current().toStringRegion() + "\n");
		wr.write("setGenome " + this.getGenomicCoordsHistory().current().getSamSeqDictSource() + "\n");
		for(Track tr : this.getTrackSet().getTrackList()){
			if( ! (tr instanceof TrackReads) ){ // Track reads excluded as they are part of TrackCoverage. 
				                                // This is wrong as you could have one without the other!
				wr.write(tr.settingsToString() + "\n");
			}
		}
		wr.close();
	} 

	
	private String getMemoryStat(){
		float mem= (float) ((float)Runtime.getRuntime().totalMemory() / 1000000d);
		String memStats= "Mem: " +  Math.round(mem * 10)/10 + " MB";
		return memStats;
	}

	protected TrackSet getTrackSet() {
		return trackSet;
	}

	protected void setTrackSet(TrackSet trackSet) {
		this.trackSet = trackSet;
	}

	protected boolean isNoFormat() {
		return noFormat;
	}

	protected void setNoFormat(boolean noFormat) {
		this.noFormat = noFormat;
	}

	protected GenomicCoordsHistory getGenomicCoordsHistory() {
		return genomicCoordsHistory;
	}

	protected void setGenomicCoordsHistory(GenomicCoordsHistory genomicCoordsHistory) {
		this.genomicCoordsHistory = genomicCoordsHistory;
	}

	protected String getSnapshotFile() {
		return snapshotFile;
	}

	protected void setSnapshotFile(String snapshotFile) {
		this.snapshotFile = snapshotFile;
	}

	protected int getWindowSize() throws InvalidGenomicCoordsException, IOException {
		return this.genomicCoordsHistory.current().getUserWindowSize();
	}

}
