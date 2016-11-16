package samTextViewer;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.sql.SQLException;

import org.apache.commons.lang3.StringUtils;

import com.itextpdf.text.DocumentException;

import coloring.Pdf;
import exceptions.InvalidCommandLineException;
import exceptions.InvalidGenomicCoordsException;
import exceptions.InvalidRecordException;
import tracks.Track;
import tracks.TrackBookmark;
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

	public void iterateTracks() throws IOException, InvalidGenomicCoordsException, InvalidRecordException, ClassNotFoundException, SQLException, InvalidCommandLineException, DocumentException{
		
		GenomicCoords currentGC= this.genomicCoordsHistory.current();
		
		if(currentGC.getChromIdeogram(20, this.noFormat) != null){
			Utils.printer(currentGC.getChromIdeogram(20, this.noFormat) + "\n", this.snapshotFile);
		}			

		this.getTrackSet().setAutoYLimits();
		
		for(Track track : trackSet.getTrackList()){
			
			track.setNoFormat(this.noFormat);
			if( ! track.getGc().equalCoordsAndWindowSize(currentGC) ){
				track.setGc(currentGC);
			}
			if(track.getyMaxLines() > 0 && !track.isHideTrack()){
				Utils.printer(track.getTitle(), this.snapshotFile);
				Utils.printer(track.printToScreen() + "\n", this.snapshotFile);
				Utils.printer(track.getPrintableConsensusSequence(), this.snapshotFile);
				Utils.printer(track.printFeaturesToFile(), this.snapshotFile);
			}
		}

		// Ruler and sequence
		// ------------------
		Utils.printer(currentGC.printableRefSeq(noFormat), snapshotFile);
		String ruler= currentGC.printableRuler(10, noFormat);
		Utils.printer(ruler + "\n", snapshotFile);

		// Position, memory, etc
		// ---------------------
		String footer= this.getFooter(currentGC);
		if(!noFormat){
			Utils.printer("\033[48;5;231;34m" + footer + "\033[48;5;231;30m\n", this.snapshotFile);
		} else {
			Utils.printer(footer + "\n", this.snapshotFile);
		}

		// Optionally convert to pdf
		// -------------------------
		if(this.snapshotFile != null && this.snapshotFile.endsWith(".pdf")){
			(new Pdf(new File(this.snapshotFile))).convert(new File(this.snapshotFile), 10);
		} 	
	}
	
	private String getFooter(GenomicCoords currentGC) throws InvalidGenomicCoordsException, IOException {
		String footer= currentGC.toString() + "; " + 
				Math.rint(currentGC.getBpPerScreenColumn() * 10d)/10d + " bp/char; " + this.getMemoryStat();
        // Add midpoint marker
        int mid= (int) Math.rint(this.getWindowSize() / 2.0);
        int npad= mid - footer.length();
        if(npad > 1){
            String pad= StringUtils.repeat(" ", npad - 1);
            footer += pad;
            footer += "/\\";
        }
		 return footer;
	}

	public void exportTrackSetSettings(String filename) throws IOException{
		
		BufferedWriter wr= new BufferedWriter(new FileWriter(new File(filename)));
		if(this.getTrackSet().getTrackList().size() == 0){
			System.err.println("No track found. Nothing exported.");
			wr.close();
			return;
		}
		// First we write the bookmark track, if any. Doing it first avoids moving around heavy files.
		for(Track tr : this.getTrackSet().getTrackList()){
			if( tr instanceof TrackBookmark ){
				wr.write(tr.settingsToString() + "\n");
			}
		}
		
		wr.write("goto " + this.getGenomicCoordsHistory().current().toStringRegion() + "\n");
		wr.write("setGenome " + this.getGenomicCoordsHistory().current().getSamSeqDictSource() + "\n");
		String orderTracks= "orderTracks ";
		for(Track tr : this.getTrackSet().getTrackList()){
			orderTracks += tr.getTrackTag().replaceAll("(#|@)\\d+$", "") + " ";
			if(tr instanceof TrackBookmark){
				continue; // Already done.
			}
			if( ! (tr instanceof TrackReads) ){ // Track reads excluded as they are part of TrackCoverage. 
				                                // This is wrong as you could have one without the other!
				wr.write(tr.settingsToString() + "\n");
			}
		}
		wr.write(orderTracks + "\n"); // Note that re-ordering may be different from original as the ID part was stripped! 
		wr.close();
	} 
	
	private String getMemoryStat() throws InvalidGenomicCoordsException, IOException{
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
