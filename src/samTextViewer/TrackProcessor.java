package samTextViewer;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.sql.SQLException;
import java.util.Iterator;

import org.apache.commons.lang3.StringUtils;

import com.google.common.base.Splitter;
import com.itextpdf.text.DocumentException;

import coloring.Config;
import coloring.ConfigKey;
import coloring.Pdf;
import exceptions.InvalidColourException;
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
	private boolean appendToSnapshotFile= false;
	// int windowSize= 160;
	
	/* C O N S T R U C T O R S */
	
	public TrackProcessor(TrackSet trackSet, GenomicCoordsHistory genomicCoordsHistory) throws IOException, InvalidRecordException {
		this.trackSet= trackSet;
		this.genomicCoordsHistory= genomicCoordsHistory;
	}
	
	/* M E T H O D S */

	public void iterateTracks() throws IOException, InvalidGenomicCoordsException, InvalidRecordException, ClassNotFoundException, SQLException, InvalidCommandLineException, DocumentException, InvalidColourException{
		
		GenomicCoords currentGC= this.genomicCoordsHistory.current();
		
		StringBuilder outputString= new StringBuilder();
		
		if(currentGC.getChromIdeogram(20, this.noFormat) != null){
			outputString.append(currentGC.getChromIdeogram(20, this.noFormat) + "\n");
		}			

		// Update tracks to new genomic coords
		for(Track track : trackSet.getTrackList()){
			if( ! track.getGc().equalCoordsAndWindowSize(currentGC) ){
				track.setGc(currentGC);
			}			
		}
		// Set new y limits as required. This step has to come after the positioning to new coordinates because
		// we may need to autoscale to global min or max.
		this.getTrackSet().setAutoYLimits();

		// Visualize as required
		for(Track track : trackSet.getTrackList()){
			
			track.setNoFormat(this.noFormat);
			if(track.getyMaxLines() > 0 && !track.isHideTrack()){
				outputString.append(track.getTitle());
				outputString.append(track.printToScreen() + "\n");
				outputString.append(track.getPrintableConsensusSequence());
				outputString.append(track.printFeaturesToFile());
			}
		}

		// Ruler and sequence
		// ------------------
		outputString.append(currentGC.printableRefSeq(noFormat)); // this.printer(currentGC.printableRefSeq(noFormat), snapshotFile);
		outputString.append(currentGC.printableRuler(10, noFormat) + "\n"); // this.printer(ruler + "\n", snapshotFile);

		// Position, memory, etc
		// ---------------------
		String footer= this.getFooter(currentGC);
		if(!noFormat){
			outputString.append("\033[48;5;");
			outputString.append(Config.getColor(ConfigKey.background));
			outputString.append(";38;5;");
			outputString.append(Config.getColor(ConfigKey.footer));
			outputString.append("m");
			outputString.append(footer);
			outputString.append("\033[38;5;");
			outputString.append(Config.getColor(ConfigKey.foreground));
			outputString.append("m");
		} else {
			outputString.append(footer);
		}

		String printable= this.fillUpLines(outputString.toString());
		
		// Print to screen
		System.out.println(printable);
		
		// Optionally save to file
		// -----------------------
		if(this.snapshotFile != null && this.snapshotFile.endsWith(".pdf")){
			(new Pdf(printable)).convert(new File(this.snapshotFile), 10, this.appendToSnapshotFile);
		
		} else if(this.snapshotFile != null){
			BufferedWriter wr= new BufferedWriter(new FileWriter(new File(this.snapshotFile), this.appendToSnapshotFile));
			wr.write(Utils.stripAnsiCodes(printable));
			wr.write("\n-------8<-------------[ cut here ]----------------------\n\n");
			wr.close();
		}
		this.snapshotFile= null;
	}
	
	/** String screenshot is the string almost ready to be printed to screen or file.
	 * What is missing is to fill up lines with whitespaces until the end of the screen so
	 * that terminals like tmux do not show lines of mixed colours. 
	 * @throws IOException 
	 * @throws InvalidGenomicCoordsException 
	 * */
	private String fillUpLines(String screenshot) throws InvalidGenomicCoordsException, IOException{
		StringBuilder full= new StringBuilder();
		Iterator<String> iter = Splitter.on("\n").split(screenshot).iterator();
		int screenSize= this.getWindowSize();
		while(iter.hasNext()){
			String line= iter.next();
			int npad= screenSize - Utils.stripAnsiCodes(line).length();
// This check throws false positives when `print -full` finds lines longer than screenwidth
//			if(npad < 0){
//				throw new RuntimeException("Line width greater then screenwidth! This should not happen!\n"
//						+ "Please report this as a bug.");
//			}
			if(npad < 0){
				npad= 0; // Nothing to be done.
			}
			String filler= StringUtils.repeat(' ', npad);
			full.append(line);
			full.append(filler);
			if(iter.hasNext()){
				full.append("\n");
			}
		}
		return full.toString();
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
		
		filename= Utils.tildeToHomeDir(filename);
		
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

//	private void printer(String xprint, String filename) throws IOException{
//		System.out.print(xprint);
//		if(filename == null){
//			return;
//		}
//		if(! filename.toLowerCase().endsWith(".pdf")){
//			// We write file as plain text so strip ansi codes.
//			xprint= Utils.stripAnsiCodes(xprint);
//		}
//		BufferedWriter wr= new BufferedWriter(new FileWriter(new File(filename), true));
//		wr.write(xprint);
//		wr.close();
//	}

	
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

	protected boolean isAppendToSnapshotFile() {
		return appendToSnapshotFile;
	}

	protected void setAppendToSnapshotFile(boolean appendToSnapshotFile) {
		this.appendToSnapshotFile = appendToSnapshotFile;
	}

}
