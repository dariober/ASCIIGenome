package tracks;

import java.util.ArrayList;
import java.util.List;
import htsjdk.samtools.filter.SamRecordFilter;
import samTextViewer.GenomicCoords;

public class Track {

	private String title= "";
	protected int yMaxLines= 10;
	private String filename= "N/A"; // File name as given in input
	private String fileTag= "N/A"; // File name for title
	private List<Double> screenScores= new ArrayList<Double>();
	private GenomicCoords gc;
	private boolean noFormat= false; 
	private double yLimitMin= Double.NaN; // Same as R ylim()
	private double yLimitMax= Double.NaN;
	private List<SamRecordFilter> filters= new ArrayList<SamRecordFilter>();
	private boolean bs; 
	/** Max size of genomic region before the track shuts down to prevent excessive slow down */
	protected final int MAX_REGION_SIZE= 100000;   
	
//	public Track(){}

	/* Min value of screen scores. Not to be confused with the y limit **/
	public double getMinScreenScores(){
		Double ymin= Double.NaN;
	 	for(Double x : this.screenScores){
	 		if(ymin.isNaN() && !x.isNaN()){
	 			ymin= x;
	 		} else if (x < ymin){
	 			ymin= x;
	 		}
	 	}
	 	return ymin;
	}
	
	/* Max value of screen scores. Not to be confused with the y limit **/
	public double getMaxScreenScores(){
		Double ymax= Double.NaN;
	 	for(Double x : this.screenScores){
	 		if(ymax.isNaN() && !x.isNaN()){
	 			ymax= x;
	 		} else if (x > ymax){
	 			ymax= x;
	 		}
	 	}
	 	return ymax;
	}
	
	/* Printers */
	public String printToScreen(){
		return null;
	}

	public String printFeatures(int windowSize){
		return "";
	};
	
	public String toString(){
		return this.getFilename();
	}
	
	/* Setters and getters */
	public void setTitle(String title){
		this.title= title;
	}
	public String getTitle(){
		return this.title;
	}
	public int getyMaxLines() {
		return yMaxLines;
	}
	public void setyMaxLines(int yMaxLines) {
		this.yMaxLines = yMaxLines;
	}
	public String getFilename() {
		return filename;
	}
	public void setFilename(String filename) {
		this.filename = filename;
	}

	public String getFileTag() { return fileTag; }
	public void setFileTag(String fileTag) { this.fileTag = fileTag; }
	
	protected List<Double> getScreenScores() {
		return screenScores;
	}
	protected void setScreenScores(List<Double> screenScores) {
		this.screenScores = screenScores;
	}

	public GenomicCoords getGc() {
		return gc;
	}

	public void setGc(GenomicCoords gc) {
		this.gc = gc;
	}

	public boolean isNoFormat() { return noFormat; }
	public void setNoFormat(boolean noFormat) { this.noFormat = noFormat; }

	public double getYLimitMin() { return yLimitMin;}
	public void setYLimitMin(double ymin) { this.yLimitMin = ymin;}

	public double getYLimitMax() { return yLimitMax; }
	public void setYLimitMax(double ymax) { this.yLimitMax = ymax; }

	public List<SamRecordFilter> getFilters() { return filters; }
	public void setFilters(List<SamRecordFilter> filters) { this.filters = filters; }

	public boolean isBs() { return bs; }
	public void setBs(boolean bs) { this.bs = bs; }

	public void setHideRegex(String hideRegex) { }
	public String getHideRegex() { return ""; }
	
	public void setShowRegex(String showRegex) { }
	public String getShowRegex() { return ""; }
	
}

