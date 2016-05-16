package tracks;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.lang3.StringUtils;

import samTextViewer.GenomicCoords;
import samTextViewer.Utils;

public class TrackIntervalFeature extends Track {
 
	private List<IntervalFeature> intervalFeatureList= new ArrayList<IntervalFeature>();  
	private IntervalFeatureSet intervalFeatureSet;
	
	/* C o n s t r u c t o r */

	public TrackIntervalFeature(String filename, GenomicCoords gc) throws IOException{
		this.setGc(gc);
		this.setFilename(filename);
		this.intervalFeatureSet= new IntervalFeatureSet(filename);
		this.update();
	}
	
	/* Methods */
	
	public void update() throws IOException{
		this.intervalFeatureList = this.intervalFeatureSet.getFeaturesInInterval(
				this.getGc().getChrom(), this.getGc().getFrom(), this.getGc().getTo());
		for(IntervalFeature ift : intervalFeatureList){
			ift.mapToScreen(this.getGc().getMapping());
		}
	}
		
	@Override
	public void setHideRegex(String hideRegex) {
		this.intervalFeatureSet.setHideRegex(hideRegex);
	}
	@Override
	public String getHideRegex() {
		return this.intervalFeatureSet.getHideRegex();
	}

	@Override
	public void setShowRegex(String showRegex) {
		this.intervalFeatureSet.setShowRegex(showRegex);
	}
	@Override
	public String getShowRegex() {
		return this.intervalFeatureSet.getShowRegex();
	}
	
	@Override
	public String printToScreen() {
	
		List<String> printable= new ArrayList<String>();
		
		int nLines= 0;
		for(List<IntervalFeature> listToPrint : this.stackFeatures()){
			nLines++;
			if(nLines > this.yMaxLines){
				// Limit the number of lines in output
				break;
			}
			printable.add(this.printToScreenOneLine(listToPrint));
		}
		return StringUtils.join(printable, "\n");
	}
	
	/** Return a string of a single line of (typically de-stacked) reads
	 * */
	private String printToScreenOneLine(List<IntervalFeature> listToPrint) {
		
		List<String> printable= new ArrayList<String>();
		for(int i= 0; i < this.getGc().getMapping().size(); i++){ // First create empty track
			printable.add(" ");
		}
		for(IntervalFeature intervalFeature : listToPrint){
			if(intervalFeature.getScreenFrom() == -1){
				continue; // Feature doesn't map to screen, this shouldn't happen though
			}
			String x= intervalFeature.assignTextToFeature(this.isNoFormat());
			for(int j= intervalFeature.getScreenFrom(); j <= intervalFeature.getScreenTo(); j++){
				printable.set(j, x);
			}
		}
		return StringUtils.join(printable, "");
	}
	
	@Override
	public String getTitle(){
		return this.getFileTag() + "\n";
	}

	/**		
	 * Put in the same list reads that will go in the same line of text. 
	 * This method separates features touching or overlapping each other, useful for visualization.
	 * Each item of the output list is an IntervalFeatureSet going on its own line.
	 * 
	 * See also TrackReads.stackReads();
	 */
	private List<List<IntervalFeature>> stackFeatures(){
		
		// Make a copy of the IntervalFeature list. Items will be popped out as they are 
		// added to individual lines. 
		List<IntervalFeature> flatList= new ArrayList<IntervalFeature>();
		for(IntervalFeature x : this.intervalFeatureList){
			flatList.add(x);
		}
		
		List<List<IntervalFeature>> listOfLines= new ArrayList<List<IntervalFeature>>();
		if(flatList.size() == 0){
			return listOfLines;
		}
		List<IntervalFeature> line= new ArrayList<IntervalFeature>();
		line.add(flatList.get(0)); 
		flatList.remove(0);
		listOfLines.add(line);

		while(true){
			ArrayList<IntervalFeature> trToRemove= new ArrayList<IntervalFeature>();
			// Find a read in input whose start is greater then end of current
			for(int i=0; i < flatList.size(); i++){
				IntervalFeature intervalFeature= flatList.get(i);
				int gap= 1; // Add a space between book-end features
				if(intervalFeature.getScreenFrom() > line.get(line.size()-1).getScreenTo()+gap){ // +2 because we want some space between adjacent reads
					listOfLines.get(listOfLines.size()-1).add(intervalFeature); // Append to the last line. 
					trToRemove.add(intervalFeature);
				}
			} // At the end of the loop you have put in line as many reads as you can. 
			for(IntervalFeature intervalFeature : trToRemove){ 
				flatList.remove(flatList.indexOf(intervalFeature));
			}
			// Create a new line, add the first IntervalFeature in list
			if(flatList.size() > 0){
				line= new ArrayList<IntervalFeature>();
				line.add(flatList.get(0));
				listOfLines.add(line);
				flatList.remove(0);
			} else {
				break;
			}
		}
		return listOfLines;
	}

	@Override
	public String printFeatures(int windowSize){
		List<String> featureList= new ArrayList<String>();
		for(IntervalFeature ift : intervalFeatureList){
			featureList.add(ift.getRaw());
		}
		List<String> tabList= Utils.tabulateList(featureList);
		StringBuilder sb= new StringBuilder();
		for(String x : tabList){
			if(x.length() > windowSize){
				x= x.substring(0, windowSize);
			}			
			sb.append(x + "\n");
		}
		return sb.toString(); // NB: Leave last trailing /n
	}

	protected IntervalFeatureSet getIntervalFeatureSet() { return intervalFeatureSet; }

}
