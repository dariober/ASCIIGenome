package tracks;

import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.lang3.StringUtils;

import exceptions.InvalidGenomicCoordsException;
import exceptions.InvalidRecordException;
import samTextViewer.GenomicCoords;
import samTextViewer.Utils;

public class TrackIntervalFeature extends Track {
 
	private List<IntervalFeature> intervalFeatureList= new ArrayList<IntervalFeature>();  
	/**For GTF/GFF data: Use this attribute to get the feature names 
	 * */
	protected IntervalFeatureSet intervalFeatureSet;
	
	/* C o n s t r u c t o r */

	public TrackIntervalFeature(String filename, GenomicCoords gc) throws IOException, InvalidGenomicCoordsException, ClassNotFoundException, InvalidRecordException, SQLException{
		this.setGc(gc);
		this.setFilename(filename);
		this.intervalFeatureSet= new IntervalFeatureSet(filename);
		this.update();
	}
	
	public TrackIntervalFeature(IntervalFeatureSet intervalFeatureSet, GenomicCoords gc) throws IOException, InvalidGenomicCoordsException{
		this.setGc(gc);
		this.intervalFeatureSet= intervalFeatureSet;
		this.update();
	}
	
	/* M e t h o d s */
	
	public void update() throws IOException, InvalidGenomicCoordsException{
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
	public String printToScreen() throws InvalidGenomicCoordsException {
	
		List<String> printable= new ArrayList<String>();		
		int nLines= 0;
		try {
			for(List<IntervalFeature> listToPrint : this.stackFeatures()){
				nLines++;
				if(nLines > this.yMaxLines){
					// Limit the number of lines in output
					break;
				}
				printable.add(this.printToScreenOneLine(listToPrint));
			}
		} catch (Exception e) {
			e.printStackTrace();
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
			intervalFeature.setGtfAttributeForName(this.getGtfAttributeForName());
			String nameOnFeature= intervalFeature.getName().trim() + "_";
			int relPos= 0;
			for(int j= intervalFeature.getScreenFrom(); j <= intervalFeature.getScreenTo(); j++){
				
				// Default is to use the feature type as printable text, e.g. 'E', 'T', '>', '<', '|', etc.
				String text= intervalFeature.assignTextToFeature(this.isNoFormat()); 
 
				if((intervalFeature.getScreenTo() - intervalFeature.getScreenFrom() + 1) > 4
						// (j - intervalFeature.getScreenFrom()) > 0 // First char is feature type (E, C, etc.)
						&& j < intervalFeature.getScreenTo() // Last char is feature type (E, C, etc.)
						&& relPos < nameOnFeature.length() 
						&& !nameOnFeature.equals("._")){
					// If these conds are satisfied, use the chars in the name as printable chars. 
					Character x= nameOnFeature.charAt(relPos); 
					if(this.isNoFormat()){
						text= Character.toString(x); 	
					} else {
						text= FormatGTF.format(x, intervalFeature.getStrand());
					}
					relPos += 1;
				}

				printable.set(j, text);
			}
		}
		return StringUtils.join(printable, "");
	}
	
	@Override
	public String getTitle(){
		String sq= "";
		if(this.getFeatureDisplayMode().equals(FeatureDisplayMode.SQUASHED)){
			sq= "; squashed";
		} else if (this.getFeatureDisplayMode().equals(FeatureDisplayMode.MERGED)){
			sq= "; merged";
		}
		String gapped= "";
		if(this.getGap() == 0){
			gapped= "; ungapped";
		}
		String title=  this.getTrackTag() + "; " 
	                 + "Incl " + this.getShowRegex()
	                 + " Excl " + this.getHideRegex()
	                 + " N: " + this.intervalFeatureList.size()
	                 + sq
	                 + gapped;
		this.getHideRegex();
		return this.formatTitle(title) + "\n";
	}

	
	
	/** Remove positional duplicates from list of interval features for more compact visualization. 
	 * Squashing is done according to feature field which should be applicable to GTF/GFF only.*/
	private List<IntervalFeature> squashFeatures(List<IntervalFeature> intervalList){

		List<IntervalFeature> stack= new ArrayList<IntervalFeature>();
		List<IntervalFeature> squashed= new ArrayList<IntervalFeature>();
		for(IntervalFeature interval : intervalList){
			if(stack.size() == 0 || stack.get(0).equalStranded(interval)){
				// Accumulate features with same coords.
				stack.add(interval);
			} else {
				squashed.add(stack.get(0));
				stack.clear();
				stack.add(interval);
			}
		}
		if(stack.size() > 0){
			squashed.add(stack.get(0));
		}
		return squashed;
	}
	
	/**		
	 * Put in the same list reads that will go in the same line of text. 
	 * This method separates features touching or overlapping each other, useful for visualization.
	 * Each item of the output list is an IntervalFeatureSet going on its own line.
	 * @param space Space between text feature that touch each other on screen. Use 1 to have at least one space
	 * so that distinct features never look merged with adjacent ones. 0 will not put any space and will give more
	 * compact view.  
	 * See also TrackReads.stackReads();
	 * @throws InvalidGenomicCoordsException 
	 */
	private List<List<IntervalFeature>> stackFeatures() throws InvalidGenomicCoordsException{
		
		List<IntervalFeature> intervals; 
		if(this.getFeatureDisplayMode().equals(FeatureDisplayMode.SQUASHED)){
			intervals = this.squashFeatures(this.intervalFeatureList);
		} else if(this.getFeatureDisplayMode().equals(FeatureDisplayMode.MERGED)) {
			intervals = Utils.mergeIntervalFeatures(this.intervalFeatureList);
		} else {
			intervals = this.intervalFeatureList;
		}
		
		// Make a copy of the IntervalFeature list. Items will be popped out as they are 
		// added to individual lines. 
		List<IntervalFeature> flatList= new ArrayList<IntervalFeature>();
		for(IntervalFeature x : intervals){ // this.intervalFeatureList
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
				// int gap= 1; // Add a space between book-end features
				if(intervalFeature.getScreenFrom() > line.get(line.size()-1).getScreenTo()+this.getGap()){ // +2 because we want some space between adjacent reads
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
	/**Print raw features under track. 
	 * windowSize size the number of characters before clipping occurs. This is 
	 * typically the window size for plotting. windowSize is used only by CLIP mode.  
	 * */
	public String printFeatures(int windowSize){
		
		if(this.getPrintMode().equals(PrintRawLine.FULL)){
			windowSize= Integer.MAX_VALUE;
		} else if(this.getPrintMode().equals(PrintRawLine.CLIP)){
			// Keep windowSize as it is
		} else {
			return "";
		} 
		
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

	protected void setIntervalFeatureSet(IntervalFeatureSet intervalFeatureSet2) {
		
	}

}
