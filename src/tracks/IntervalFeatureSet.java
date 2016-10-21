package tracks;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Set;
import java.util.regex.Pattern;

import org.apache.commons.io.FileUtils;

import exceptions.InvalidGenomicCoordsException;
import exceptions.InvalidRecordException;
import htsjdk.tribble.readers.TabixReader;
import htsjdk.tribble.readers.TabixReader.Iterator;
import samTextViewer.GenomicCoords;
import samTextViewer.Utils;
import sortBgzipIndex.MakeTabixIndex;

/**
 * Class containing sets of IntervalFeature objects. Essentially a representation of
 * a bed or gtf file. Similar to pybedtools BedTool object.
 * Essentially this is a HashMap with chrom as keys and List<IntervalFeature> as values.
 * Implementation is pretty naive but should be good enough.
 * @author berald01
 */
public class IntervalFeatureSet {
	
	/** Key is chromosome */
	private TabixReader tabixReader= null;
	private TrackFormat type;
	private String hideRegex= "^$"; // Regex to capture feature to hide/show
	private String showRegex= ".*";
	
	/* C o n s t r u c t o r */
	
	/** Construct from bed or gtf file.
	 * @throws IOException 
	 * @throws InvalidGenomicCoordsException 
	 * @throws SQLException 
	 * @throws InvalidRecordException 
	 * @throws ClassNotFoundException */
	public IntervalFeatureSet(String infile) throws IOException, InvalidGenomicCoordsException, ClassNotFoundException, InvalidRecordException, SQLException{
		
		this.type= Utils.getFileTypeFromName(new File(infile).getName());
		String sourceFile= infile;
		
		if( ! Utils.hasTabixIndex(new File(infile).getAbsolutePath())){
			// Tabix index not found for this file. Sort and index input to tmp.

			//final File tmpdir = Files.createTempDirectory("asciigenome_tmp_").toFile();	
			//sourceFile= new File(tmpdir, new File(infile).getName()).getAbsolutePath();
			
			if( ! infile.endsWith(".gz")){
				sourceFile += ".gz";
			}
			sourceFile= File.createTempFile("asciigenome_intFeatSet.", "." + new File(infile).getName()).getAbsolutePath();
			new File(sourceFile).deleteOnExit();
			new File(sourceFile + ".tbi").deleteOnExit();

			new MakeTabixIndex(infile, new File( sourceFile ), Utils.trackFormatToTabixFormat(this.type));
			
			//Runtime.getRuntime().addShutdownHook(new Thread(new Runnable() {
			//	// Delete on exit non-empty dir. See http://stackoverflow.com/questions/11165253/deleting-a-directory-on-exit-in-java
			//	@Override
			//	public void run() {
			//		FileUtils.deleteQuietly(tmpdir);
			//	}
			//}));
			
		} 

		this.tabixReader= new TabixReader(new File(sourceFile).getAbsolutePath());
			
	}
	
	/** Initialize directly from map of IntervalFeatures
	 * */
//	public IntervalFeatureSet(Map <String, List<IntervalFeature>> intervalMap, TrackFormat type){
//		this.intervalMap= intervalMap;
//		this.sortIntervalsWithinChroms();
//		this.type= type;
//	}
	
	/* M e t h o d s */

	/** Return true if a feature in the IntervalFeatureSet is visible, i.e. it
	 * passes the regex filters. Note that regex filters are applied to the raw string.
	 * */
	private boolean featureIsVisible(IntervalFeature x){
		boolean showIt= Pattern.compile(this.showRegex).matcher(x.getRaw()).find();
		boolean hideIt= false;
		if(!this.hideRegex.isEmpty()){
			hideIt= Pattern.compile(this.hideRegex).matcher(x.getRaw()).find();	
		}
		if(showIt && !hideIt){
			return true;
		} else {
			return false;
		} 
	}
	
	public List<IntervalFeature> getFeaturesInInterval(String chrom, int from, int to) throws IOException, InvalidGenomicCoordsException{

		if(from < 1){
			System.err.println("from < 1: " + from + "; resetting to 1."); 
			from= 1;
		}
		
		if(from > to || to < 1){
			System.err.println("Invalid coordinates: from: " + from + ";to: " + to 
					+ "Resetting to initial 1-" + Integer.MAX_VALUE);
			from= 1;
			to= Integer.MAX_VALUE;
			throw new InvalidGenomicCoordsException();
		}		
		List<IntervalFeature> xFeatures= new ArrayList<IntervalFeature>();
		//if(isTabix){
		Iterator qry = this.tabixReader.query(chrom,  from-1, to);
		while(true){
			String q = qry.next();
			if(q == null){
				break;
			}
			IntervalFeature intervalFeature= new IntervalFeature(q, this.type);
			xFeatures.add(intervalFeature);
		} 
		
		// Remove hidden features
		List<IntervalFeature> xFeaturesFiltered= new ArrayList<IntervalFeature>();
		for(IntervalFeature x : xFeatures){
			if(this.featureIsVisible(x)){
				xFeaturesFiltered.add(x);
			}
		}
		return xFeaturesFiltered;
	}
	
	/** add introns to the set of features in input. 
	 * */
	//private List<IntervalFeature> addIntrons(List<IntervalFeature> features){
	//	// * Collect features belonging to the same transcript_id
	//	LinkedHashMap<String, IntervalFeature> transcripts= new LinkedHashMap<String, IntervalFeature>(); 
	//	
	//	// * Walk along features and mark as introns all the gaps
	//	return features;
	//
	//}
	
	public static boolean isValidBedLine(String line){
		String[] bdg= line.split("\t");
		if(bdg.length < 3){
			return false;
		}
		try{
			Integer.parseInt(bdg[1]);
			Integer.parseInt(bdg[2]);
		} catch(NumberFormatException e){
			return false;
		}
		return true;
	}
	
//	private void sortIntervalsWithinChroms(){
//		for(String chrom : this.intervalMap.keySet()){
//			List<IntervalFeature> interalList = this.intervalMap.get(chrom);
//			Collections.sort(interalList);
//		}
//	}
		
	/** Get the next feature on chrom after "from" position or null if no 
	 * feature found 
	 * @throws IOException 
	 * @throws InvalidGenomicCoordsException */
	private IntervalFeature getNextFeatureOnChrom(String chrom, int from) throws IOException, InvalidGenomicCoordsException{
		
		Iterator iter = this.tabixReader.query(chrom, from-1, Integer.MAX_VALUE);
		while(true){
			String line= iter.next();
			if(line == null){
				return null;
			} 
			IntervalFeature x= new IntervalFeature(line, this.type);
			if(x.getFrom() > from && this.featureIsVisible(x)){
				return x;
			}
		}
	}
	
	/** Return the set chroms sorted but with and first chrom set to startChrom.
	 *  */
	protected List<String> getChromListStartingAt(Set<String> chroms, String startChrom){
		// Set:            chr1 chr2 chr3 chr4 chr5 chr6 chr7
		// StartChrom:     chr3
		// Ordered chroms: chr3 chr4 chr5 chr6 chr7 chr1 chr2
		List<String> orderedChroms= new ArrayList<String>();
		orderedChroms.addAll(chroms);
		Collections.sort(orderedChroms);
		int idx= orderedChroms.indexOf(startChrom);
		if(idx == -1){ // If startChrom is not present at all in the bed/gtf file.
			return orderedChroms;
		}
		List<String> chromsStartingAt= new ArrayList<String>();
		chromsStartingAt.addAll(orderedChroms.subList(idx, orderedChroms.size()));
		chromsStartingAt.addAll(orderedChroms.subList(0, idx));
		return chromsStartingAt;
	}
	
	/** Searching the current chrom starting at "from" to find the *next* feature matching the given string. 
	 * If not found, search the other chroms, if not found restart from the beginning of
	 * the current chrom until the "from" position is reached. 
	 * @throws InvalidGenomicCoordsException */
	protected IntervalFeature findNextRegexInGenome(String query, String chrom, int from) throws IOException, InvalidGenomicCoordsException{
		
		int startingPoint= from-1; // -1 becouse tabix.query from is 0 based (seems so at least)
		List<String> chromSearchOrder = getChromListStartingAt(this.tabixReader.getChromosomes(), chrom);
		chromSearchOrder.add(chrom);
		for(String curChrom : chromSearchOrder){
			
			Iterator iter = this.tabixReader.query(curChrom , startingPoint, Integer.MAX_VALUE);
			while(true){
				String line= iter.next();
				if(line == null) break;
				boolean matched= Pattern.compile(query).matcher(line).find();
				if(matched){
					IntervalFeature x= new IntervalFeature(line, this.type);
					if(x.getFrom() > startingPoint && this.featureIsVisible(x)){
						return x;
					}
				} 
			}
			startingPoint= 0;
		} return null; // Not found anywhere
	}
	
	/** Find all the feature matching regex.
	 * Only the feature on one chromosome are returned and this chromsome is the first one to have a match.
	 * The search starts from the beginning of the current chrom and if nothing is found continues
	 * to the other chroms. 
	 * @throws InvalidGenomicCoordsException */
	private List<IntervalFeature> findAllChromMatchInGenome(String query, GenomicCoords currentGc) throws IOException, InvalidGenomicCoordsException{
		
		// Accumulate features here
		List<IntervalFeature> matchedFeatures= new ArrayList<IntervalFeature>(); 

		// We start search from input chrom
		List<String> chromSearchOrder= null;
		chromSearchOrder = getChromListStartingAt(this.tabixReader.getChromosomes(), currentGc.getChrom());
		
		chromSearchOrder.add(currentGc.getChrom());		
		for(String curChrom : chromSearchOrder){
		
			Iterator iter = this.tabixReader.query(curChrom , 0, Integer.MAX_VALUE);
			while(true){
				String line= iter.next();
				if(line == null) break;
				boolean matched= Pattern.compile(query).matcher(line).find();
				if(matched){
					IntervalFeature x= new IntervalFeature(line, this.type);
					if(this.featureIsVisible(x)){
						matchedFeatures.add(x);
					}
				}
			}
			if(matchedFeatures.size() > 0){
				// At least one feature matching regex found on this chrom.
				// Chech we are at the same position as the beginning. if so, continue to other chroms
				if(matchedFeatures.get(0).getChrom().equals(currentGc.getChrom()) && 
				   matchedFeatures.get(0).getFrom() == currentGc.getFrom() &&
				   matchedFeatures.get(matchedFeatures.size()-1).getTo() == currentGc.getTo()){
				   // Discard results and keep searching other chroms.
					matchedFeatures= new ArrayList<IntervalFeature>();
				} else {
					break;
				}
			}
		} // Loop chrom
		return matchedFeatures;
	}
	
	/** Execute findAllChromRegexInGenome() and return the extreme coordinates of the matched features */
	protected GenomicCoords genomicCoordsAllChromMatchInGenome(String query, GenomicCoords currentGc) throws IOException, InvalidGenomicCoordsException{

		List<IntervalFeature> matchedFeatures = findAllChromMatchInGenome(query, currentGc);
		
		if(matchedFeatures.size() == 0){
			return currentGc;
		}
		
		// Now get the coords of the first and last feature matched.
		String chrom= matchedFeatures.get(0).getChrom();
		int startFrom= matchedFeatures.get(0).getFrom();
		int endTo= matchedFeatures.get(matchedFeatures.size()-1).getTo();
		GenomicCoords allMatchesGc= new GenomicCoords(
				chrom, 
				startFrom, 
				endTo, 
				currentGc.getSamSeqDict(),
				currentGc.getUserWindowSize(),
				currentGc.getFastaFile());
		return allMatchesGc;
		
	}
	
	public GenomicCoords findNextMatch(GenomicCoords currentGc, String query) throws IOException, InvalidGenomicCoordsException{

		IntervalFeature nextFeature= findNextRegexInGenome(query, currentGc.getChrom(), currentGc.getTo());
		if(nextFeature == null){
			return currentGc;
		}
		GenomicCoords nextGc= new GenomicCoords(
				nextFeature.getChrom(), 
				nextFeature.getFrom(), 
				nextFeature.getFrom() + currentGc.getGenomicWindowSize() - 1, 
				currentGc.getSamSeqDict(),
				currentGc.getUserWindowSize(),
				currentGc.getFastaFile());
		return nextGc;
	}

	protected GenomicCoords startEndOfNextFeature(GenomicCoords currentGc) throws InvalidGenomicCoordsException, IOException {
		IntervalFeature nextFeature= getNextFeatureOnChrom(currentGc.getChrom(), currentGc.getTo());
		if(nextFeature == null){
			return currentGc;
		}
		GenomicCoords nextGc= new GenomicCoords(
				nextFeature.getChrom(), 
				nextFeature.getFrom(), 
				nextFeature.getTo(), 
				currentGc.getSamSeqDict(),
				currentGc.getUserWindowSize(),
				currentGc.getFastaFile());
		return nextGc;		
	}
	
	/**Return the coordinates of the next feature so that the start coincide with the start of the feature and
	 * the end is the start + windowSize.  
	 * */
	public GenomicCoords coordsOfNextFeature(GenomicCoords currentGc) throws InvalidGenomicCoordsException, IOException {
		IntervalFeature nextFeature= getNextFeatureOnChrom(currentGc.getChrom(), currentGc.getTo());
		if(nextFeature == null){
			return currentGc;
		}
		GenomicCoords nextGc= new GenomicCoords(
				nextFeature.getChrom(), 
				nextFeature.getFrom(), 
				nextFeature.getFrom() + currentGc.getGenomicWindowSize() -1, 
				currentGc.getSamSeqDict(),
				currentGc.getUserWindowSize(),
				currentGc.getFastaFile());
		return nextGc;
	}
	
	/* S e t t e r s  and  G e t t e r s */
    
//	protected Map<String, List<IntervalFeature>> getIntervalMap() {
//		return intervalMap;
//	}

	protected void setHideRegex(String hideRegex) {
		this.hideRegex = hideRegex;
	}
	protected String getHideRegex() {
		return hideRegex;
	}

	protected void setShowRegex(String showRegex) {
		this.showRegex = showRegex;
	}
	protected String getShowRegex() {
		return showRegex;
	}

	
}
