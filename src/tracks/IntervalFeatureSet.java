package tracks;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.Reader;
import java.net.URL;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.zip.GZIPInputStream;

import org.apache.commons.validator.routines.UrlValidator;

import exceptions.InvalidGenomicCoordsException;
import htsjdk.tribble.readers.TabixReader;
import htsjdk.tribble.readers.TabixReader.Iterator;
import samTextViewer.GenomicCoords;
import samTextViewer.Utils;

/**
 * Class containing sets of IntervalFeature objects. Essentially a representation of
 * a bed or gtf file. Similar to pybedtools BedTool object.
 * Essentially this is a HashMap with chrom as keys and List<IntervalFeature> as values.
 * Implementation is pretty naive but should be good enough.
 * @author berald01
 */
public class IntervalFeatureSet {
	
	private Map <String, List<IntervalFeature>> intervalMap; // new HashMap <String, List<IntervalFeature>>(); 
	private TabixReader tabixReader= null;
	private boolean isTabix= false;
	private TrackFormat type;
	private String hideRegex= ""; // Regex to capture feature to hide/show
	private String showRegex= ".*";
	
	/* C o n s t r u c t o r */
	
	/** Construct from bed or gtf file.
	 * @throws IOException */
	public IntervalFeatureSet(String infile) throws IOException{
		
		this.type= Utils.getFileTypeFromName(new File(infile).getName());
		
		if(Utils.hasTabixIndex(new File(infile).getAbsolutePath())){
			this.tabixReader= new TabixReader(new File(infile).getAbsolutePath());
			this.isTabix= true;
		} else {
			this.intervalMap= loadFileIntoIntervalMap(infile);
			this.sortIntervalsWithinChroms();
			this.isTabix= false;
		}
	}
	
	/* M e t h o d s */

	/** Return true if a feature in the IntervalFeatureSet is visible, i.e. it
	 * passes the regex filters. Note that regex filters are applied to the raw string.
	 * */
	private boolean featureIsVisible(IntervalFeature x){
		if(x.getRaw().matches(this.showRegex) && !x.getRaw().matches(this.hideRegex)){
			return true;
		} else {
			return false;
		} 
	}
	
	public List<IntervalFeature> getFeaturesInInterval(String chrom, int from, int to) throws IOException{
		if(from > to || from < 1 || to < 1){
			throw new RuntimeException("Invalid range: " + from + "-" + to);
		}		
		List<IntervalFeature> xFeatures= new ArrayList<IntervalFeature>();
		if(isTabix){
			Iterator qry = this.tabixReader.query(chrom,  from, to);
			while(true){
				String q = qry.next();
				if(q == null){
					break;
				}
				IntervalFeature intervalFeature= new IntervalFeature(q, this.type);
				xFeatures.add(intervalFeature);
			}
		} else {
			List<IntervalFeature> thisChrom= this.intervalMap.get(chrom);		
			if(thisChrom == null){
				return xFeatures;
			}
			/* Overlap scenarios
			             from      to
			               |-------|
			     --------************-----  4.
			     ---------****------------- 1.
			     ------------****---------- 3.
			 	 -----------------****----- 2.
			     ----------------------***- Break iterating */
			for(IntervalFeature x : thisChrom){
				if( (x.getFrom() >= from && x.getFrom() <= to)  // 2. 3.
					|| (x.getTo() >= from && x.getTo() <= to)   // 1. 3.
					|| (x.getFrom() <= from && x.getTo() >= to)   // 4.
					){
					xFeatures.add(x);
				}
				if(x.getFrom() > to){
					break;
				}
			}
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
	
	private Map <String, List<IntervalFeature>> loadFileIntoIntervalMap(String infile) throws IOException{
		
		System.err.print("Reading file '" + infile + "'...");
		
		Map <String, List<IntervalFeature>> intervalMap= new HashMap<String, List<IntervalFeature>>(); 
		
		BufferedReader br= null;
		InputStream gzipStream= null;
		UrlValidator urlValidator = new UrlValidator();
		if(infile.endsWith(".gz")){
			if(urlValidator.isValid(infile)) {
				gzipStream = new GZIPInputStream(new URL(infile).openStream());
			} else {
				gzipStream = new GZIPInputStream(new FileInputStream(infile));
			}
			Reader decoder = new InputStreamReader(gzipStream, "UTF-8");
			br = new BufferedReader(decoder);
		} else if(urlValidator.isValid(infile)) {
			InputStream instream= new URL(infile).openStream();
			Reader decoder = new InputStreamReader(instream, "UTF-8");
			br = new BufferedReader(decoder);
		} else {
			br = new BufferedReader(new FileReader(infile));
		}
		String line;
		boolean isFirst= true;
		while ((line = br.readLine()) != null) {
			if(line.trim().equals("##FASTA")){
				// Stop reading once you hit the optional sequence section in gff3 format
				// See http://gmod.org/wiki/GFF#GFF3_Sequence_Section
				break;
			}
			if(line.trim().startsWith("#") || line.trim().isEmpty()){
				continue;
			}
			if(isFirst && this.type.equals(TrackFormat.BED)){
				isFirst= false;
				if(!isValidBedLine(line)){ // Allow first line to fail: Might be header.
					System.err.print("First line skipped. ");
					continue;
				}
			}
			IntervalFeature f= new IntervalFeature(line, this.type);
			if(intervalMap.containsKey(f.getChrom())){
				intervalMap.get(f.getChrom()).add(f);
			} else {
				List<IntervalFeature> il= new ArrayList<IntervalFeature>();
				il.add(f);
				intervalMap.put(f.getChrom(), il);
			}
		}
		br.close();
		if(gzipStream != null){
			gzipStream.close();
		}
		System.err.println(" Done");
		return intervalMap;
	} 
	
	private void sortIntervalsWithinChroms(){
		for(String chrom : this.intervalMap.keySet()){
			List<IntervalFeature> interalList = this.intervalMap.get(chrom);
			Collections.sort(interalList);
		}
	}
		
	/** Get the next feature on chrom after "from" position or null if no 
	 * feature found 
	 * @throws IOException */
	private IntervalFeature getNextFeatureOnChrom(String chrom, int from) throws IOException{
		
		if(this.intervalMap != null){
			List<IntervalFeature> featuresList = this.intervalMap.get(chrom);
			if(featuresList == null){ // No features at all on this chrom.
				return null;
			}
			for(IntervalFeature x : featuresList){
				if(x.getFrom() > from && this.featureIsVisible(x)){
					return x;
				}
			}
			return null;
		} else if(this.isTabix){
			Iterator iter = this.tabixReader.query(chrom, from, Integer.MAX_VALUE);
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
		return null;
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
	 * the current chrom until the "from" position is reached. */
	protected IntervalFeature findNextRegexInGenome(String regex, String chrom, int from) throws IOException{
		
		if(this.intervalMap != null){
	
			// We start search from input chrom and starting position. We'll return to 
			// the start of this chrom only after having searched all the other chroms.
			int startingPoint= from;
			List<String> chromSearchOrder = getChromListStartingAt(this.intervalMap.keySet(), chrom);
			chromSearchOrder.add(chrom);
			for(String curChrom : chromSearchOrder){
				
				List<IntervalFeature> featuresList = this.intervalMap.get(curChrom);
				for(IntervalFeature x : featuresList){
					//if(x.getFrom() > startingPoint && x.getRaw().toLowerCase().contains(regex.toLowerCase())){
					if(x.getFrom() > startingPoint && x.getRaw().matches(regex) && this.featureIsVisible(x)){ 
						return x;
					}
				}
				startingPoint= 0;
				
			} return null; // Not found anywhere

		} else if(this.isTabix){
		
			int startingPoint= from;
			List<String> chromSearchOrder = getChromListStartingAt(this.tabixReader.getChromosomes(), chrom);
			chromSearchOrder.add(chrom);
			for(String curChrom : chromSearchOrder){
				
				Iterator iter = this.tabixReader.query(curChrom , startingPoint, Integer.MAX_VALUE);
				while(true){
					String line= iter.next();
					if(line == null) break;
					if(line.matches(regex)){
						IntervalFeature x= new IntervalFeature(line, this.type);
						if(x.getFrom() > startingPoint && this.featureIsVisible(x)){
							return x;
						}
					}
				}
				startingPoint= 0;
			} return null; // Not found anywhere
		}
		return null;
	}
	
	/** Find all the feature matching regex.
	 * Only the feature on one chromosome are returned and this chromsome is the first one to have a match.
	 * The search starts from the beginning of the current chrom and if nothing is found continues
	 * to the other chroms. 
	 * @throws InvalidGenomicCoordsException */
	private List<IntervalFeature> findAllChromRegexInGenome(String regex, GenomicCoords currentGc) throws IOException, InvalidGenomicCoordsException{
		
		// Accumulate features here
		List<IntervalFeature> matchedFeatures= new ArrayList<IntervalFeature>(); 

		// We start search from input chrom
		List<String> chromSearchOrder= null;
		if(this.intervalMap != null){
			chromSearchOrder = getChromListStartingAt(this.intervalMap.keySet(), currentGc.getChrom());
		} else if (this.isTabix){
			chromSearchOrder = getChromListStartingAt(this.tabixReader.getChromosomes(), currentGc.getChrom());
		} else {
			throw new RuntimeException("Cannot init chroms");
		}
		
		chromSearchOrder.add(currentGc.getChrom());		
		for(String curChrom : chromSearchOrder){
		
			if(this.intervalMap != null){
				List<IntervalFeature> featuresList = this.intervalMap.get(curChrom);
				for(IntervalFeature x : featuresList){
					if(x.getRaw().matches(regex) && this.featureIsVisible(x)){ 
						matchedFeatures.add(x);
					}
				}
			} else if(this.isTabix){
				Iterator iter = this.tabixReader.query(curChrom , 0, Integer.MAX_VALUE);
				while(true){
					String line= iter.next();
					if(line == null) break;
					if(line.matches(regex)){
						IntervalFeature x= new IntervalFeature(line, this.type);
						if(this.featureIsVisible(x)){
							matchedFeatures.add(x);
						}
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
	protected GenomicCoords genomicCoordsAllChromRegexInGenome(String regex, GenomicCoords currentGc) throws IOException, InvalidGenomicCoordsException{

		List<IntervalFeature> matchedFeatures = findAllChromRegexInGenome(regex, currentGc);
		
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
	
	public GenomicCoords findNextRegex(GenomicCoords currentGc, String regex) throws IOException, InvalidGenomicCoordsException{

		IntervalFeature nextFeature= findNextRegexInGenome(regex, currentGc.getChrom(), currentGc.getTo());
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
    
	protected Map<String, List<IntervalFeature>> getIntervalMap() {
		return intervalMap;
	}

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
