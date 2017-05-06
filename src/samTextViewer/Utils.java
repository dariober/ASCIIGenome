package samTextViewer;

import java.io.BufferedReader;
import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.io.Reader;
import java.net.HttpURLConnection;
import java.net.MalformedURLException;
import java.net.URL;
import java.nio.charset.StandardCharsets;
import java.nio.file.FileSystems;
import java.nio.file.FileVisitResult;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.PathMatcher;
import java.nio.file.Paths;
import java.nio.file.SimpleFileVisitor;
import java.nio.file.attribute.BasicFileAttributes;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;

import org.apache.commons.lang.StringUtils;
import org.apache.commons.lang3.text.StrMatcher;
import org.apache.commons.lang3.text.StrTokenizer;
import org.apache.commons.validator.routines.UrlValidator;
import org.broad.igv.bbfile.BBFileReader;
import org.broad.igv.bbfile.BedFeature;
import org.broad.igv.bbfile.BigBedIterator;
import org.broad.igv.bbfile.BigWigIterator;
import org.broad.igv.bbfile.WigItem;
import org.broad.igv.tdf.TDFReader;

import com.google.common.base.Joiner;
import com.google.common.base.Splitter;
import com.google.common.base.Strings;
import com.google.common.collect.Lists;
import com.google.gson.JsonArray;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.google.gson.JsonParser;

import exceptions.InvalidCommandLineException;
import exceptions.InvalidGenomicCoordsException;
import exceptions.InvalidRecordException;
import faidx.Faidx;
import faidx.UnindexableFastaFileException;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.filter.AggregateFilter;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.tribble.index.tabix.TabixFormat;
import htsjdk.tribble.readers.TabixReader;
import tracks.IntervalFeature;
import tracks.Track;
import tracks.TrackFormat;
import ucsc.UcscGenePred;

/**
 * @author berald01
 *
 */
public class Utils {
	
	/** Returns true of the list of arguments argList contains the given flag
	 * IMPORTANT SIDE EFFECT: If found, the argument flag is removed from the list. 
	 * */
	public static boolean argListContainsFlag(List<String> argList, String flag){
		boolean hasFlag= false;
		if(argList.contains(flag)){
			argList.remove(flag);
			hasFlag= true;
		}
		return hasFlag;
	}

	/** Check if the list of arguments contains "param" and if so return the argument. 
	 * Returns null if arglist does not contain the parameter
	 * IMPORTANT SIDE EFFECT: If found, the parameter and its argument are removed from argList. 
	 * @throws InvalidCommandLineException 
	 * */
	public static String getArgForParam(List<String> argList, String param, String defArg) throws InvalidCommandLineException{
		
		int idx= argList.indexOf(param);
		if(idx == -1){
			return defArg;
		}
		if(idx == (argList.size() - 1)){
			// If param is the last item in the list you cannot get its argument!
			throw new InvalidCommandLineException();
		}
		String arg= argList.get(idx+1);
		
		// Remove param and arg
		argList.remove(idx);
		argList.remove(idx);
		
		return arg;

	}

	public static void checkFasta(String fasta, int debug) throws IOException, UnindexableFastaFileException {
		if(fasta == null){
			return;
		}
		File fafile= new File(fasta);
		if(!fafile.isFile()){
			System.err.println("Fasta file '" + fasta + "' not found.");
			if(debug == 0 || debug == 1){
				System.exit(1);
			} else if (debug == 2){
				throw new IOException();
			}
		} 
		if(!fafile.canRead()){
			System.err.println("Fasta file '" + fasta + "' is not readable.");
			if(debug == 0 || debug == 1){
				System.exit(1);
			} else if (debug == 2){
				throw new IOException();
			}
		}
		
		IndexedFastaSequenceFile faSeqFile = null;
		try {
			faSeqFile= new IndexedFastaSequenceFile(fafile);
			faSeqFile.close();
		} catch (FileNotFoundException e) {
			System.err.println("\nIndexing '" + fasta + "'.");
			new Faidx(new File(fasta));
			(new File(fasta + ".fai")).deleteOnExit();
		}
	}
	
    public static long getAlignedReadCount(String bam) throws IOException{

		/*  ------------------------------------------------------ */
		/* This chunk prepares SamReader from local bam or URL bam */
		UrlValidator urlValidator = new UrlValidator();
		SamReaderFactory srf=SamReaderFactory.make();
		srf.validationStringency(ValidationStringency.SILENT);
		SamReader samReader;
		if(urlValidator.isValid(bam)){
			samReader = SamReaderFactory.makeDefault().open(
					SamInputResource.of(new URL(bam)).index(new URL(bam + ".bai"))
			);
		} else {
			samReader= srf.open(new File(bam));
		}
		/*  ------------------------------------------------------ */

		List<SAMSequenceRecord> sequences = samReader.getFileHeader().getSequenceDictionary().getSequences();
		long alnCount= 0;
		for(SAMSequenceRecord x : sequences){
			alnCount += samReader.indexing().getIndex().getMetaData(x.getSequenceIndex()).getAlignedRecordCount();
		}
		samReader.close();
		return alnCount;
    }

	/** Merge overlapping features. If screenCoords is true, merge is based on the screen coordinates. Otherwise use
	 * genomic coordinates.  
	 * */
	public static List<IntervalFeature> mergeIntervalFeatures(List<IntervalFeature> intervalList, boolean screenCoords) throws InvalidGenomicCoordsException{
		List<IntervalFeature> mergedList= new ArrayList<IntervalFeature>();		 
		if(intervalList.size() == 0){
			return mergedList;
		}
		
		String mergedChrom = null;
		int mergedFrom= -1;
		int mergedTo= -1;
		int mergedScreenFrom= -1;
		int mergedScreenTo= -1;
		int numMrgIntv= 1; // Number of intervals in the merged one. 
		Set<Character> strand= new HashSet<Character>(); // Put here all the different strands found in the merged features.		

		for(int i= 0; i < (intervalList.size()+1); i++){
			// We do an additional loop to add to the mergedList the last interval.
			// The last loop has interval == null so below you need to account for it
			IntervalFeature interval= null;
			if(i < intervalList.size()){
				interval= intervalList.get(i); 
			}
			
			if(mergedFrom < 0){ // Init interval
				mergedChrom= interval.getChrom(); 
				mergedFrom= interval.getFrom();
				mergedTo= interval.getTo();
				mergedScreenFrom= interval.getScreenFrom();
				mergedScreenTo= interval.getScreenTo();
				continue;
			}
			// Sanity check: The list to be merged is on the same chrom and sorted by start pos.
			if(i < intervalList.size() && (!mergedChrom.equals(interval.getChrom()) || mergedFrom > interval.getFrom() || mergedFrom > mergedTo)){
				System.err.println(mergedChrom + " " + mergedFrom + " " + mergedTo);
				throw new RuntimeException();
			} 
			
			boolean overlap= false; 
			if( screenCoords && i < intervalList.size() && (mergedScreenFrom <= interval.getScreenTo() && mergedScreenTo >= (interval.getScreenFrom()-1)) ){ 
				// Overlap: Extend <to> coordinate. See also http://stackoverflow.com/questions/325933/determine-whether-two-date-ranges-overlap
				overlap= true;
			} else if(i < intervalList.size() && (mergedFrom <= interval.getTo() && mergedTo >= (interval.getFrom()-1) )){ 
				// Overlap on genomic coordinates
				overlap= true;
			}
			if(overlap){ // Overlap found in screen or genomic coords.  
			    mergedTo= Math.max(interval.getTo(), mergedTo);
			    mergedScreenTo= Math.max(interval.getScreenTo(), mergedScreenTo);
				strand.add(interval.getStrand());
				numMrgIntv++;
				overlap= false;
			} 
			else {
				// No overlap add merged interval to list and reset new merged interval
				IntervalFeature x= new IntervalFeature(mergedChrom + "\t" + (mergedFrom-1) + "\t" + mergedTo, TrackFormat.BED);
				x.setScreenFrom(mergedScreenFrom);
				x.setScreenTo(mergedScreenTo);
				if(strand.size() == 1){
					x.setStrand(strand.iterator().next());
				} 
				strand.clear();

				if(x.equals(intervalList.get(i-1)) && numMrgIntv == 1){
					mergedList.add(intervalList.get(i-1));
				} else {
					mergedList.add(x);
				}
				
				if(i < intervalList.size()){
					// Do not reset from/to if you are in extra loop.
					mergedFrom= interval.getFrom();
					mergedTo= interval.getTo();
					mergedScreenFrom= interval.getScreenFrom();
					mergedScreenTo= interval.getScreenTo();
					strand.add(interval.getStrand());
					numMrgIntv= 1;
				}
			}
		}
		return mergedList;
	}
	
//	public static LinkedHashMap<String, Integer> xterm256ColorCodes(){
//		// See http://misc.flogisoft.com/bash/tip_colors_and_formatting
//		// From http://jonasjacek.github.io/colors/
//		LinkedHashMap<String, Integer> colourCodes= new LinkedHashMap<String, Integer>();
//		colourCodes.put("default", 39);
//		colourCodes.put("black", 30);
//		colourCodes.put("red", 31);
//		colourCodes.put("green", 32);
//		colourCodes.put("yellow", 33);
//		colourCodes.put("blue", 34);
//		colourCodes.put("magenta", 35);
//		colourCodes.put("cyan", 36);
//		colourCodes.put("light_grey", 37);
//		colourCodes.put("grey", 90);
//		colourCodes.put("light_red", 91);
//		colourCodes.put("light_green", 92);
//		colourCodes.put("light_yellow", 93);
//		colourCodes.put("light_blue", 94);
//		colourCodes.put("light_magenta", 95);
//		colourCodes.put("light_cyan", 96);
//		colourCodes.put("white", 97);
//		// To be continued
//		return colourCodes;
//	}
	
//	public static Color ansiColourToGraphicsColor(int ansiColor) throws InvalidColourException{
//		
//		if(!xterm256ColorCodes().entrySet().contains(ansiColor)){
//			throw new InvalidColourException();
//		}
//		if(ansiColor == 30){ return Color.BLACK; }
//		if(ansiColor == 31 || ansiColor == 91){ return Color.RED; }
//		if(ansiColor == 32 || ansiColor == 92){ return Color.GREEN; }
//		if(ansiColor == 33 || ansiColor == 93){ return Color.YELLOW; }
//		if(ansiColor == 34 || ansiColor == 94){ return Color.BLUE; }
//		if(ansiColor == 35 || ansiColor == 95){ return Color.MAGENTA; }
//		if(ansiColor == 36 || ansiColor == 96){ return Color.CYAN; }
//		if(ansiColor == 37){ return Color.LIGHT_GRAY; }
//		if(ansiColor == 90){ return Color.DARK_GRAY; }
//		if(ansiColor == 97){ return Color.WHITE; }
//		return Color.BLACK;
//	}
	
	/** Return true if fileName has a valid tabix index. 
	 * @throws IOException 
	 * */
	public static boolean hasTabixIndex(String fileName) throws IOException{
		
		if((new UrlValidator()).isValid(fileName) && fileName.startsWith("ftp")){
			// Because of issue #51
			return false;
		}
		
		try{
			TabixReader tabixReader= new TabixReader(fileName);
			tabixReader.readLine();
			tabixReader.close();
			return true;
		} catch (Exception e){
			return false;
		}
	}
	
	/** Get the first chrom string from first line of input file. As you add support for more filetypes you should update 
	 * this function. This method is very dirty and shouldn't be trusted 100% 
	 * @throws InvalidGenomicCoordsException 
	 * @throws SQLException 
	 * @throws InvalidRecordException 
	 * @throws InvalidCommandLineException 
	 * @throws ClassNotFoundException */
	public static String initRegionFromFile(String x) throws IOException, InvalidGenomicCoordsException, ClassNotFoundException, InvalidCommandLineException, InvalidRecordException, SQLException{
		UrlValidator urlValidator = new UrlValidator();
		String region= "";
		TrackFormat fmt= Utils.getFileTypeFromName(x); 
		
		if(fmt.equals(TrackFormat.BAM)){
			
			SamReader samReader;
			if(urlValidator.isValid(x)){
				samReader = SamReaderFactory.makeDefault().open(SamInputResource.of(new URL(x)));
			} else {			
				SamReaderFactory srf=SamReaderFactory.make();
				srf.validationStringency(ValidationStringency.SILENT);
				samReader = srf.open(new File(x));
			}
			region= samReader.getFileHeader().getSequence(0).getSequenceName();
			samReader.close();
			return region;
		
		} else if(fmt.equals(TrackFormat.BIGWIG) && !urlValidator.isValid(x)){
			// Loading from URL is painfully slow so do not initialize from URL
			return initRegionFromBigWig(x);
			
		} else if(fmt.equals(TrackFormat.BIGBED) && !urlValidator.isValid(x)){
			// Loading from URL is painfully slow so do not initialize from URL
			return initRegionFromBigBed(x);

		} else if(urlValidator.isValid(x) && (fmt.equals(TrackFormat.BIGWIG) || fmt.equals(TrackFormat.BIGBED))){
			System.err.println("Refusing to initialize from URL");
			throw new InvalidGenomicCoordsException();
			
		} else if(fmt.equals(TrackFormat.TDF)){
			Iterator<String> iter = TDFReader.getReader(x).getChromosomeNames().iterator();
			while(iter.hasNext()){
				region= iter.next();
				if(!region.equals("All")){
					return region;
				}
			} 
			System.err.println("Cannot initialize from " + x);
			throw new RuntimeException();
		
		} else if(Utils.isUcscGenePredSource(x)){
			return initRegionFromUcscGenePredSource(x);
			
		} else {
			// Input file appears to be a generic interval file. We expect chrom to be in column 1
			BufferedReader br;
			GZIPInputStream gzipStream;
			if(x.toLowerCase().endsWith(".gz") || x.toLowerCase().endsWith(".bgz")){
				if(urlValidator.isValid(x)) {
					gzipStream = new GZIPInputStream(new URL(x).openStream());
				} else {
					InputStream fileStream = new FileInputStream(x);
					gzipStream = new GZIPInputStream(fileStream);
				}
				Reader decoder = new InputStreamReader(gzipStream, "UTF-8");
				br = new BufferedReader(decoder);
			} else {
				if(urlValidator.isValid(x)) {
					InputStream instream= new URL(x).openStream();
					Reader decoder = new InputStreamReader(instream, "UTF-8");
					br = new BufferedReader(decoder);
				} else {
					br = new BufferedReader(new FileReader(x));
				}
			}
			String line;
			while ((line = br.readLine()) != null){
				line= line.trim();
				if(line.startsWith("#") || line.isEmpty() || line.startsWith("track ")){
					continue;
				}
				IntervalFeature feature= new IntervalFeature(line, fmt);
				region= feature.getChrom() + ":" + feature.getFrom(); 
				br.close();
				return region;
			}
		} 
		System.err.println("Cannot initialize from " + x);
		throw new RuntimeException();
	}
	
	private static String initRegionFromBigBed(String bigBedFile) throws IOException{
		
		BBFileReader reader= new BBFileReader(bigBedFile);
		if(! reader.isBigBedFile()){
			System.err.println("File " + bigBedFile + " is not bigBed.");
			throw new RuntimeException();
		}
		String region= reader.getChromosomeNames().get(0); // Just get chrom to start with
		
		for(String chrom : reader.getChromosomeNames()){
			BigBedIterator iter = reader.getBigBedIterator(chrom, 0, chrom, Integer.MAX_VALUE, false);
			if(iter.hasNext()){
				BedFeature x= (BedFeature) iter.next();
				region= x.getChromosome() + ":" + (x.getStartBase() + 1);
				reader.close();
				return region;
			}
		}
		reader.close();
		return region;
	}
	
	private static String initRegionFromBigWig(String bigWigFile) throws IOException{
		
		BBFileReader reader= new BBFileReader(bigWigFile);
		if(! reader.isBigWigFile()){
			System.err.println("File " + bigWigFile + " is not bigWig.");
			throw new RuntimeException();
		}
		String region= reader.getChromosomeNames().get(0); // Just get chrom to start with
		
		for(String chrom : reader.getChromosomeNames()){
			BigWigIterator iter = reader.getBigWigIterator(chrom, 0, chrom, Integer.MAX_VALUE, false);
			if(iter.hasNext()){
				WigItem x = iter.next();
				region= x.getChromosome() + ":" + (x.getStartBase() + 1);
				reader.close();
				return region;
			}
		}
		reader.close();
		return region;
	}
	
	private static String initRegionFromUcscGenePredSource(String x) throws ClassNotFoundException, IOException, InvalidCommandLineException, InvalidGenomicCoordsException, InvalidRecordException, SQLException {
		
		String xfile= new UcscGenePred(x, 1).getTabixFile();
		GZIPInputStream gzipStream;
		InputStream fileStream = new FileInputStream(xfile);
		gzipStream = new GZIPInputStream(fileStream);
		Reader decoder = new InputStreamReader(gzipStream, "UTF-8");
		BufferedReader br = new BufferedReader(decoder);
		String line= br.readLine();
		br.close();
		List<String> xlist= Lists.newArrayList(Splitter.on("\t").split(line));
		String region= xlist.get(0) + ":" + xlist.get(3) + "-" + xlist.get(4);
		return region;
		
	}

	public static boolean bamHasIndex(String bam) throws IOException{

		/*  ------------------------------------------------------ */
		/* This chunk prepares SamReader from local bam or URL bam */
		UrlValidator urlValidator = new UrlValidator();
		SamReaderFactory srf=SamReaderFactory.make();
		srf.validationStringency(ValidationStringency.SILENT);
		SamReader samReader;
		if(urlValidator.isValid(bam)){
			samReader = SamReaderFactory.makeDefault().open(
					SamInputResource.of(new URL(bam)).index(new URL(bam + ".bai"))
			);
		} else {
			samReader= srf.open(new File(bam));
		}
		/*  ------------------------------------------------------ */

		// SamReaderFactory srf=SamReaderFactory.make();
		// srf.validationStringency(ValidationStringency.SILENT);
		// SamReader samReader = srf.open(new File(bam));
		boolean hasIndex= samReader.hasIndex();
		samReader.close();
		return hasIndex;
		
	}
	
	/*public static GenomicCoords findNextStringOnFile(String string, String filename, GenomicCoords curGc,
			Map<String, IntervalFeatureSet> intervalFiles ) throws InvalidGenomicCoordsException, IOException{

		String chosenFn= "";
		if(filename.isEmpty() && intervalFiles.size() == 1){ // Only one file to chose from: Get that one
			chosenFn= new ArrayList<String>(intervalFiles.keySet()).get(0);
		} else {
			// Try to match file perfectly as it was added from cli, including path if any
			for(String fn : intervalFiles.keySet()){
				if(fn.equals(filename)){
					chosenFn = fn;
				}
			}
			if(chosenFn.isEmpty()){
				// Or try to match only file name.
				for(String fn : intervalFiles.keySet()){ // Do not look for a perfect match since the original input might contain path. 
					String onlyName= new File(fn).getName();
					if(onlyName.equals(filename)){
						chosenFn = fn;
					}
				}
			}
		}
		if(chosenFn.isEmpty()){
			System.err.println("File " + filename + " not found in file set\n" + intervalFiles.keySet());
			return curGc;
		}
		return intervalFiles.get(chosenFn).findNextString(curGc, string);
	} */
	
	public static TrackFormat getFileTypeFromName(String fileName){
		fileName= fileName.toLowerCase();
		
		if(    fileName.endsWith(".bed") 
		    || fileName.endsWith(".bed.gz") 
		    || fileName.endsWith(".bed.gz.tbi")
		    || fileName.endsWith(".bed.bgz")
		    || fileName.endsWith(".bed.bgz.tbi")){
			return TrackFormat.BED;
		
		} else if( fileName.endsWith(".gtf") 
				|| fileName.endsWith(".gtf.gz")
				|| fileName.endsWith(".gtf.gz.tbi")
				|| fileName.endsWith(".gtf.bgz")
				|| fileName.endsWith(".gtf.bgz.tbi")){
			return TrackFormat.GTF;
			
	    } else if(
				fileName.endsWith(".gff") 
				|| fileName.endsWith(".gff.gz") 
				|| fileName.endsWith(".gff.gz.tbi")
				|| fileName.endsWith(".gff.bgz") 
				|| fileName.endsWith(".gff.bgz.tbi")
				|| fileName.endsWith(".gff3")
				|| fileName.endsWith(".gff3.gz") 
				|| fileName.endsWith(".gff3.gz.tbi")
				|| fileName.endsWith(".gff3.bgz") 
				|| fileName.endsWith(".gff3.bgz.tbi")){
			return TrackFormat.GFF;
			
		} else if(fileName.endsWith(".bam") || fileName.endsWith(".cram")){
			return TrackFormat.BAM;
			
		} else if(fileName.endsWith(".bigwig") || fileName.endsWith(".bw")) {
			return TrackFormat.BIGWIG;
		} else if(fileName.endsWith(".bigbed") || fileName.endsWith(".bb")) {
			return TrackFormat.BIGBED;
		} else if(fileName.endsWith(".tdf")) {
			return TrackFormat.TDF;
		} else if(fileName.endsWith(".bedgraph.gz") || fileName.endsWith(".bedgraph")) {
			return TrackFormat.BEDGRAPH;
		} else if(fileName.endsWith(".vcf.gz") 
				|| fileName.endsWith(".vcf")
				|| fileName.endsWith(".vcf.bgz")){
			return TrackFormat.VCF;
		} else {
			// System.err.println("Unsopported file: " + fileName);
			return TrackFormat.BED;
		}
	}
	
//	public static LinkedHashMap<String, IntervalFeatureSet> createIntervalFeatureSets(List<String> fileNames) throws IOException, InvalidGenomicCoordsException, ClassNotFoundException, InvalidRecordException, SQLException{
//		LinkedHashMap<String, IntervalFeatureSet> ifsets= new LinkedHashMap<String, IntervalFeatureSet>();
//		for(String x : fileNames){
//			String f= x;
//			if(getFileTypeFromName(x).equals(TrackFormat.BED) || getFileTypeFromName(x).equals(TrackFormat.GFF)){
//				if(!ifsets.containsKey(x)){ // If the input has duplicates, do not reload duplicates!
//					IntervalFeatureSet ifs= new IntervalFeatureSet(f);
//					ifsets.put(x, ifs);
//				}
//			}
//		}
//		return ifsets;
//	}
	
    /** 
     * Transpose list of list as if they were a table. No empty cells should be present. 
     * See http://stackoverflow.com/questions/2941997/how-to-transpose-listlist
     * FROM
     * [[a, b, c, d], [a, b, c, d], [a, b, c, d]] 
     * TO
     * [[a, a, a, a],
     *  [b, b, b, b],
     *  [c, c, c, c]]
     * @param table Table like list of lists, no empty cells.
     * @return
     */
    public static <T> List<List<T>> transpose(List<List<T>> table) {
    	List<List<T>> ret = new ArrayList<List<T>>();
        final int N = table.get(0).size();
        for (int i = 0; i < N; i++) {
            List<T> col = new ArrayList<T>();
            for (List<T> row : table) {
                col.add(row.get(i));
            }
            ret.add(col);
        }
        return ret;
    }
	
    /**
     * Map list of values mapping top genomic positions to a smaller list of positions by averaging 
     * values mapped to the same reference position. These averages are the values that will be used on the 
     * y-axis.
     * @param values Values to collapse
     * @param valuePositions Genomic position of the values
     * @param referencePositions Arrival positions. I.e. where the valuePositions should be mapped to. 
     * Typically this is obtained from GenomicCoords.getMapping();  
     * @return
     */
    public static List<Double> collapseValues(List<Double> values, List<Integer> valuePositions, 
    		List<Double> referencePositions){

    	// First store here all the values mapping to each reference position, then take the average.
    	LinkedHashMap<Integer, List<Double>> zlist= new LinkedHashMap<Integer, List<Double>>();
    	for(int i= 0; i < referencePositions.size(); i++){
    		zlist.put(i, new ArrayList<Double>());
    	} 
    	
    	for(int i= 0; i < valuePositions.size(); i++){
    		if(values.get(i) == null){
    			continue;
    		}
    		int pos= valuePositions.get(i);
    		if(pos >= referencePositions.get(0) && pos <= referencePositions.get(referencePositions.size()-1)){
	    	// Do not consider data points outside screenMap.
    			int j= Utils.getIndexOfclosestValue(pos, referencePositions);
    			zlist.get(j).add((double)values.get(i));
    		}
    	}
    	
    	List<Double> compressed= new ArrayList<Double>();
    	for(int i= 0; i < referencePositions.size(); i++){
    		compressed.add(Utils.calculateAverage(zlist.get(i)));
    	}
    	return compressed;
    }

	/**
	 * Get sequence as byte[] for the given genomic coords.
	 * @param fasta
	 * @param gc
	 * @return
	 * @throws IOException
	 */
	public static byte[] prepareRefSeq(String fasta, GenomicCoords gc) throws IOException{

		byte[] faSeq= null;
		if(fasta != null){
			IndexedFastaSequenceFile faSeqFile = null;
			try {
				faSeqFile = new IndexedFastaSequenceFile(new File(fasta));
				try{
					faSeq= faSeqFile.getSubsequenceAt(gc.getChrom(), gc.getFrom(), gc.getTo()).getBases();
				} catch (NullPointerException e){
					System.err.println("Cannot fetch sequence " + gc.toString());
					e.printStackTrace();
				}
				faSeqFile.close();
			} catch (FileNotFoundException e) {
				e.printStackTrace();
			}
		}
		return faSeq;
	}     
	
	private static int parseStringToIntWithUnits(String x){
		x= x.trim();
		int multiplier= 0;
		if(x.endsWith("k") || x.endsWith("K")){
			multiplier= 1000;
			x= x.substring(0, x.length()-1).trim();
		} else if(x.endsWith("m") || x.endsWith("M")){
			multiplier= 1000000;
			x= x.substring(0, x.length()-1).trim();
		} else if(x.matches("^\\-{0,1}\\d+$") || x.matches("^\\+{0,1}\\d+$")){
			multiplier= 1;
		} else {
			throw new RuntimeException("Invalid string to convert to int: " + x);
		}
		long pos= (long) (Double.parseDouble(x) * multiplier);
		if( pos >= Integer.MAX_VALUE ){
			System.err.println("Invalid end coordinate: " + pos);
			pos= Integer.MAX_VALUE;
		}
		return (int)pos;
	}
	
	/**
	 * Parse user to modify the current genomics coordinates in input to new ones
	 * to move.
	 * @param bam 
	 * @return
	 * @throws IOException 
	 * @throws InvalidGenomicCoordsException 
	 */
	public static String parseConsoleInput(List<String> tokens, GenomicCoords gc) throws InvalidGenomicCoordsException, IOException{
		
		String region= "";
		String chrom= gc.getChrom();
		Integer from= gc.getFrom();
		Integer to= gc.getTo();
		
		int windowSize= to - from + 1;
		int halfWindow= (int)Math.rint(windowSize / 2d);
		if(tokens.get(0).equals("ff")){				
			from += halfWindow; 
			to += halfWindow;
			if(gc.getSamSeqDict() != null && !gc.getSamSeqDict().isEmpty()){
				int chromLen= gc.getSamSeqDict().getSequence(chrom).getSequenceLength();
				if(to > chromLen){
					to= chromLen;
					from= to - gc.getUserWindowSize() + 1;
				}
			}			
			return chrom + ":" + from + "-" + to;
		} else if(tokens.get(0).equals("bb")) {
			from -= halfWindow;
			to -= halfWindow; 
			if(from < 1){
				from= 1;
				to= from + gc.getUserWindowSize() - 1;
			}
			return chrom + ":" + from + "-" + to;
		} else if(tokens.get(0).equals("f")){
			int step= (int)Math.rint(windowSize / 10d);
			if(tokens.size() > 1){
				try{
					step= (int)Math.rint(windowSize * Double.parseDouble(tokens.get(1)));
				} catch(NumberFormatException e){
					System.err.println("Cannot parse " + tokens.get(1) + " to numeric.");
					step= 0;
				} 
			}			
			from += step; 
			to += step;
			if(gc.getSamSeqDict() != null && !gc.getSamSeqDict().isEmpty()){
				int chromLen= gc.getSamSeqDict().getSequence(chrom).getSequenceLength();
				if(to > chromLen){
					to= chromLen;
					from= to - gc.getUserWindowSize() + 1;
				}
			}			
			return chrom + ":" + from + "-" + to;
			
		} else if(tokens.get(0).equals("b")){
			int step= (int)Math.rint(windowSize / 10d);
			if(tokens.size() > 1){
				try{
					step= (int)Math.rint(windowSize * Double.parseDouble(tokens.get(1)));
				} catch(NumberFormatException e){
					System.err.println("Cannot parse " + tokens.get(1) + " to numeric.");
					step= 0;
				} 
			}						
			from -= step;
			to -= step; 
			if(from < 1){
				from= 1;
				to= from + gc.getUserWindowSize() - 1;
			}
			return chrom + ":" + from + "-" + to;
		
		} else if(tokens.get(0).startsWith("+") 
				|| tokens.get(0).startsWith("-")){
			int offset= parseStringToIntWithUnits(tokens.get(0));
			from += offset;
			if(from <= 0){
				from= 1;
				to= gc.getGenomicWindowSize();
			} else {
				to += offset;
			}
			return chrom + ":" + from + "-" + to;
		}else if (tokens.equals("q")) {
			System.exit(0);	
		} else {
			throw new RuntimeException("Invalid input for " + tokens);
		}
		return region;
	}
	
	/**Parse the rawInput string in the form ':123-456' to return either
	 * the first int or both ints. 
	 * */
	protected static String parseGoToRegion(String rawInput){
		
		String[] fromToRaw= rawInput.trim().
				replaceAll(",", "").  // Remove thousands sep if any
				replace(GenomicCoords.TICKED, " ").  // From chrom ideogram 
				replaceAll("-", " "). // Replace - with space for splitting 
				split(" +");
		
		List<String> fromTo= new ArrayList<String>();
		for(String x : fromToRaw){
			if(!x.trim().isEmpty()){	
				fromTo.add(x);
			}
		}
		
		if(fromTo.size() == 1){
			return String.valueOf(Utils.parseStringToIntWithUnits(fromTo.get(0)));
			// Integer.parseInt(fromTo[0].trim()); // Check you actually got an int.
			// return fromTo[0].trim();
		} else {
			String xfrom= String.valueOf(Utils.parseStringToIntWithUnits(fromTo.get(0)));
			String xto= String.valueOf(Utils.parseStringToIntWithUnits(fromTo.get(fromTo.size()-1)));
			return xfrom + "-" + xto;
			// Integer.parseInt(fromTo[0].trim()); // Check you actually got an int.
			// Integer.parseInt(fromTo[fromTo.length - 1].trim());
			// return fromTo[0].trim() + "-" + fromTo[fromTo.length - 1].trim();
		}
	}
	
	public static boolean isInteger(String s) {
	 
		s= s.replaceFirst("\\+", ""); // Allow first char to be +
		
		try { 
	        Integer.parseInt(s); 
	    } catch(NumberFormatException e) { 
	        return false; 
	    } catch(NullPointerException e) {
	        return false;
	    }
	    // only got here if we didn't return false
	    return true;
	}

	/**
	 * Average of ints in array x. Adapted from:
	 * http://stackoverflow.com/questions/10791568/calculating-average-of-an-array-list
	 * null values are ignored, like R mean(..., na.rm= TRUE). 
	 * Returns Float.NaN if input list is empty or only nulls. You can check for Float.NaN
	 * with Float.isNaN(x); 
	 * @param marks
	 * @return
	 */
	public static Double calculateAverage(List<Double> list) {
		double sum = 0;
		long  N= 0;  
		if(!list.isEmpty()) {
			for (Double z : list) {
				if(z  != null && !Double.isNaN(z)){
					sum += z;
					N++;
				}
			}
			return (double)sum / N;
		}
		return Double.NaN;
	}

	/**
	 * Binary search to get the index position of the value in list closest to a given value.
	 * The searched list is expected to be sorted, there is no check whether this is the case.
	 * @param genomePos
	 * @param mapping
	 * @return
	 */
	public static int getIndexOfclosestValue(double genomePos, List<Double> mapping){

		// Position is before or after the extremes of the list.
		if(genomePos <= mapping.get(0)){
			return 0;
		}
		if(genomePos >= mapping.get(mapping.size()-1)){
			return mapping.size()-1;
		}
		
		int closest= Arrays.binarySearch(mapping.toArray(new Double[mapping.size()]), (double)genomePos);
		if(closest < 0){
			// If < 0 the value is not found in the list and binarySearch returns the insertion point.
			// See binarySearch docs. We need to convert the insertion point to the closest value.
			int insertionPoint= -(closest + 1);
			double leftDiff= genomePos - mapping.get(insertionPoint - 1);
			double rightDiff= mapping.get(insertionPoint) - genomePos;
			if(leftDiff < rightDiff){
				return insertionPoint - 1;
			} else {
				return insertionPoint;
			}
		}
		return closest;
		
//		double bestDiff= Integer.MAX_VALUE;
//		int closest= -1;
//		for(int idx= 0; idx < mapping.size(); idx++){ 
//			// Iterate through entire list to find closest position on screen, it's a bit wasteful since
//			// the list is ordered, but it's ok.
//			double candidate= mapping.get(idx);
//			double diff= Math.abs(genomePos - candidate);
//			if(diff < bestDiff){
//				closest= idx;
//				bestDiff= diff;
//			}
//		}
//		if(closest < 0){
//			throw new RuntimeException("Invalid index position: " + closest);
//		}
//		return closest;
	}
	
	/** Return true  */
	public static boolean allIsNaN(Iterable<Double> x){
		for(Double y : x){
			if(!y.isNaN()){
				return false;
			}
		}
		return true;
	}

	
	public static List<Double> seqFromToLenOut(double from, double to, int lengthOut){
		
		if(lengthOut < 0){
			String msg= "Sequence from " + from + " to " + to + ". Invalid lenght of sequence: Cannot be < 1. Got " + lengthOut;
			throw new RuntimeException(msg);
		} 
		if(lengthOut == 0){ // Consistent with R seq(): return 0-length vector
			return new ArrayList<Double>();
		}		
		
		List<Double> mapping= new ArrayList<Double>();
		double span= to - from + 1;
		double step= ((double)span - 1)/(lengthOut - 1);
		mapping.add((double)from);
		for(int i= 1; i < lengthOut; i++){
			mapping.add((double)mapping.get(i-1)+step);
		}

		if(lengthOut == 1){ // Consistent with R seq(from, to, length.out= 1) -> from
		//	mapping.add((double)from);
			return mapping;
		}
		
		double diffTo= Math.abs(mapping.get(mapping.size() - 1) - to);
		if(diffTo > Math.abs((to + 1e-9))){
			String msg= "Error generating sequence from " + from + " to " + to + " length " + lengthOut + "\n" +
					     "Last point: " + mapping.get(mapping.size() - 1) + "\n" +
					     "To diff: " + diffTo + "\n" +
					     "Step: " + step;
			throw new RuntimeException(msg);
		} else {
			mapping.set(mapping.size()-1, (double)to);
		}
		
		double diffFrom= Math.abs(mapping.get(0) - from);		
		if(diffFrom > 0.01 || mapping.size() != lengthOut){
			String msg= "Error generating sequence:\n" +
					    "Expected size: " + lengthOut + "; Effective: " + mapping.size() + "\n" + 
					    "From diff: " + diffFrom;
			throw new RuntimeException(msg);
		}
		return mapping;
	}

	/*Nicely tabulate list of rows. Each row is tab separated **/
	public static List<String> tabulateList(List<String> rawList) {
		
		// * Split each row in a list of strings. I.e. make list of lists
		List<ArrayList<String>> rawTable= new ArrayList<ArrayList<String>>();
		int ncol= 0;
		for(String x : rawList){
			List<String> row = new ArrayList<String>();
			for(String item : x.split("\t")){
				row.add(item);
			}
			rawTable.add((ArrayList<String>) row);
			// * Get max number of columns (max list length)
			if(row.size() > ncol){
				ncol= row.size(); 
			}
		}
				
		// * Iterate through each column 
		List<ArrayList<String>> paddedTable= new ArrayList<ArrayList<String>>();
		
		for(int i= 0; i < ncol; i++){
			// Collect all items in column i in a list. I.e. select column i
			List<String> col= new ArrayList<String>();
			for(ArrayList<String> row : rawTable){
				if(row.size() > i){
					col.add(row.get(i));
				} else { // If a row doesn't have enough columns add a dummy field 
					col.add("");
				}
			}
			// Get the longest string in this column
			int maxStr= 0;
			for(String x : col){
				if(x.length() > maxStr){
					maxStr= x.length();
				}
			}
			// ** Pass thorugh the column again and pad with spaces to match length of longest string
			// maxStr+=1; // +1 is for separating
			for(int j= 0; j < col.size(); j++){
				String padded= String.format("%-" + maxStr + "s", col.get(j));
				col.set(j, padded);
			}
			paddedTable.add((ArrayList<String>) col);
		}
		// Each list in padded table is a column. We need to create rows as strings
		List<String> outputTable= new ArrayList<String>();
		for(int r= 0; r < rawList.size(); r++){
			StringBuilder row= new StringBuilder();
			for(int c= 0; c < paddedTable.size(); c++){
				row.append(paddedTable.get(c).get(r) + " ");
			}
			outputTable.add(row.toString().trim());
		}
		return outputTable;
	}

	/** Function to round x and y to a number of digits enough to show the difference in range
	 * This is for pretty printing only.
	 * */
	public static Double[] roundToSignificantDigits(double x, double y, int nSignif) {

		Double[] rounded= new Double[2];
		
	    double diff= Math.abs(x - y);
	    if (diff < 1e-16){
	    	rounded[0]= x;
	    	rounded[1]= y;
	    	return rounded;
	    }
	    if(diff > 1){
	    	// Round to 2 digits regardless of how large is the diff
	    	rounded[0]= Math.rint(x * Math.pow(10.0, nSignif))/Math.pow(10.0, nSignif);
	    	rounded[1]= Math.rint(y * Math.pow(10.0, nSignif))/Math.pow(10.0, nSignif);
			return rounded;
	    } else {
	    	// Diff is small, < 1. See how many digits you need to approximate
	    	// Get number of leading zeros
	    	int nzeros= (int) (Math.ceil(Math.abs(Math.log10(diff))) + nSignif);
	    	rounded[0]= Math.rint(x * Math.pow(10.0, nzeros))/Math.pow(10.0, nzeros);
	    	rounded[1]= Math.rint(y * Math.pow(10.0, nzeros))/Math.pow(10.0, nzeros);
	    	return rounded;
	    }   	        		
	}
	
	/** From http://stackoverflow.com/questions/202302/rounding-to-an-arbitrary-number-of-significant-digits */
	public static double roundToSignificantFigures(double num, int n) {
	    if(num == 0) {
	        return 0;
	    }

	    final double d = Math.ceil(Math.log10(num < 0 ? -num: num));
	    final int power = n - (int) d;

	    final double magnitude = Math.pow(10, power);
	    final long shifted = Math.round(num*magnitude);
	    return shifted/magnitude;
	}
	
	/** Convert 000 and 000,000 to k and M suffix. E.g. 1000 -> "1k"; 123000000 -> 123M
	 * See also roundToSignificantFigures() to round numbers.  
	 * */
	public static String parseIntToMetricSuffix(int x){
		String xint= String.valueOf(x);
		if(xint.endsWith("000000")){
			xint= xint.replaceAll("000000$", "M");
		} else if(xint.endsWith("000")){
			xint= xint.replaceAll("000$", "k");
		}
		return xint;
	}

	
	/** Returns true if URL file exists. 
	 * NB: Returns true also if the url path exists but it's just a directory and not a file! 
	 * From http://stackoverflow.com/questions/4596447/check-if-file-exists-on-remote-server-using-its-url
	 * */
	public static boolean urlFileExists(String URLName){

		try{ // For ftp files
			InputStream ftp = new URL(URLName).openStream();
			ftp.close();
			return true;
		} catch(Exception e){
			//
		}
		
	    try {
	        HttpURLConnection.setFollowRedirects(false);
	        // note : you may also need
	        //        HttpURLConnection.setInstanceFollowRedirects(false)
	        HttpURLConnection con =
	           (HttpURLConnection) new URL(URLName).openConnection();
	        con.setRequestMethod("HEAD");
	        return (con.getResponseCode() == HttpURLConnection.HTTP_OK);
	      }
	      catch (Exception e) {
	         return false;
	      }
	}

//	public static List<String> checkAvailableInput(List<String> inputFileList){
//		
//	}
	
	/** Add track(s) to list of input files 
	 * @param inputFileList Existing list of files to be extended
	 * @param newFileNames List of files to append
	 * @throws InvalidCommandLineException 
	 */
	public static void addSourceName(List<String> inputFileList, List<String> newFileNames, int debug) throws InvalidCommandLineException {

		List<String> dropMe= new ArrayList<String>();
		List<String> addMe= new ArrayList<String>();
		for(String x : newFileNames){
			x= x.trim();
			if(!new File(x).isFile() && !Utils.urlFileExists(x) && !Utils.isUcscGenePredSource(x)){
				dropMe.add(x);
				System.err.println("Unable to add " + x);
				if(debug == 2){
					throw new InvalidCommandLineException();
				}
			} 
		}
		for(String x : dropMe){
			//System.err.println("\nWarning: File " + x + " is not a local file.\n");
			newFileNames.remove(x);
		}
		for(String x : addMe){
			newFileNames.add(x);
		}
		inputFileList.addAll(newFileNames);
		
	}

	public static String printSamSeqDict(SAMSequenceDictionary samSeqDict, int graphSize){
		
		if(samSeqDict == null || samSeqDict.isEmpty()){
			System.err.println("Sequence dictionary not available.");
			return "";
		}
		
		// Prepare a list of strings. Each string is a row tab separated
		List<String> tabList= new ArrayList<String>();
		int maxChromLen= -1;
		for(SAMSequenceRecord x : samSeqDict.getSequences()){
			String row= x.getSequenceName() + "\t" + x.getSequenceLength();			
			tabList.add(row);
			if(x.getSequenceLength() > maxChromLen){
				maxChromLen= x.getSequenceLength(); 
			}
		}
		double bpPerChar= (double)maxChromLen / (double)graphSize;
		for(int i= 0; i < samSeqDict.getSequences().size(); i++){
			SAMSequenceRecord x= samSeqDict.getSequences().get(i);
			int n= (int)Math.rint(x.getSequenceLength()/bpPerChar);
			String bar= StringUtils.join(Collections.nCopies(n, "|"), ""); // String.join("", Collections.nCopies(n, "|"));
			String row= tabList.get(i) + "\t" + bar;
			tabList.set(i, row);
		}
		List<String> table= Utils.tabulateList(tabList);
		StringBuilder out= new StringBuilder();
		for(String x : table){
			out.append(x).append("\n");
		}
		return out.toString().trim();
	}
	
	/** Parse cmdInput to extract the integer after the arg. (see tests)
	 * @param defaultInt Default value if parsing returns nonsense
	 * @throws InvalidCommandLineException 
	 * */
	public static int parseZoom(String cmdInput, int defaultInt) throws InvalidCommandLineException {
		String[] zz= cmdInput.trim().split(" +");
		int nz= defaultInt;
		if(zz.length >= 2){
			nz= Integer.parseInt(zz[1]);
		} 
		if(nz < 0){
			throw new InvalidCommandLineException();
		}
		return nz;
	}

	/** Split string x in tokens. Effectively just a friendly wrapper around StrTokenizer.
	 * Use *single* quotes for avoiding splitting. 
	 */
	public static ArrayList<String> tokenize(String x, String delimiterString){
		
		if(x == null){
			return null;
		}
		
		// See also http://stackoverflow.com/questions/38161437/inconsistent-behaviour-of-strtokenizer-to-split-string
		StrTokenizer str= new StrTokenizer(x);
    	str.setTrimmerMatcher(StrMatcher.spaceMatcher());
		str.setDelimiterString(delimiterString);
		str.setQuoteChar('\'');
		ArrayList<String> tokens= (ArrayList<String>) str.getTokenList();
		for(int i= 0; i < tokens.size(); i++){
			String tok= tokens.get(i).trim();
			tokens.set(i, tok);
		}
		return tokens;
	
	}
	
	public static String stripAnsiCodes(String x){
		return x.replaceAll("\\033\\[[;\\d]*m", "");
	}
	
	/** Get a filaname to write to. GenomicCoords obj is used to get current position and 
	 * create a suitable filename from it, provided a filename is not given.
	 * The string '%r' in the file name, is replaced with the current position. Useful to construct
	 * file names like myPeaks.%r.pdf -> myPeaks.chr1_1234-5000.pdf.
	 * */
	public static String parseCmdinputToGetSnapshotFile(String cmdInput, GenomicCoords gc) throws IOException{
		
		final String REGVAR= "%r";
		
		String snapshotFile= cmdInput.trim().replaceAll("^save", "").trim();
		
		String region= gc.getChrom() + "_" + gc.getFrom() + "-" + gc.getTo();
		
		if(snapshotFile.isEmpty()){
			snapshotFile= REGVAR + ".txt"; 
		} else if(snapshotFile.equals(".pdf")){
			snapshotFile= REGVAR + ".pdf";
		} 
		snapshotFile= snapshotFile.replace(REGVAR, region); // Special string '%r' is replaced with the region 
		snapshotFile= Utils.tildeToHomeDir(snapshotFile);
		
		File file = new File(snapshotFile);
		if(file.exists() && !file.canWrite()){
			System.err.println("Cannot write to " + snapshotFile);
			snapshotFile= null;
			return snapshotFile;
		}
		if(!file.exists()){
			try{
				file.createNewFile();
				file.delete();
			} catch(IOException e) {
				System.err.println("Cannot create file " + snapshotFile);
				snapshotFile= null;
				return snapshotFile;
			}
		}
		return snapshotFile;
	}
	
	/** Expand ~/ to user's home dir in file path name. See tests for behaviour
	 * */
	public static String tildeToHomeDir(String path){
		return path.replaceAll("^~" + File.separator, System.getProperty("user.home") + File.separator);
	}
	
	/**
	 * Count reads in interval using the given filters.
	 * @param bam
	 * @param gc
	 * @param filters List of filters to apply
	 * @return
	 * @throws MalformedURLException 
	 */
	public static long countReadsInWindow(String bam, GenomicCoords gc, List<SamRecordFilter> filters) throws MalformedURLException {

		/*  ------------------------------------------------------ */
		/* This chunk prepares SamReader from local bam or URL bam */
		UrlValidator urlValidator = new UrlValidator();
		SamReaderFactory srf=SamReaderFactory.make();
		srf.validationStringency(ValidationStringency.SILENT);
		SamReader samReader;
		if(urlValidator.isValid(bam)){
			samReader = srf.open(SamInputResource.of(new URL(bam)).index(new URL(bam + ".bai")));
		} else {
			samReader= srf.open(new File(bam));
		}
		/*  ------------------------------------------------------ */
		
		long cnt= 0;
		
		Iterator<SAMRecord> sam= samReader.query(gc.getChrom(), gc.getFrom(), gc.getTo(), false);
		AggregateFilter aggregateFilter= new AggregateFilter(filters);
		while(sam.hasNext()){
			SAMRecord rec= sam.next();
			if( !aggregateFilter.filterOut(rec) ){
				cnt++;
			}
		}
		try {
			samReader.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return cnt;
	}

	/** Sort Map by value in descending order. See
	 * http://stackoverflow.com/questions/109383/sort-a-mapkey-value-by-values-java
	 *  */
	public static <K, V extends Comparable<? super V>> Map<K, V> sortByValue( Map<K, V> map ){
	    List<Map.Entry<K, V>> list = new LinkedList<>( map.entrySet() );
	    Collections.sort( list, new Comparator<Map.Entry<K, V>>() {
	        @Override
	        public int compare( Map.Entry<K, V> o1, Map.Entry<K, V> o2 ){
	            return -( o1.getValue() ).compareTo( o2.getValue() );
	        }
	    });
	    
	    Map<K, V> result = new LinkedHashMap<>();
	    for (Map.Entry<K, V> entry : list){
	        result.put( entry.getKey(), entry.getValue() );
	    }
	    return result;
	}

//	public static List<SamRecordFilter> cleanInappropriateCallIfNotPairedRead(List<SamRecordFilter> filter){
//		List<SamRecordFilter> cleanfilter= new ArrayList<SamRecordFilter>(); 
//		for(SamRecordFilter x : filter){
//			if(x.equals(new FirstOfPairFilter(true)) ||
//			   x.equals(new FirstOfPairFilter(false))){
//			   //
//			} else {
//				cleanfilter.add(x);
//			}
//		}
//		return cleanfilter;
//	} 
	
	public static TabixFormat trackFormatToTabixFormat(TrackFormat fmt){
		
		TabixFormat tbx= null;
		if(fmt.equals(TrackFormat.BAM)){
			tbx= TabixFormat.SAM; 
		} else if (fmt.equals(TrackFormat.BED) || fmt.equals(TrackFormat.BEDGRAPH)){
			tbx= TabixFormat.BED; 
		} else if (fmt.equals(TrackFormat.GFF) || fmt.equals(TrackFormat.GTF)){
			tbx= TabixFormat.GFF;
		} else if (fmt.equals(TrackFormat.VCF)){
			tbx= TabixFormat.VCF;
		} else {
			throw new RuntimeException();
		}
		return tbx;
		
	}
	
	/** Same as R range(..., na.rm= TRUE) function: Return min and max of 
	 * list of values ignoring NaNs.
	 * */
	public static Double[] range(List<Double> y){
		Double[] range= new Double[2];
		
		Double ymin= Double.NaN;
		Double ymax= Double.NaN;
		for(Double x : y){
			if(!x.isNaN()){
				if(x > ymax || ymax.isNaN()){
					ymax= x;
				} 
				if(x < ymin || ymin.isNaN()){
					ymin= x;
				}
			}
		}
		range[0]= ymin;
		range[1]= ymax;
		return range;
	}
	
	/** Convert coordinates to string suitable for initializing GenomicCoords object.
	 * See tests for handling odd cases.
	 * */
	public static String coordinatesToString(String chrom, Integer from, Integer to){
		
		if(from == null){
			from = -1;
		}
		if(to == null){
			to = -1;
		}
		
		String xfrom;
		if(from <= 0){
			xfrom= "1";
		} else {
			xfrom= from.toString();
		}
		String xto;
		if(to <= 0 || to < from){
			xto= "";
		} else {
			xto= "-" + to.toString();
		}
		return chrom + ":" + xfrom + xto;
	} 
	
	/** True if filename is a UCSC genePred file or a valid connection to database. 
	 * */
	public static boolean isUcscGenePredSource(String filename) {
		try{
			new UcscGenePred(filename, 1000);
			return true;
		} catch(Exception e) {
			return false;
		}
	}

	/** Return list of files matching the list of glob expressions. This method should behave
	 * similarly to GNU `ls` command. 
	 * 
	 * * Directories are not returned 
	 * 
	 * Files are returned as they are matched, 
	 * so not necessarily in alphabetical order.
	 *  
	 * */
	public static List<String> globFiles(List<String> cmdInput) throws IOException {
		
		List<String> globbed= new ArrayList<String>();
		
		for(String x : cmdInput){
			
			if(Utils.urlFileExists(x)){
				globbed.add(x);
				continue;
			}
			x= Utils.tildeToHomeDir(x);
			x= x.replaceAll(File.separator + "+$", ""); // Remove trailing dir sep
			x= x.replaceAll(File.separator + "+", File.separator); // Remove double dir sep like "/foo//bar" -> /foo/bar 

			String location;
			if(new File(x).isDirectory()){
				location= x;
				x= x + File.separator + "*"; // From "my_dir" to "my_dir/*" 
			} else {
				location= new File(x).getParent();
				if(location == null){
					location= "";
				}
			} 
			for(Path p : match(x, location)){
				globbed.add(p.toString());
			}
		}		
		return globbed;
	}

	/** Search for a glob pattern in a given directory and its sub directories
	 * See http://javapapers.com/java/glob-with-java-nio/
	 * */
	private static List<Path> match(String glob, String location) throws IOException {
		
		final List<Path> globbed= new ArrayList<Path>();
		
		final PathMatcher pathMatcher = FileSystems.getDefault().getPathMatcher("glob:" + glob);
		
		Files.walkFileTree(Paths.get(location), new SimpleFileVisitor<Path>() {
			
			@Override
			public FileVisitResult visitFile(Path path, BasicFileAttributes attrs) throws IOException {
				if (pathMatcher.matches(path)) {
					globbed.add(path);
				}
				return FileVisitResult.CONTINUE;
			}

			@Override
			public FileVisitResult visitFileFailed(Path file, IOException exc)
					throws IOException {
				return FileVisitResult.CONTINUE;
			}
		});
		return globbed;
	}

	/** Query github repo to check if a version newer then this one is available.
	 * Returns list of length 2: ["this version", "latest version on github"]
	 * @param timeout Return if no response is received after so many milliseconds.
	 * @throws IOException 
	 * */
	protected static List<String> checkUpdates(long timeout) throws IOException {
		
		List<String> thisAndGitVersion= new ArrayList<String>();
		
		// Get version of this ASCIIGenome
		thisAndGitVersion.add(ArgParse.VERSION);
		
		BufferedReader br= null;
		timeout= timeout + System.currentTimeMillis();
		// Get github versions
		URL url = new URL("https://api.github.com/repos/dariober/ASCIIGenome/tags");
		br = new BufferedReader(new InputStreamReader(url.openStream()));
//		if(br == null){
//			System.err.println("Note: Couldn't check for updates.");
//			thisAndGitVersion.add(ArgParse.VERSION); // If timed out assume up to date.
//		}

        String line;
        StringBuilder sb= new StringBuilder();
        while ((line = br.readLine()) != null) {
        	sb.append(line + '\n');
        }

        JsonElement jelement = new JsonParser().parse(sb.toString());
        JsonArray  jarr = jelement.getAsJsonArray(); //.getAsJsonObject();
        JsonObject jobj = jarr.get(0).getAsJsonObject();
        String tag= jobj.get("name").getAsString().replaceFirst("v", "");
        
        thisAndGitVersion.add(tag);
        return thisAndGitVersion;
	}
	
	/**
	 * Compares two version strings. 
	 * 
	 * Use this instead of String.compareTo() for a non-lexicographical 
	 * comparison that works for version strings. e.g. "1.10".compareTo("1.6").
	 * 
	 * From 
	 * http://stackoverflow.com/questions/6701948/efficient-way-to-compare-version-strings-in-java
	 * 
	 * @note It does not work if "1.10" is supposed to be equal to "1.10.0".
	 * 
	 * @param str1 a string of ordinal numbers separated by decimal points. 
	 * @param str2 a string of ordinal numbers separated by decimal points.
	 * @return The result is a negative integer if str1 is _numerically_ less than str2. 
	 *         The result is a positive integer if str1 is _numerically_ greater than str2. 
	 *         The result is zero if the strings are _numerically_ equal.
	 *        
	 */
	public static int versionCompare(String str1, String str2) {
	    String[] vals1 = str1.split("\\.");
	    String[] vals2 = str2.split("\\.");
	    int i = 0;
	    // set index to first non-equal ordinal or length of shortest version string
	    while (i < vals1.length && i < vals2.length && vals1[i].equals(vals2[i])) {
	      i++;
	    }
	    // compare first non-equal ordinal number
	    if (i < vals1.length && i < vals2.length) {
	        int diff = Integer.valueOf(vals1[i]).compareTo(Integer.valueOf(vals2[i]));
	        return Integer.signum(diff);
	    }
	    // the strings are equal or one string is a substring of the other
	    // e.g. "1.2.3" = "1.2.3" or "1.2.3" < "1.2.3.4"
	    return Integer.signum(vals1.length - vals2.length);
	}

	/** Stream the raw line through awk and return true if the output of awk is the same as the input
	 * line. If awk returns empty output then the line didn't pass the awk filter.
	 * If output is not empty and not equal to input, return null.
	 * See tests for behaviour. 
	 * */
	public static Boolean passAwkFilter(String rawLine, String awkScript) throws IOException {
		
		awkScript= awkScript.trim();

		if(awkScript.isEmpty()){
			return true;
		}
		
		// We need to separate the awk script from the arguments. The arguments could contain single quotes:
		// -v var=foo '$1 == var && $2 == baz'
		// For this, reverse the string and look for the first occurrence of "' " which corresponds to 
		// the opening of the awk script.
		int scriptStartsAt= awkScript.length() - new StringBuilder(awkScript).reverse().toString().indexOf("' ");
		if(scriptStartsAt == awkScript.length() + 1){
			scriptStartsAt= 1;
		}
		// Now add the functions to the script right at the start of the script, after the command args
		awkScript= awkScript.substring(0, scriptStartsAt) + Track.awkFunc + awkScript.substring(scriptStartsAt);
				
		InputStream is= new ByteArrayInputStream(rawLine.getBytes(StandardCharsets.US_ASCII));
		
		String[] args= Utils.tokenize(awkScript, " ").toArray(new String[0]); 
		
		PrintStream stdout = System.out;
		ByteArrayOutputStream baos = new ByteArrayOutputStream();
		try{
			PrintStream os= new PrintStream(baos);
			new org.jawk.Main(args, is, os, System.err);
		} catch(Exception e){
			throw new IOException();
		} finally{
			System.setOut(stdout);
			is.close();
		}
		String output = new String(baos.toByteArray(), StandardCharsets.UTF_8);
		if(output.trim().isEmpty()){
			return false;
		} else if(output.trim().equals(rawLine.trim())){
			return true;
		} else {
			// Awk output is not empty or equal to input line. Reset awk script and return null
			// to signal this condition.
			return null;
		}
		
	}

	/** Stream the raw line through awk and return true if the output of awk is the same as the input
	 * line. If awk returns empty output then the line didn't pass the awk filter.
	 * If output is not empty and not equal to input, return null.
	 * See tests for behaviour. 
	 * */
	public static boolean[] passAwkFilter(String[] rawLines, String awkScript) throws IOException {
		
		boolean[] results= new boolean[rawLines.length];

		awkScript= awkScript.trim();

		if(awkScript.isEmpty()){
			for(int i= 0; i < rawLines.length; i++){
				results[i]= true;
			}
			return results;
		}
		
		// We need to separate the awk script from the arguments. The arguments could contain single quotes:
		// -v var=foo '$1 == var && $2 == baz'
		// For this, reverse the string and look for the first occurrence of "' " which corresponds to 
		// the opening of the awk script.
		int scriptStartsAt= awkScript.length() - new StringBuilder(awkScript).reverse().toString().indexOf("' ");
		if(scriptStartsAt == awkScript.length() + 1){
			scriptStartsAt= 1;
		}
		// Now add the functions to the script right at the start of the script, after the command args
		awkScript= awkScript.substring(0, scriptStartsAt) + Track.awkFunc + awkScript.substring(scriptStartsAt);
		
		ByteArrayOutputStream baosIn = new ByteArrayOutputStream();
		for (String line : rawLines) {
			baosIn.write((line+"\n").getBytes());
		}
		InputStream is= new ByteArrayInputStream(baosIn.toByteArray());
		
		String[] args= Utils.tokenize(awkScript, " ").toArray(new String[0]); 
		
		PrintStream stdout = System.out;
		ByteArrayOutputStream baos = new ByteArrayOutputStream();
		try{
			PrintStream os= new PrintStream(baos);
			new org.jawk.Main(args, is, os, System.err);
		} catch(Exception e){
			throw new IOException();
		} finally{
			System.setOut(stdout);
			is.close();
		}
		
		String output[] = new String(baos.toByteArray(), StandardCharsets.US_ASCII).split("\n");
		int j= 0;
		for(int i=0; i < rawLines.length; i++){
			String inLine= rawLines[i];
			if(output.length > j){
				String outLine= output[j];				
				if(inLine.equals(outLine)){
					results[i]= true;
					j++;
				} else {
					results[i]= false;
				}
			} else {
				results[i]= false;
			}
		}
		return results;
	}

	
	/** Right-pad each line in string x with whitespaces. Each line is 
	 * defined by the newline. Each line is padded to become at leas of length size.   
	 * */
	public static String padEndMultiLine(String x, int size) {

		if(x.isEmpty()){
			return x;
		}
		
		List<String> split= Splitter.on("\n").splitToList(x);
		List<String> mline= new ArrayList<String>(split);
		
		for(int i= 0; i < mline.size(); i++){
			mline.set(i, Strings.padEnd(mline.get(i), size, ' '));
		}
		return Joiner.on("\n").join(mline);
	}

	/**Parse string region in the form <chrom>:[start[-end]] to a list containing the
	 * the three elements. See tests for behavior.
	 * <start> and <end> if present are guaranteed to be parsable to positive int.
	 * @throws InvalidGenomicCoordsException 
	 * */
	public static List<String> parseStringCoordsToList(String region) throws InvalidGenomicCoordsException {
		
		List<String> coords= new ArrayList<String>(3);
		coords.add(null);
		coords.add("1");
		coords.add("536870912"); // Max size of binning index for tabix is 2^29. See also https://github.com/samtools/htslib/issues/435

		String chrom= StringUtils.substringBeforeLast(region, ":").trim();
		coords.set(0, chrom);
		
		String fromTo= StringUtils.substringAfterLast(region, ":").replaceAll(",", "").replaceAll("\\s", "");
		if(fromTo.isEmpty()){
			// Only chrom given
			return coords;
		}
		
		if( ! fromTo.replaceFirst("-", "").matches("[0-9]+")){
			// If the from-to part does not contain only digits with the exception of the - separator,
			// we assume this is a chrom name containing : and missing the from-to part.
			coords.set(0, region);
			return coords;
		}

		int nsep= StringUtils.countMatches(fromTo, "-");
		Integer from= null;
		Integer to= null;
		if(nsep == 0){ 
			// Only start position given
			from= Integer.parseInt(StringUtils.substringBefore(fromTo, "-").trim());
			to= from;
		} else if(nsep == 1){ // From and To positions given.
			from= Integer.parseInt(StringUtils.substringBefore(fromTo, "-").trim());
			to= Integer.parseInt(StringUtils.substringAfter(fromTo, "-").trim());
			if(from > to || from <= 0 || to <= 0 || (to-from+1) > 536870912){
				throw new InvalidGenomicCoordsException();	
			}
		} else {
			throw new InvalidGenomicCoordsException();
		}
		coords.set(1, from.toString());
		coords.set(2, to.toString());

		return coords;
	}

	/** Prepare SamReader from local bam or URL bam */
	public static SamReader getSamReader(String workFilename) throws MalformedURLException {
		UrlValidator urlValidator = new UrlValidator();
		SamReaderFactory srf=SamReaderFactory.make();
		srf.validationStringency(ValidationStringency.SILENT);
		SamReader samReader;
		if(urlValidator.isValid(workFilename)){
			samReader = srf.open(SamInputResource.of(new URL(workFilename)).index(new URL(workFilename + ".bai")));
		} else {
			samReader= srf.open(new File(workFilename));
		}
		return samReader;
	}

	/**Return the indexes of the printable characters interspersed in ansi formatting  
	 * MEMO: The sequence \033 is NOT four characters. It's just one!
	 * */
	public static List<Integer> indexOfCharsOnFormattedLine(String fmtString) {

		List<Integer> idx= new ArrayList<Integer>();
		// Start assuming that each char in input string is a printable char, not part of ANSI code.
		for(int i= 0; i < fmtString.length(); i++){
			idx.add(i);
		}
		
		Matcher m= Pattern.compile("\\033\\[[;\\d]*m").matcher(fmtString);
		while(m.find()){
			// Now remove from the index list the characters that are part of the escape
			for(int i= m.start(); i < m.end(); i++){
				int del= idx.indexOf(i);
				idx.remove(del);
			}
		}
		return idx;
	}
	
}
