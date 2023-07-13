package samTextViewer;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.Reader;
import java.math.BigDecimal;
import java.math.RoundingMode;
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
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.EnumSet;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;

import org.apache.commons.lang3.StringUtils;
import org.apache.commons.lang3.text.StrMatcher;
import org.apache.commons.lang3.text.StrTokenizer;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.apache.commons.validator.routines.UrlValidator;
import org.biojava.nbio.genome.parsers.gff.FeatureI;
import org.biojava.nbio.genome.parsers.gff.FeatureList;
import org.biojava.nbio.genome.parsers.gff.GFF3Reader;
import org.broad.igv.bbfile.BBFileReader;
import org.broad.igv.bbfile.BedFeature;
import org.broad.igv.bbfile.BigBedIterator;
import org.broad.igv.bbfile.BigWigIterator;
import org.broad.igv.bbfile.WigItem;
import org.broad.igv.tdf.TDFReader;

import com.google.common.base.Joiner;
import com.google.common.base.Splitter;
import com.google.common.base.Strings;
import com.google.gson.JsonArray;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.google.gson.JsonParser;

import exceptions.InvalidColourException;
import exceptions.InvalidCommandLineException;
import exceptions.InvalidGenomicCoordsException;
import exceptions.InvalidRecordException;
import faidx.Faidx;
import faidx.UnindexableFastaFileException;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.filter.AggregateFilter;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.index.tabix.TabixFormat;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.tribble.readers.TabixReader;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderVersion;
import tracks.IntervalFeature;
import tracks.Track;
import tracks.TrackFormat;
import utils.Tokenizer;

/**
 * @author berald01
 *
 */
public class Utils {
	
	public static double round(double value, int places) {
	    if (places < 0) throw new IllegalArgumentException();
	 
	    BigDecimal bd = new BigDecimal(Double.toString(value));
	    bd = bd.setScale(places, RoundingMode.HALF_EVEN);
	    return bd.doubleValue();
	}
	
	/**Create temp file in the current working directory or in the system's
	 * tmp dir if failing to create in cwd.*/
	public static File createTempFile(String prefix, String suffix, boolean deleteOnExit){
		File tmp = null;
		try{
			tmp= File.createTempFile(prefix, suffix, new File(System.getProperty("user.dir")));
		} catch(IOException e){
			try {
				tmp= File.createTempFile(prefix, suffix);
			} catch (IOException e1) {
				e1.printStackTrace();
			}
		}
		if(deleteOnExit) {
			tmp.deleteOnExit();	
		}
		return tmp;
	}
	
	/**Return a buffered reader which iterate through the file or URL. Input may be
	 * compressed. 
	 * */
	public static BufferedReader reader(String fileOrUrl) throws MalformedURLException, IOException{
		
		BufferedReader br= null;
		InputStream gzipStream= null;
		UrlValidator urlValidator = new UrlValidator();
		if(fileOrUrl.endsWith(".gz") || fileOrUrl.endsWith(".bgz")){
			if(urlValidator.isValid(fileOrUrl)) {
				gzipStream = new GZIPInputStream(new URL(fileOrUrl).openStream());
			} else {
				gzipStream = new GZIPInputStream(new FileInputStream(fileOrUrl));
			}
			Reader decoder = new InputStreamReader(gzipStream, "UTF-8");
			br = new BufferedReader(decoder);
		} else if(urlValidator.isValid(fileOrUrl)) {
			InputStream instream= new URL(fileOrUrl).openStream();
			Reader decoder = new InputStreamReader(instream, "UTF-8");
			br = new BufferedReader(decoder);
		} else {
			br = new BufferedReader(new FileReader(fileOrUrl));
		}
		return br;
	}
	
	/**Extract template name from sam record read name. I.e. remove
	 * everything after first backspace and possible /1 /2 suffixes.
	 * Avoid using replaceAll() for this as it is much slower.  
	 * */
	public static String templateNameFromSamReadName(String readName){
		if(readName.contains(" ")){
			readName= readName.substring(0, readName.indexOf(" "));
		}
		if(readName.endsWith("/1") || readName.endsWith("/2")){
			readName= readName.substring(0, readName.length() - 2);
		}
		return readName;
	}
	
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

	/** Check if the list of arguments contains "param" and if so return nargs argument after it. 
	 * Returns null if arglist does not contain the parameter.
	 * IMPORTANT SIDE EFFECT: If found, the parameter and its arguments are removed from input list. 
	 * @throws InvalidCommandLineException 
	 * */
	public static List<String> getNArgsForParam(List<String> argList, String param, int nargs) throws InvalidCommandLineException {
		if(nargs < 1){
			System.err.println("narg must be >= 1. Got " + nargs);
			throw new InvalidCommandLineException(); 
		}
//		List<String> args= new ArrayList<String>();
		
		int idx= argList.indexOf(param);
		if(idx == -1){
			return null;
		}
		if(idx == (argList.size() - 1)){
			// If param is the last item in the list you cannot get its argument!
			throw new InvalidCommandLineException();
		}
		if(idx+1+nargs > argList.size()){
			System.err.println("Not enough arguments passed to parameter " + param);
			throw new InvalidCommandLineException();
		}
		
		List<String> args= new ArrayList<String>();
		for(int i= idx+1; i < idx+1+nargs; i++){
			// Do not use List.subList() because you get a view of the input, not a copy.
			args.add(argList.get(i));
		}
		// Remove param and args from list
		argList.remove(idx);
		while(nargs > 0){
			argList.remove(idx);
			nargs--;
		}
		return args;
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
	 * @throws InvalidColourException 
	 * */
	public static List<IntervalFeature> mergeIntervalFeatures(List<IntervalFeature> intervalList, boolean screenCoords) throws InvalidGenomicCoordsException, InvalidColourException{
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
				IntervalFeature x= new IntervalFeature(mergedChrom + "\t" + (mergedFrom-1) + "\t" + mergedTo, TrackFormat.BED, null, -1);
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
		for(IntervalFeature x : mergedList){
			x.getIdeogram(true, true);
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
	@SuppressWarnings("unused")
	public static String initRegionFromFile(String x, String referenceSequence) throws IOException, InvalidGenomicCoordsException, ClassNotFoundException, InvalidCommandLineException, InvalidRecordException, SQLException{
		UrlValidator urlValidator = new UrlValidator();
		String region= "";
		TrackFormat fmt= Utils.getFileTypeFromName(x); 
		
		if(Utils.isCRAM(x) && referenceSequence == null){
		    throw new RuntimeException("Failed to read file " + x + ": CRAM files require a reference genome file");
		}
		
		if(fmt.equals(TrackFormat.BAM)) {
			SamReader samReader;
			if(urlValidator.isValid(x)){
				SamReaderFactory srf = SamReaderFactory.makeDefault();
				if(referenceSequence != null) {
				    srf.referenceSequence(new File(referenceSequence));
				}
				samReader = srf.open(SamInputResource.of(new URL(x)));
			} else {
				SamReaderFactory srf=SamReaderFactory.make();
				if(referenceSequence != null) {
				    srf.referenceSequence(new File(referenceSequence));
				}
				srf.validationStringency(ValidationStringency.SILENT);
				samReader = srf.open(new File(x));
			}
			// Default: Start from the first contig in dictionary
			region= samReader.getFileHeader().getSequence(0).getSequenceName();
			SAMRecordIterator iter = samReader.iterator();
			if(iter.hasNext()){
				// If there are records in this BAM, init from first record
				SAMRecord rec = iter.next();
				if(rec.getContig() != null){
					// See issue#86 for why we need to check null
					region= rec.getContig() + ":" + rec.getAlignmentStart();					
				}
				samReader.close();
			}
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
		
		} else {
			// Input file appears to be a generic interval file. We expect chrom to be in column 1
			// VCF files are also included here since they are either gzip or plain ASCII.
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
				if(fmt.equals(TrackFormat.VCF)){
					region= line.split("\t")[0] + ":" + line.split("\t")[1]; 
				} else {
					IntervalFeature feature= new IntervalFeature(line, fmt, null, -1);
					region= feature.getChrom() + ":" + feature.getFrom();
				}
				br.close();
				return region;
			}
			if(line == null){ // This means the input has no records
				region= "Undefined_contig";
				if(fmt.equals(TrackFormat.VCF)){
					SAMSequenceDictionary seqdict = getVCFHeader(x).getSequenceDictionary();
					if(seqdict != null){
						Iterator<SAMSequenceRecord> iter = seqdict.getSequences().iterator();
						if(iter.hasNext()){
							region= iter.next().getSequenceName();
						}
					}
				}
				return region;
			}
		} 
		System.err.println("Cannot initialize from " + x);
		throw new RuntimeException();
	}

	public static String initRegionFromFile(String x) throws IOException, InvalidGenomicCoordsException, ClassNotFoundException, InvalidCommandLineException, InvalidRecordException, SQLException{	
	    return initRegionFromFile(x, null);
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
			
		} else if(fileName.endsWith(".bam") || 
		          fileName.endsWith(".sam") || 
		          fileName.endsWith(".sam.gz") || 
		          Utils.isCRAM(fileName)){
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
	 * @throws InvalidCommandLineException 
	 */
	public static String parseConsoleInput(List<String> tokens, GenomicCoords gc) throws InvalidGenomicCoordsException, IOException, InvalidCommandLineException{
		
//		String region= "";
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
					throw new InvalidCommandLineException();
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
					throw new InvalidCommandLineException();
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
		} 
//		else if (tokens.get(0).equals("q")) {
//			System.exit(0);	
//		} 
		else {
			throw new RuntimeException("Invalid input for " + tokens);
		}
//		return region;
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

	/** Nicely tabulate list of rows. Each row is tab separated 
	 * The terminalWidth param is used to decide whether rows should be flushed left instead of
	 * being nicely tabulated. If the amount of white space in a cell is too much relative to 
	 * terminalWidth, then flush left. With special value: -1 never flush left, with 0 always flush.
	 */
	public static List<String> tabulateList(List<String> rawList, int terminalWidth, String columnSep) {
		// This method could be streamlined to be more efficient. There are quite a few
		// Lists moved around that could be avoided. However, the size of the table is
		// typically small enough that we prefer a clearer layout over efficiency.
		
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
		List<List<String>> paddedTable= new ArrayList<List<String>>();
		
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
			int maxStr= 1;
			for(String x : col){
				x= Utils.stripAnsiCodes(x);
				if(x.length() > maxStr){
					maxStr= x.length();
				}
			}
			// ** Pass through the column again and pad with spaces to match length of longest string
			for(int j= 0; j < col.size(); j++){
				int nblanks= maxStr - Utils.stripAnsiCodes(col.get(j)).length();
				String padded= col.get(j) + StringUtils.repeat(" ", nblanks); // String.format("%-" + maxStr + "s", col.get(j));
				col.set(j, padded);
			}
			paddedTable.add((List<String>) col);
		}

		// In paddedTable each inner list is a column. Transpose to have
		// inner list as row as it is more convenient from now on.
		List<List<String>> tTable= new ArrayList<List<String>>();
		for(int r= 0; r < rawList.size(); r++){
			List<String> row= new ArrayList<String>();
			for(int c= 0; c < paddedTable.size(); c++){
				row.add(paddedTable.get(c).get(r));
			}
			tTable.add(row);
		}
		
		// If a cell (String) has too many spaces to the right, flush left, i.e. trim(),
		// that cell and all the cells to right on that row. This is to prevent odd
		// formatting where a long string in one cell makes all the other cells look very empty.
		terminalWidth= terminalWidth < 0 ? Integer.MAX_VALUE : terminalWidth;
		for(List<String> row : tTable){
			boolean flush= false;
			for(int i= 0; i < row.size(); i++){
				if(! flush){
					int whiteSize= Utils.stripAnsiCodes(row.get(i)).length() - Utils.stripAnsiCodes(row.get(i)).trim().length();
					if(whiteSize > terminalWidth/4.0){
						// The rest of this row will be flushed
						flush= true;						
					}
				}
				if(flush){
					row.set(i, row.get(i).replaceAll("^ +| +$", ""));
				}
			}
		}

		// Finally, join row into a list of single, printable strings:		
		List<String> outputTable= new ArrayList<String>();
		for(List<String> lstRow : tTable){
			String row= Joiner.on(columnSep).join(lstRow); // Here you decide what separates columns.
			outputTable.add(row.toString()); // Do not use .trim() here otherwise you could strip ANSI formatting
		}
		return outputTable;
	}

	/** Function to round x and y to a number of digits enough to show the difference in range
	 * This is for pretty printing only.
	 * */
	public static String[] roundToSignificantDigits(double x, double y, int nSignif) {

		Double[] rounded= new Double[2];
		
	    double diff= Math.abs(x - y);
	    if (diff < 1e-16){
	    	rounded[0]= x;
	    	rounded[1]= y;
	    }
	    else if(diff > 1){
	    	// Round to 2 digits regardless of how large is the diff
	    	rounded[0]= Math.rint(x * Math.pow(10.0, nSignif))/Math.pow(10.0, nSignif);
	    	rounded[1]= Math.rint(y * Math.pow(10.0, nSignif))/Math.pow(10.0, nSignif);
		} else {
	    	// Diff is small, < 1. See how many digits you need to approximate
	    	// Get number of leading zeros
	    	int nzeros= (int) (Math.ceil(Math.abs(Math.log10(diff))) + nSignif);
	    	rounded[0]= Math.rint(x * Math.pow(10.0, nzeros))/Math.pow(10.0, nzeros);
	    	rounded[1]= Math.rint(y * Math.pow(10.0, nzeros))/Math.pow(10.0, nzeros);
	    }
	    String[] out= new String[2];
	    out[0]= rounded[0].toString();
	    out[1]= rounded[1].toString();
	    return out;
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
	
	/** Add track(s) to list of input files 
	 * @param inputFileList Existing list of files to be extended
	 * @param newFileNames List of files to append
	 * @throws InvalidCommandLineException 
	 * @throws IOException 
	 */
	public static void addSourceName(List<String> inputFileList, List<String> newFileNames, int debug) throws InvalidCommandLineException, IOException {

		List<String> dropMe= new ArrayList<String>();
		List<String> addMe= new ArrayList<String>();
		for(int i= 0; i < newFileNames.size(); i++){
			String x= newFileNames.get(i).trim();
			if(!new File(x).isFile() && !Utils.urlFileExists(x)){
				dropMe.add(x);
				System.err.println("Unable to add " + x);
				throw new InvalidCommandLineException();
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

	/** Read a sample of file x and return its track format. 
	 * This method is not generic so keep it private. Input file must be uncompressed,
	 * must have header if SAM or VCF.
	 * @throws IOException 
	 * */
	public static TrackFormat sniffFile(File x) throws IOException{

		try{
			VCFFileReader vcf= new VCFFileReader(x, false);
			for(@SuppressWarnings("unused") VariantContext rec : vcf){
				//
			}
			vcf.close();
			return TrackFormat.VCF;
		} catch(Exception e){
			//
		}
		
		try{
			SamReader sam = SamReaderFactory.make().open(x);
			for(@SuppressWarnings("unused") SAMRecord rec : sam){
				//
			}
			sam.close();
			return TrackFormat.BAM;
		} catch(Exception e){
			//
		}

		BufferedReader br= new BufferedReader(new FileReader(x));
		int maxLines= 100000;
		boolean firstLine= true;
		String line= null;
		List<String[]> sample= new ArrayList<String[]>();
		while((line= br.readLine()) != null){
			line= line.trim();
			if(line.isEmpty()){
				continue;
			}
			if(firstLine && line.startsWith("##") && line.contains("gff-version")){
				br.close();
				return TrackFormat.GFF;
			}
			firstLine= false;
			if(line.startsWith("#")){
				continue;
			}
			sample.add(line.split("\t"));
			maxLines--;
			if(maxLines == 0){
				br.close();
				break;
			}
		}
		br.close();
		// Try GTF format. We don't distiguish here between GTF and GFF.
		boolean isGtf= true;
		for(String[] s : sample){
			if(s.length < 8){
				isGtf= false;
				break;
			}
			try{
				int start= Integer.valueOf(s[3]);
				int end= Integer.valueOf(s[4]);
				if(start > end || start < 1 || end < 1){
					isGtf= false;
					break;
				}
			} catch(NumberFormatException e){
				isGtf= false;
				break;				
			}
			if( ! (s[6].equals("+") || s[6].equals("-") || s[6].equals("."))){
				isGtf= false;
				break;								
			}
		}
		if(isGtf){
			return TrackFormat.GTF;
		}
		
		// Try bedgraph
		boolean isBedgraph= true;
		for(String[] bdg : sample){
			if(bdg.length < 4){
				isBedgraph= false;
				break;
			}
			
			try{
				Integer.valueOf(bdg[1]);
				Integer.valueOf(bdg[2]);
				Double.valueOf(bdg[3]);
			} catch(NumberFormatException e){
				isBedgraph= false;
				break;				
			}
			if(Integer.valueOf(bdg[1]) < 0 || Integer.valueOf(bdg[2]) < 0 || Integer.valueOf(bdg[1]) > Integer.valueOf(bdg[2])){
				isBedgraph= false;
				break;
			}
		}
		if(isBedgraph){
			return TrackFormat.BEDGRAPH;
		}
		
		// Last option: BED
		boolean isBed= true;
		for(String[] bed : sample){
			if(bed.length < 3){
				isBed= false;
				break;
			}
			int start;
			int end;
			try{
				start= Integer.valueOf(bed[1]);
				end= Integer.valueOf(bed[2]);
			} catch(NumberFormatException e){
				isBed= false;
				break;				
			}
			if(start < 0 || end < 0 || start > end){
				isBed= false;
				break;
			}
		}
		if(isBed){
			return TrackFormat.BED;	
		} else {
			throw new IOException("Input format cannot be determined.");
		}
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
		List<String> table= Utils.tabulateList(tabList, -1, " ");
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
		
		// This is a hack to allow empty tokens to be passed at the command line. 
		// An empty 
		x= x.replace("''", "' '");
		
		// See also http://stackoverflow.com/questions/38161437/inconsistent-behaviour-of-strtokenizer-to-split-string
		StrTokenizer str= new StrTokenizer(x);
    	str.setTrimmerMatcher(StrMatcher.spaceMatcher());
		str.setDelimiterString(delimiterString);
		str.setQuoteChar('\'');
		// str.setIgnoreEmptyTokens(false);
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
		
		String region= gc.getChrom() + "_" + gc.getFrom() + "_" + gc.getTo();
		
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
		return path.replaceAll("^~" + Pattern.quote(File.separator), System.getProperty("user.home") + File.separator);
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
	public static Float[] range(List<Float> y){
		Float[] range= new Float[2];
		
		Float ymin= Float.NaN;
		Float ymax= Float.NaN;
		for(Float x : y){
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
			x= x.replaceAll(Pattern.quote(File.separator) + "+$", ""); // Remove trailing dir sep
			x= x.replaceAll(Pattern.quote(File.separator) + "+", File.separator); // Remove double dir sep like "/foo//bar" -> /foo/bar 

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

		String line;
        StringBuilder sb= new StringBuilder();
        while ((line = br.readLine()) != null) {
        	sb.append(line + '\n');
        }

        JsonElement jelement = new JsonParser().parse(sb.toString());
        JsonArray  jarr = jelement.getAsJsonArray();
        Iterator<JsonElement> iter = jarr.iterator();

        List<String> tag= new ArrayList<String>();
        while(iter.hasNext()){
        	// Get all tags
            JsonObject jobj = (JsonObject) iter.next();
            tag.add(jobj.get("name").getAsString().replaceAll("^[^0-9]*", "").trim());
        }
        if(tag.size() == 0){
        	thisAndGitVersion.add("0");
        } else {
        	thisAndGitVersion.add(tag.get(0));
        }
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

    public static ArrayList<String> execSystemCommand(String[] inputList, List<String> cmd) throws IOException, InterruptedException {

        ProcessBuilder pb = new ProcessBuilder().command(cmd);
        Process p = pb.start();

        ArrayList<String> results= new ArrayList<String>();
        
        Thread readerThread = new Thread(() -> {
            try {
                try (BufferedReader reader = new BufferedReader(new 
                        InputStreamReader(p.getInputStream()))) {
                    String line;
                    while ((line = reader.readLine()) != null) {
                        results.add(line);
                    }
                }
            } catch (Exception e) {
                throw new RuntimeException("Unhandled", e);
            }
        });
        readerThread.start();
        
        OutputStream stdin = p.getOutputStream();

        for(String line : inputList){
            stdin.write(line.getBytes(StandardCharsets.UTF_8));
            stdin.write('\n');
            try {
                stdin.flush();
            } catch (IOException e) {
                throwCmdException(p);
            }
        }

        stdin.close();
        readerThread.join();
        p.waitFor();
        
        if(p.exitValue() != 0) {
            throwCmdException(p);
        }
        return results;
    }
    
    private static void throwCmdException(Process p) throws IOException, InterruptedException {
        p.waitFor();
        System.err.println("Command returned with non-zero exit value " + p.exitValue());
        BufferedReader err = new BufferedReader(new InputStreamReader(p.getErrorStream()));
        String errline = "";
        while ((errline = err.readLine())!= null) {
            System.err.println(errline);
        }
        err.close();
        throw new IOException();  
    }

    /**
     * Convert the awk script to a list suitable for processBuilder
     * @param awkScript e.g.: -v=VAR5 '$2 == 0'
     * @return
     */
    private static ArrayList<String> prepareAwkScript(String awkScript){
        awkScript= awkScript.trim();
        // * Parse command string into arguments and awk script
        ArrayList<String> args= (ArrayList<String>) new Tokenizer(awkScript).tokenize();

        args.add(0, "awk");
        String script= args.remove(args.size()-1); // 'remove' returns the element removed
        script = Track.awkFunc + "\n" + script;
        args.add(script);
        
        return args;
    }
    
	/** Stream the raw line through awk and return true if the output of awk is the same as the input
	 * line. If awk returns empty output then the line didn't pass the awk filter.
	 * If output is not empty and not equal to input, return null.
	 * See tests for behaviour. 
     * This function uses the operating system's awk
	 * @throws InterruptedException 
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

        ArrayList<String> args = prepareAwkScript(awkScript);
        
        ArrayList<String> output = new ArrayList<String>();
        try {
            output = execSystemCommand(rawLines, args);
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
        
        // Check input and output. If an input line is found in output at the
        // same line number where it should be, add True else False.
        int j= 0;
        for(int i=0; i < rawLines.length; i++){
            String inLine= rawLines[i];
            if(output.size() > j){
                String outLine= output.get(j);              
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
    
	/** Stream the raw line through awk and return true if the output of awk is the same as the input
	 * line. If awk returns empty output then the line didn't pass the awk filter.
	 * If output is not empty and not equal to input, return null.
	 * See tests for behaviour. 
     * This function uses built in Java Jawk
	 * */
    /*
    public static boolean[] passAwkFilter(String[] rawLines, String awkScript) throws IOException {
		
		boolean[] results= new boolean[rawLines.length];

		awkScript= awkScript.trim();

		if(awkScript.isEmpty()){
			for(int i= 0; i < rawLines.length; i++){
				results[i]= true;
			}
			return results;
		}

		// * Parse command string into arguments and awk script
		List<String> args= new Tokenizer(awkScript).tokenize();
		
		String script= args.remove(args.size()-1); // 'remove' returns the element removed
		File awkFile= Utils.createTempFile(".asciigenome.", ".awk", true);
		
		BufferedWriter wr= new BufferedWriter(new FileWriter(awkFile));
		wr.write(Track.awkFunc);
		wr.write(script);
		wr.close();
		args.add("-f");
		args.add(awkFile.getAbsolutePath());

		ByteArrayOutputStream baosIn = new ByteArrayOutputStream();
		for (String line : rawLines) {
			baosIn.write((line+"\n").getBytes());
		}
		InputStream is= new ByteArrayInputStream(baosIn.toByteArray());
		PrintStream stdout = System.out;
		ByteArrayOutputStream baos = new ByteArrayOutputStream();

		try{
			PrintStream os= new PrintStream(baos);
			new org.jawk.Main(args.toArray(new String[0]), is, os, System.err);
		} catch(Exception e){
			// e.printStackTrace();
			throw new IOException();
		} finally{
			System.setOut(stdout);
			is.close();
			awkFile.delete();
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
	*/

	/** Right-pad each line in string x with whitespaces. Each line is 
	 * defined by the newline. Each line is padded to become at least of length size.   
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
	public static SamReader getSamReader(String workFilename, String referenceSequence) throws MalformedURLException {
		UrlValidator urlValidator = new UrlValidator();
		SamReaderFactory srf=SamReaderFactory.make();
		if(referenceSequence != null) {
		    srf.referenceSequence(new File(referenceSequence));
		}
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

	/**Get terminal width. 
	 * */
	public static int getTerminalWidth() throws IOException {
		int terminalWidth= jline.TerminalFactory.get().getWidth(); 
		if(terminalWidth <= 0){
			terminalWidth= 80;
		}
		return terminalWidth;
	}

	/**Get VCFHeader from the given source which could be URL or local file.*/
	public static VCFHeader getVCFHeader(String source) throws MalformedURLException{
		VCFHeader vcfHeader;
		if( Utils.urlFileExists(source) ){
			URL url= new URL(source);
			AbstractFeatureReader<VariantContext, LineIterator> reader = AbstractFeatureReader.getFeatureReader(url.toExternalForm(), new VCFCodec(), false);
			vcfHeader = (VCFHeader) reader.getHeader();
		} else {
			VCFFileReader reader = new VCFFileReader(new File(source), false); // Set requiredIndex false!
			vcfHeader= reader.getFileHeader();
			reader.close();
		}
		return vcfHeader;
	}
	
	public static VCFHeaderVersion getVCFHeaderVersion(VCFHeader vcfHeader){
		Iterator<VCFHeaderLine> iter = vcfHeader.getMetaDataInInputOrder().iterator();
		while(iter.hasNext()){
			VCFHeaderLine hl = iter.next();
			if(hl.getKey().equals("fileformat")){
				return VCFHeaderVersion.toHeaderVersion(hl.getValue());
			}
		}
		return null;
	}

	/** Sort and index input sam or bam or cram
	 * @throws IOException 
	 * */
	public static void sortAndIndexSamOrBam(String inSamOrBam, String sortedBam, boolean deleteOnExit, String referenceSequence) throws IOException {

	    if(Utils.isCRAM(sortedBam)) {
	        throw new IOException("CRAM output is not supported");
	    }
	    
		/*  ------------------------------------------------------ */
		/* This chunk prepares SamReader from local bam or URL bam */
		UrlValidator urlValidator = new UrlValidator();
		SamReaderFactory srf=SamReaderFactory.make();
		if(Utils.isCRAM(inSamOrBam)) {
		    srf.referenceSequence(new File(referenceSequence));
		}
		srf.validationStringency(ValidationStringency.SILENT);
		SamReader samReader;
		if(urlValidator.isValid(inSamOrBam)){
			samReader = SamReaderFactory.makeDefault().open(SamInputResource.of(new URL(inSamOrBam)));
		} else {
			samReader= srf.open(new File(inSamOrBam));
		}
		/*  ------------------------------------------------------ */
		
		samReader.getFileHeader().setSortOrder(SortOrder.coordinate);
		
		File out= new File(sortedBam);
		if(deleteOnExit){
			out.deleteOnExit();
			File idx= new File(out.getAbsolutePath().replaceAll("\\.bam$", "") + ".bai");
			idx.deleteOnExit();
		}
		
		SAMFileWriter outputSam= new SAMFileWriterFactory()
				.setCreateIndex(true)
				.makeSAMOrBAMWriter(samReader.getFileHeader(), false, out);

		for (final SAMRecord samRecord : samReader) {
			outputSam.addAlignment(samRecord);
        }
		samReader.close();
		outputSam.close();
	}

	/**True if SAM read names are equal. Read name strings are parsed to remove
	 * parts that are not part of the name.  
	 * */
	public static boolean equalReadNames(String readName, String readName2) {
		return cleanSamReadName(readName).equals(cleanSamReadName(readName2));
	}
	private static String cleanSamReadName(String readName){
		int blank= readName.indexOf(" ");
		if(blank >= 0){
			readName= readName.substring(0, blank);
		}
		if(readName.length() > 2 && (readName.endsWith("/1") || readName.endsWith("/2"))){
			readName= readName.substring(0, readName.length()-2);
		}
		return readName;
	}

	public static List<String> vcfHeaderToStrings(VCFHeader header) {

	    File fakeVCFFile;
	    List<String> str= new ArrayList<String>();
	    try {
	    	fakeVCFFile = Utils.createTempFile(".vcfheader", ".vcf", true);
		    final VariantContextWriter writer = new VariantContextWriterBuilder()
		            .setOutputFile(fakeVCFFile)
		            .setReferenceDictionary(header.getSequenceDictionary())
		            .setOptions(EnumSet.of(Options.ALLOW_MISSING_FIELDS_IN_HEADER, Options.WRITE_FULL_FORMAT_FIELD))
		            .build();
		    writer.writeHeader(header);
		    writer.close();
		    BufferedReader br= new BufferedReader(new FileReader(fakeVCFFile));
		    String line= br.readLine();
		    while(line != null){
		    	str.add(line);
		    	line= br.readLine();
		    }
		    br.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return str;
	}

	/**Parse string x and round the numbers it contains to d decimal places.
	 * Typically, x is a raw line read from BED, GFF whatever. What makes a number
	 * is guessed from context without too much sophistication.
	 * With d < 0, do nothing and return x as it is.  
	 * */
	public static String roundNumbers(final String x, int d, TrackFormat trackFormat) {
		if(d < 0){
			return x;
		}
		// Capture any substring that looks like a number with decimals. Integers are left untouched.
		Pattern p = Pattern.compile("-?\\d+\\.\\d+"); 
		String xm= x;
		if(trackFormat.equals(TrackFormat.GTF)){
			// We consider "123.123" in double quotes as a number to comply with cufflinks/StringTie. 
			xm= xm.replaceAll("\"", " ");
		}
		
		Matcher m = p.matcher(xm);
		List<Integer> starts= new ArrayList<Integer>();
		List<Integer> ends= new ArrayList<Integer>();
		List<String> repls= new ArrayList<String>();
		while (m.find()) {
			if(m.group().replace("-", "").startsWith("00")){
				// We treat something like "000.123" as NaN. 
				continue;
			}
			// If these chars precede the captured string, then the string may be an actual number 
			if(m.start() == 0 || xm.charAt(m.start()-1) == '=' 
					          || xm.charAt(m.start()-1) == ' '
					          || xm.charAt(m.start()-1) == ','
					          || xm.charAt(m.start()-1) == '\t'){
				// If these chars follow the captured string, then the string is an actual number
				if(m.end() == xm.length() || xm.charAt(m.end()) == ';'
							|| xm.charAt(m.end()) == ' '      
							|| xm.charAt(m.end()) == ','
							|| xm.charAt(m.end()) == '\t'){
					DecimalFormat format;
					if(d == 0){
						format = new DecimalFormat("0");
					} else {
						format = new DecimalFormat("0." + StringUtils.repeat("#", d));
					}
					String newval= format.format(Double.valueOf(m.group()));
					starts.add(m.start());
					ends.add(m.end());
					repls.add(newval);
				}
			}
		}
		if(starts.size() == 0){
			return x;
		}
		StringBuilder formattedX= new StringBuilder();
		for(int i= 0; i < starts.size(); i++){
			if(i == 0){
				formattedX.append(x.substring(0, starts.get(i)));	
			} else {
				formattedX.append(x.substring(ends.get(i-1), starts.get(i)));
			}
			formattedX.append(repls.get(i));
		}
		formattedX.append(x.substring(ends.get(ends.size()-1)));
		return formattedX.toString();
	}

	/**Winsorise vector x. Adapted from https://www.r-bloggers.com/winsorization/.
	 * */
	public static List<Float> winsor2(List<Float> x, double multiple) {
				/*
				winsor2<- function (x, multiple=3)
				{
				   med <- median(x)
				   y <- x - med
				   sc <- mad(y, center=0) * multiple
				   y[ y > sc ] <- sc
				   y[ y < -sc ] <- -sc
				   y + med
				}
				*/
		if(multiple <= 0){
			throw new ArithmeticException(); 
		}
		DescriptiveStatistics stats = new DescriptiveStatistics();
		for(float z : x){
			stats.addValue(z);
		}
		float median = (float)stats.getPercentile(50);
		List<Float> y= new ArrayList<Float>(x);
		for(int i= 0; i < x.size(); i++){
			y.set(i, x.get(i) - median);
		}
		float sc= (float) (Utils.medianAbsoluteDeviation(y, 0) * multiple);
		for(int i= 0; i < y.size(); i++){
			if(y.get(i) > sc){
				y.set(i, sc);
			}
			else if(y.get(i) < -sc){
				y.set(i, -sc);
			}
			y.set(i, y.get(i) + median);
		}
		return y;
	}
	/** Translated from R function mad(x, center, constant= 1.4826, na.rm= FALSE, low= FALSE, high= FALSE). 
	 * */
	private static float medianAbsoluteDeviation(List<Float> x, float center) {
		DescriptiveStatistics stats = new DescriptiveStatistics();
		for(int i= 0; i < x.size(); i++){
			stats.addValue(Math.abs(x.get(i) - center));
		}
		float median= (float)stats.getPercentile(50);
		return (float)1.4826 * median;
	}

	/**Interpret x to return boolean. Similar to Boolean.valueOf() but more
	 * flexible. See code and tests for exact behaviour. Case insensitive. 
	 * Throw exception otherwise.
	 * See tests. 
	 * */
	public static Boolean asBoolean(String x) {

		if(x == null || x.trim().isEmpty()){
			throw new IllegalArgumentException();
		}

		x= x.trim().toLowerCase();
		if("true".matches("^" + x + ".*") || "yes".matches("^" + x + ".*") || x.equals("on")){
			return true;
		}
		if("false".matches("^" + x + ".*") || "no".matches("^" + x + ".*") || x.equals("off")){
			return false;
		}
		throw new IllegalArgumentException();
	}
	
	public static List<String> suggestCommand(String hint, List<String> options){
		if(hint.length() < 3){
			return new ArrayList<>();
		}
		hint= hint.toLowerCase();

		Map<Integer, List<String>> candidates= new TreeMap<Integer, List<String>>();
		
		for(String candidate : options){
			if(candidate.length() < 3){
				continue;
			}
			String x= candidate.toLowerCase();
			int nm;
			if(x.length() >= 3 && hint.length() >= 3 && (x.contains(hint) || hint.contains(x))){
				nm= 0;
			} else {
				nm= levenshtein(hint, candidate.toLowerCase());
			}
			if(candidates.containsKey(nm)){
				candidates.get(nm).add(candidate);
			} else {
				List<String> tier= new ArrayList<String>();
				tier.add(candidate);
				candidates.put(nm, tier);
			}
		}
		
		java.util.Map.Entry<Integer, List<String>> out = candidates.entrySet().iterator().next();
		return out.getValue();
	}
	
    /**
     * Calculates Damerau–Levenshtein distance between string {@code a} and
     * {@code b} with given costs.
     * 
     * @param a
     *            String
     * @param b
     *            String
     * @return Damerau–Levenshtein distance between {@code a} and {@code b}
     * 
     * Method from argparse4j
     */
    private static int levenshtein(String a, String b) {
    	
    	 final int SUBSTITUTION_COST = 2;
    	 final int SWAP_COST = 0;
    	 final int DELETION_COST = 4;
    	 final int ADDITION_COST = 1;
    	
        int aLen = a.length();
        int bLen = b.length();
        int[][] dp = new int[3][bLen + 1];
        for (int i = 0; i <= bLen; ++i) {
            dp[1][i] = i;
        }
        for (int i = 1; i <= aLen; ++i) {
            dp[0][0] = i;
            for (int j = 1; j <= bLen; ++j) {
                dp[0][j] = dp[1][j - 1]
                        + (a.charAt(i - 1) == b.charAt(j - 1) ? 0 : SUBSTITUTION_COST);
                if (i >= 2 && j >= 2 && a.charAt(i - 1) != b.charAt(j - 1)
                        && a.charAt(i - 2) == b.charAt(j - 1)
                        && a.charAt(i - 1) == b.charAt(j - 2)) {
                    dp[0][j] = Math.min(dp[0][j], dp[2][j - 2] + SWAP_COST);
                }
                dp[0][j] = Math.min(dp[0][j],
                        Math.min(dp[1][j] + DELETION_COST, dp[0][j - 1] + ADDITION_COST));
            }
            int[] temp = dp[2];
            dp[2] = dp[1];
            dp[1] = dp[0];
            dp[0] = temp;
        }
        return dp[1][bLen];
    }

	public static boolean isOverlapping(int xStart, int xEnd, int yStart, int yEnd) {
		if(xStart > xEnd || yStart > yEnd){
			throw new ArithmeticException();
		}
		if(xEnd < yStart || yEnd < xStart){
			return false;
		}
		return true;
	}

	public static String reformatFileName(String filename, boolean absolute) {
		UrlValidator urlValidator= new UrlValidator();
		if(urlValidator.isValid(filename)){
			return filename;
		}
		Path absfilename= Paths.get(filename).toAbsolutePath();
		if(absolute) {
			return absfilename.toString();
		}
		Path cwd= Paths.get(System.getProperty("user.dir"));
		Path relative= cwd.relativize(absfilename);
		return relative.toString();
	}

	/**rawrecs is an array of raw records and regex is either a regex or 
	 * an awk script (autodetermined). Return an array of booleans for whether 
	 * each record is matched by regex.
	 * @throws IOException 
	 * */
	public static boolean[] matchByAwkOrRegex(String[] rawLines, String regex) throws IOException {
		boolean isAwk= false;
		if(Pattern.compile("\\$\\d|\\$[A-Z]").matcher(regex).find()) {
			// We assume that if regex contains a '$' followed by digit or letter
			// we have an awk script since that would be a very unlikely regex.
			// Could we have an awk script not containing '$'? Unlikely but maybe possible
			isAwk= true;
		}

		boolean[] matched= new boolean[rawLines.length];
		if(isAwk) {
			regex= quote(regex);
//			if(! regex.trim().startsWith("'") && ! regex.trim().endsWith("'")) {
//				regex= "'" + regex + "'";
//			}
			matched= Utils.passAwkFilter(rawLines, regex);
		}
		else {
			for(int i= 0; i < matched.length; i++) {
				matched[i]= Pattern.compile(regex).matcher(rawLines[i]).find();
			}
		}
		return matched;
	}
	
	/**Add quotes around x, if necessary, so that it gets tokenized in a single string*/
	public static String quote(String x) {
		if((x.startsWith("'") && x.endsWith("'")) || (x.startsWith("\"") && x.endsWith("\""))) {
			// No need to quote
			return x;
		}
		// We need to quote the awk script using a quoting string not used inside the script itself.
		// This is a pretty bad hack. You should store the awk command as a list rather than a string
		// split and joined multiple times!
		String q= "";
		if( ! x.contains("'")) {
			q= "'";
		}
		else if(! x.contains("\"")) {
			q= "\"";
		}
		else if(! x.contains("'''")) {
			q= "'''";
		}
		else if(! x.contains("\"\"\"")) {
			q= "\"\"\"";
		}
		else {
			throw new RuntimeException();
		}
		return q + x + q;
	}

	// public static TrackFormat guessTrackType(String file) {
	//    return guessTrackType(file, null);
	// }
	
	public static TrackFormat guessTrackType(String file, String referenceSequence) {
		
		try {
		    SamReader samReader = Utils.getSamReader(file, referenceSequence);
			SAMRecordIterator iter = samReader.iterator();
			int n= 0;
			while(iter.hasNext() && n < 10) {
				iter.next();
				n++;
			}
			return TrackFormat.BAM;
		} catch(Exception e) {
			// e.printStackTrace();
		}
		
		try {
			FeatureList gff = GFF3Reader.read(file);
			Iterator<FeatureI> iter = gff.iterator();
			int n= 0;
			while(iter.hasNext() && n < 1000) {
				iter.next();
				n++;
			}
			return TrackFormat.GFF;
		} catch(Exception e) {
			// e.printStackTrace();
		}
		
		try {
			VCFFileReader reader = new VCFFileReader(new File(file), false);
			CloseableIterator<VariantContext> iter = reader.iterator();
			int n= 0;
			while(iter.hasNext() && n < 1000) {
				try {
					iter.next();
				} catch(IllegalArgumentException e) {
					// Forgive "-" as an allele
				}
				n++;
			}
			reader.close();
			return TrackFormat.VCF;
		} catch(Exception e) {
			e.printStackTrace();
		}
		return null;
	}

    public static boolean isCRAM(String sourceName) {
        return sourceName.toLowerCase().endsWith(".cram");
    } 
	
}
