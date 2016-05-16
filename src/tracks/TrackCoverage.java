package tracks;

import java.io.File;
import java.io.IOException;
import java.net.URL;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import org.apache.commons.lang3.StringUtils;
import org.apache.commons.validator.routines.UrlValidator;

import com.google.common.base.Joiner;

import htsjdk.samtools.BAMIndex;
import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import samTextViewer.GenomicCoords;
import samTextViewer.SamLocusIterator;
import samTextViewer.Utils;

@SuppressWarnings("deprecation")
public class TrackCoverage extends Track {

	/* A t t r i b u t e s */
	
	private List<ScreenLocusInfo> screenLocusInfoList= new ArrayList<ScreenLocusInfo>(); 

	/** Each dot in the screen output track corresponds to this many units of 
	 * score in the input. Typically this "reads per dot". */
	// private double scorePerDot;   
	private boolean rpm= false;
	
	
	/* C o n s t r u c t o r */
	
	/**
	 * Construct coverage track from bam alignment in the provided interval. 
	 * Loci will be sampled according to the size of the interval and the size of the printable screen. 
	 * @param bam Input bam file
	 * @param gc Interval to sample positions from
	 * @param windowSize The size of the screen in number of characters.
	 * @param filters Record filters to apply to input sam records.
	 * @param bs Should loci be parsed also as BS-Seq data? 
	 * @throws IOException 
	 */
	public TrackCoverage(String bam, GenomicCoords gc,
			List<SamRecordFilter> filters, boolean bs) throws IOException{
		
		this.setGc(gc);
		this.setFilename(bam);
		this.setFilters(filters);
		this.setBs(bs);
		this.update();
	}
	
	/* M e t h o d s */

    //private SeekableStream myIndexSeekableStream() {
    //   throw new UnsupportedOperationException();
    // }
	
	public void update() throws IOException{
		
		this.screenLocusInfoList= new ArrayList<ScreenLocusInfo>();
		if(this.getGc().getGenomicWindowSize() < this.MAX_REGION_SIZE){
			
			/*  ------------------------------------------------------ */
			/* This chunk prepares SamReader from local bam or URL bam */
			UrlValidator urlValidator = new UrlValidator();
			SamReaderFactory srf=SamReaderFactory.make();
			srf.validationStringency(ValidationStringency.SILENT);
			SamReader samReader;
			if(urlValidator.isValid(this.getFilename())){
				samReader = srf.open(SamInputResource.of(new URL(this.getFilename())).index(new URL(this.getFilename() + ".bai")));
			} else {
				samReader= srf.open(new File(this.getFilename()));
			}
			/*  ------------------------------------------------------ */
			
			IntervalList il= new IntervalList(samReader.getFileHeader());
			il.add(new Interval(this.getGc().getChrom(), this.getGc().getFrom(), this.getGc().getTo()));
			SamLocusIterator samLocIter= new SamLocusIterator(samReader, il, true);
			samLocIter.setSamFilters(this.getFilters());
			Iterator<samTextViewer.SamLocusIterator.LocusInfo> iter= samLocIter.iterator();
		
			for(int i= 0; i < this.getGc().getMapping().size(); i++){
				this.screenLocusInfoList.add(new ScreenLocusInfo());	
			}
		
			while(iter.hasNext()){			
				samTextViewer.SamLocusIterator.LocusInfo locusInfo= iter.next();
				int screenPos= Utils.getIndexOfclosestValue(locusInfo.getPosition(), this.getGc().getMapping());
				byte refBase= '\0';
				if(this.getGc().getRefSeq() != null){
					refBase= this.getGc().getRefSeq()[screenPos];
				}
				this.screenLocusInfoList.get(screenPos).increment(locusInfo, refBase, this.isBs());
			}
			samLocIter.close();
			samReader.close();	
		} 
		this.setYLimitMin(this.getYLimitMin());
		this.setYLimitMax(this.getYLimitMax());
	}
	
	/**
	 * Printable coverage track. The height of the track in lines is `yMaxLines`.
	 * @param screenToGenomeMap List of genomic positions corresponding to each column on screen.
	 * @param yMaxLines
	 * @param rpm Should read counts be normalized by library size as Reads Per Million
	 * @return HashMapwith with keys/values the printable characteristics of the track. 
	 */
	@Override
	public String printToScreen(){
				
		if(this.getyMaxLines() == 0){
			return "";
		} else if(this.screenLocusInfoList.size() == 0){
			if(this.getGc().getGenomicWindowSize() >= this.MAX_REGION_SIZE){
				return "Track not shown: Window is too large";
			}
			return "";
		}
		
		List<Double> yValues= new ArrayList<Double>();
		for(ScreenLocusInfo x : screenLocusInfoList){
			yValues.add(x.getMeanDepth());
		}
		this.setScreenScores(yValues);
				
		// this.scorePerDot= textProfile.getScorePerDot();
		if(this.rpm){
			long libSize= getAlignedReadCount(new File(this.getFilename()));
			// this.scorePerDot= this.scorePerDot / libSize * 1000000.0;
			for(int i= 0; i < yValues.size(); i++){
				yValues.set(i, yValues.get(i)/libSize * 1000000.0);
			}
		}

		TextProfile textProfile= new TextProfile(yValues, this.getyMaxLines(), this.getYLimitMin(), this.getYLimitMax());
		ArrayList<String> lineStrings= new ArrayList<String>();
		for(int i= (textProfile.getProfile().size() - 1); i >= 0; i--){
			List<String> xl= textProfile.getProfile().get(i);
			lineStrings.add(StringUtils.join(xl, ""));
		}
		return Joiner.on("\n").join(lineStrings);
	}
        
    private long getAlignedReadCount(File bam){

    	SAMFileReader.setDefaultValidationStringency(ValidationStringency.SILENT);
		@SuppressWarnings("resource")
		BAMIndex sr= new SAMFileReader(bam).getIndex();
		long alnCount= 0; 
		int i= 0;
		while(true){
			try{
				alnCount += sr.getMetaData(i).getAlignedRecordCount();
			} catch(NullPointerException e){
				break;
			}
			i++;
		}
		sr.close();
		return alnCount;
    }
    
    /* S e t t e r s   and   G e t t e r s */
    
    /* This method makes sense calling only after having set the profile. Typically after */
    public List<ScreenLocusInfo> getScreenLocusInfoList(){
    	return screenLocusInfoList;
    }

	public void setRpm(boolean rpm) {
		this.rpm = rpm;
	}

	@Override
	public String getTitle(){
		
		double[] rounded= Utils.roundToSignificantDigits(this.getMinScreenScores(), this.getMaxScreenScores(), 2);
		
		// String s= Double.toString(Utils.roundToSignificantFigures(this.scorePerDot, 4));
		// String scoreXDot= s.indexOf(".") < 0 ? s : s.replaceAll("0*$", "").replaceAll("\\.$", "");

		return this.getFileTag() 
				+ "; ylim[" + this.getYLimitMin() + " " + this.getYLimitMax() + "]" 
				+ "; range[" + rounded[0] + " " + rounded[1] + "]\n";
				// + "; .= " + scoreXDot + ";\n";
	}
}
