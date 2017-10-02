package tracks;

import java.util.ArrayList;
import java.util.List;

import htsjdk.samtools.filter.AlignedFilter;
import htsjdk.samtools.filter.SamRecordFilter;

/**Group all the filter settings that can be applied to a feature.
 * */
class FeatureFilter {
	// Grep & Awk
	final static public String DEFAULT_HIDE_REGEX= "^$"; protected String hideRegex= DEFAULT_HIDE_REGEX;
	final static public String DEFAULT_SHOW_REGEX= ".*"; protected String showRegex= DEFAULT_SHOW_REGEX;
	final static public String DEFAULT_AWK= ""; private String awk= DEFAULT_AWK;
	// SAM specific
	public static final int DEFAULT_f_FLAG= 0; private int f_flag= DEFAULT_f_FLAG;
	public static final int DEFAULT_F_FLAG= 4; private int F_flag= DEFAULT_F_FLAG;
	public static final int DEFAULT_MAPQ= 0;   private int mapq= DEFAULT_MAPQ;
	private List<SamRecordFilter> samRecordFilter= new ArrayList<SamRecordFilter>(); 
	public static final String DEFAULT_VARIANT_CHROM= ""; String variantChrom= DEFAULT_VARIANT_CHROM;
	private int variantFrom= -1;
	private int variantTo= -1;
	private byte[] faSeq;
	
	//   S E T T E R S    A N D    G E T T E R S
	
	public void setShowHideRegex(String showRegex, String hideRegex) {
		this.showRegex= showRegex;
		this.hideRegex= hideRegex;
	}

	public String getHideRegex() { 
		return this.hideRegex; 
	}
	
	public String getShowRegex() { 
		return this.showRegex; 
	}

	public String getAwk() {
		return awk;
	}

	public void setAwk(String awk) {
		this.awk = awk;
	}

	/** This int is just a setting but is NOT translated to a filter! */
	protected int get_f_flag() {
		return f_flag;
	}

	/** This int is just a setting but is NOT translated to a filter! */
	protected void set_f_flag(int f_flag) {
		this.f_flag = f_flag;
	}

	/** This int is just a setting but is NOT translated to a filter! */
	protected int get_F_flag() {
		return F_flag;
	}

	/** This int is just a setting but is NOT translated to a filter! */
	protected void set_F_flag(int F_flag) {
		this.F_flag = F_flag;
	}

	public int getMapq() {
		return mapq;
	}

	public void setMapq(int mapq) {
		this.mapq = mapq;
	}

	/** Return filter making sure the AlignedFilter to discard unmapped is set.
	 * */
	public List<SamRecordFilter> getSamRecordFilter() { 
		AlignedFilter unmapped = new AlignedFilter(true);
		if(!this.samRecordFilter.contains(unmapped)){
			this.samRecordFilter.add(unmapped); // Unmapped reads are always discarded	
		}
		return this.samRecordFilter; 
	}

	public void setSamRecordFilter(List<SamRecordFilter> samRecordFilter) {
		this.samRecordFilter = samRecordFilter;
	}
	
	public void setVariantReadInInterval(String chrom, int from, int to, byte[] faSeq) {
		this.variantChrom= chrom;
		this.variantFrom= from;
		this.variantTo= to;
		this.faSeq= faSeq;
	}
	
	public String getVariantChrom(){
		return this.variantChrom;
	}	
	public int getVariantFrom(){
		return this.variantFrom;
	}
	public int getVariantTo(){
		return this.variantTo;
	}

	public byte[] getFaSeq() {
		return faSeq;
	}

//	public void setFaSeq(byte[] faSeq) {
//		this.faSeq = faSeq;
//	}	
}
