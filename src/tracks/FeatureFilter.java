package tracks;

import java.util.ArrayList;
import java.util.List;

import htsjdk.samtools.filter.AlignedFilter;
import htsjdk.samtools.filter.SamRecordFilter;

/**Group all the filter settings that can be applied to a feature.
 * */
class FeatureFilter {
	// Grep & Awk
	final static public String HIDE_REGEX= "^$"; protected String hideRegex= HIDE_REGEX;
	final static public String SHOW_REGEX= ".*"; protected String showRegex= SHOW_REGEX;
	private String awk= "";
	// SAM specific
	public static final int f_FLAG= 0; private int f_flag= f_FLAG;
	public static final int F_FLAG= 4; private int F_flag= F_FLAG;
	public static final int MAPQ= 0;   private int mapq= MAPQ;
	private List<SamRecordFilter> samRecordFilter= new ArrayList<SamRecordFilter>(); 
	public static final String VARIANT_CHROM= ""; String variantChrom= VARIANT_CHROM;
	private int variantFrom= -1;
	private int variantTo= -1;
	
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
	
	public void setVariantReadInInterval(String chrom, int from, int to) {
		this.variantChrom= chrom;
		this.variantFrom= from;
		this.variantTo= to;
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
}
