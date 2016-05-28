package tracks;

import htsjdk.variant.variantcontext.VariantContext;

public class IntervalFeatureVCF extends IntervalFeature {
	/* A t t r i b u t e s */
	private String chrom;       // Required
	private int from;           // Required. NB 1 based also for bed files.
	private int to;             // Required 

	private float score= Float.NaN;
	private char strand= '.'; 
	private String source= "."; // Gtf specific
	private String feature= "."; // Gtf specific

	private String raw; // Raw input string exactly as read from file.
	/** Start position of feature in screen coordinates. 
	 * -1 if the feature is not part of the screenshot. */
	private int screenFrom= -1;
	private int screenTo= -1;
	
	private VariantContext variantContext;
	
	public final TrackFormat format= TrackFormat.VCF;

	
	/* C o n s t r u c t o r  */
	
	public IntervalFeatureVCF(VariantContext variantContext){
		this.variantContext= variantContext;
		this.chrom= this.variantContext.getContig();
		this.from= this.variantContext.getStart();
		this.to= this.variantContext.getEnd();
	}
	
	/*  S e t t e r s   and   g e t t e r s */
	
	@Override
	public String getChrom() {
		return variantContext.getContig();
	}
	
	@Override
	public int getFrom() {
		return variantContext.getStart();
	}

	@Override
	public int getTo() {
		return variantContext.getEnd();
	}

	/** This is variantContext.getID() */
	@Override
	public String getName() {
		if(variantContext.hasID()){
			return variantContext.getID();
		} else{
			return null;
		}
	}

	/** This is variantContext.getPhredScaledQual() */
	@Override
	public float getScore() {
		return (float) variantContext.getPhredScaledQual();
	}

	/** Always "." */
	@Override
	public char getStrand() {
		return '.';
	}

	@Override
	public int getScreenFrom() {
		return screenFrom;
	}

	@Override
	public int getScreenTo() {
		return screenTo;
	}
	
	@Override
	public String getSource() {
		return variantContext.getSource();
	}

	/** Not really "raw". This is variantContext.toString() */
	@Override
	public String getRaw() {
		return variantContext.toString();
	}
	
	//public String getFeature() {
	//	return feature;
	//}
}
