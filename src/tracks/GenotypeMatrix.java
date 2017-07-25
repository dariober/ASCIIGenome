package tracks;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Pattern;

import exceptions.InvalidColourException;
import htsjdk.variant.variantcontext.Genotype;

class GenotypeMatrix {

	/** Each `List<Genotype>` is a row in the matrix. The String is the sample name*/
	private Map<String, List<FeatureChar>> matrix= new LinkedHashMap<String, List<FeatureChar>>(); 
    
	/**List of regexes to select samples to be displayed */
	private String selectSampleRegex= ".*"; 
	
	/** Maximum number of samples (rows) to print */
	private int nMaxSamples= 10;

	
	private Map<String, String> subSampleRegex= new HashMap<String, String>();
	
    protected GenotypeMatrix() {
    	this.subSampleRegex.put("pattern", "$^"); // Default pattern-replacement to match nothing
    	this.subSampleRegex.put("replacement", "");
    }
    
    // -------------------------------------------------------------------------

    /** Replace or create the current matrix using the provided input.
     * */
    protected void makeMatrix(List<IntervalFeature> linf, int terminalWidth){

    	this.matrix= new LinkedHashMap<String, List<FeatureChar>>();
    	
    	if(linf.size() == 0){
    		return;
    	}

    	int n= 0;
        for(String sampleName : linf.get(0).getVariantContext().getSampleNamesOrderedByName()){
        
    		boolean matched= Pattern.compile(this.getSelectSampleRegex()).matcher(sampleName).find();
    		if( ! matched){
    			continue;
			}
        	
    		if(n >= this.getnMaxSamples() && this.getnMaxSamples() >= 0){
    			break;
    		}
        	
        	List<FeatureChar> genotypeRow= new ArrayList<FeatureChar>();
        	for(int col= 0; col < terminalWidth; col++){
            	FeatureChar na= new FeatureChar();
        		na.setText(' ');
        		genotypeRow.add(na); // Initialise row. Potentially there is one genotype per screen column.
        	}
        	
        	for(IntervalFeature variant : linf){
        		int col= variant.getScreenMid();
                Genotype gt= variant.getVariantContext().getGenotype(sampleName);
                FeatureChar fmtGt= new FeatureChar();
                fmtGt.addFormatGenotype(gt);
                genotypeRow.set(col, fmtGt);
            }
        	this.matrix.put(sampleName, genotypeRow);
        	n++;
        }
    }
    
    protected String printToScreen(boolean noFormat) throws InvalidColourException {
    	StringBuilder sb= new StringBuilder();
    	
    	for(final String sample : this.matrix.keySet()){
    		String printName= sample.replaceAll(this.subSampleRegex.get("pattern"), this.subSampleRegex.get("replacement"));
    		List<FeatureChar> fmtName = this.formatName(printName); 
    		List<FeatureChar> row = this.matrix.get(sample);
    		for(int i= 0; i < row.size(); i++){
    			if(i < fmtName.size()){
    				sb.append(fmtName.get(i).format(noFormat));
    			} else {
    				FeatureChar gt = this.matrix.get(sample).get(i);
        			sb.append(gt.format(noFormat));	
    			}
    		}
    		sb.append("\n");
    		// Limit by number of samples.
    	}
    	return sb.toString().replaceAll("\n$", "");
    }
    
    private List<FeatureChar> formatName(String name){
    	List<FeatureChar> fmt= new ArrayList<FeatureChar>();
    	for(int i= 0; i < name.length(); i++){
    		FeatureChar c= new FeatureChar();
    		c.setText(name.charAt(i));
    		fmt.add(c);
    	}
    	return fmt;
    }
    
	private String getSelectSampleRegex() {
		return selectSampleRegex;
	}

	/**Set list of regexes to capture sample names. If null or 0-length,
	 * reset to default ".*"
	 * */
	protected void setSelectSampleRegex(String selectSampleRegex) {
		if(selectSampleRegex == null){
			this.selectSampleRegex= ".*";
		} else {
			this.selectSampleRegex = selectSampleRegex;
		}
	}

	private int getnMaxSamples() {
		return nMaxSamples;
	}

	protected void setnMaxSamples(int nMaxSamples) {
		this.nMaxSamples = nMaxSamples;
	}

	protected void setSubSampleRegex(String pattern, String replacement) {
		this.subSampleRegex.put("pattern", pattern);
		this.subSampleRegex.put("replacement", replacement);
	}

}
