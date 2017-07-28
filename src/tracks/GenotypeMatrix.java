package tracks;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Pattern;

import javax.script.ScriptEngine;
import javax.script.ScriptEngineManager;
import javax.script.ScriptException;

import com.google.common.base.Joiner;

import exceptions.InvalidColourException;
import exceptions.InvalidGenomicCoordsException;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;

class GenotypeMatrix {

	/** Each `List<Genotype>` is a row in the matrix. The String is the sample name*/
	private Map<String, List<FeatureChar>> matrix= new LinkedHashMap<String, List<FeatureChar>>(); 
    
	/**List of regexes to select samples to be displayed */
	private String selectSampleRegex= ".*"; 
	
	/** Maximum number of samples (rows) to print */
	private int nMaxSamples= 10;

	
	private Map<String, String> subSampleRegex= new HashMap<String, String>();

	private String jsScriptFilter; 

	private ScriptEngine engine; // Leave this null. Only if required use the getter to create it.
	                             // The engine maybe called 1000s of times. So if created only once.
	
    protected GenotypeMatrix() {
    	this.subSampleRegex.put("pattern", "$^"); // Default pattern-replacement to match nothing
    	this.subSampleRegex.put("replacement", "");
    }
    
    // -------------------------------------------------------------------------

    /** Replace or create the current matrix using the provided input.
     * @throws InvalidGenomicCoordsException 
     * */
    private void makeMatrix(List<IntervalFeature> variantList, int terminalWidth, VCFHeader vcfHeader) throws InvalidGenomicCoordsException{

    	this.matrix= new LinkedHashMap<String, List<FeatureChar>>();
    	
    	if(variantList.size() == 0){
    		return;
    	}

    	int n= 0;
        for(String sampleName : variantList.get(0).getVariantContext().getSampleNamesOrderedByName()){

    		if(n >= this.getnMaxSamples() && this.getnMaxSamples() >= 0){
    			break;
    		}
        	
    		boolean matched= Pattern.compile(this.getSelectSampleRegex()).matcher(sampleName).find();
    		if( ! matched){
    			continue;
			}
    		
    		boolean keep= true;
    		if(vcfHeader != null && this.getJsScriptFilter() != null){
    			keep= this.isPassedFilter(variantList, sampleName, vcfHeader);
        		if( ! keep){
        			continue;
        		}
    		}
    		
        	List<FeatureChar> genotypeRow= new ArrayList<FeatureChar>();
        	for(int col= 0; col < terminalWidth; col++){
            	FeatureChar na= new FeatureChar();
        		na.setText(' ');
        		genotypeRow.add(na); // Initialise row. Potentially there is one genotype per screen column.
        	}
        	
        	for(IntervalFeature variant : variantList){
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
    
    /**Return true if ANY of the variants in the given sample pass the filters in javascript.
     * @throws InvalidGenomicCoordsException */
    private boolean isPassedFilter(List<IntervalFeature> variantList, String sampleName, VCFHeader vcfHeader) throws InvalidGenomicCoordsException {

		StringBuilder concatJS= new StringBuilder(); 
    	for(IntervalFeature vcf : variantList){
    		String js= this.getJsScriptFilter();
    		VariantContext ctx= vcf.getVariantContext();
    		
    		js= this.formatJsScriptWithFixedFields(js, ctx);
    		js= this.formatJsScriptWithFormat(js, sampleName, ctx, vcfHeader);
    		js= this.formatJsScriptWithInfo(js, ctx, vcfHeader);
    		concatJS.append(js);
    		concatJS.append(" || ");
    	}

    	String script= concatJS.toString().replaceAll(" \\|\\| $", ""); 
		
		Object b = null;
		try {
			b = this.getEngine().eval(script);
		} catch (ScriptException e) {
			e.printStackTrace();
		}
		try{
    		if((boolean)b){
				return true;
			}
		} catch(ClassCastException e){
			throw new InvalidGenomicCoordsException();
		}

    	return false;
	}

	/**Replace the JS script string with the actual values in the fixed VCF fields. */
    private String formatJsScriptWithFixedFields(String js, VariantContext ctx) {
    	if(js.contains("{CHROM}")){
    		js= js.replace("{CHROM}", '"' + ctx.getContig() + '"');
    	}
    	if(js.contains("{POS}")){
    		js= js.replace("{POS}", Integer.toString(ctx.getStart()));
    	}
    	if(js.contains("{ID}")){
    		js= js.replace("{ID}",  '"' + ctx.getID() + '"');
    	}
    	if(js.contains("{REF}")){
    		js= js.replace("{REF}", '"' + ctx.getReference().getBaseString() + '"');
    	}
    	if(js.contains("{ALT}")){
    		StringBuilder alleles= new StringBuilder();
    		alleles.append("[");
    		for(Allele a : ctx.getAlternateAlleles()){
    			alleles.append('"' + a.getBaseString() + '"' + ", ");
    		}
    		alleles.append("]");
    		js= js.replace("{ALT}", alleles.toString());
    	}
    	if(js.contains("{QUAL}")){
    		js= js.replace("{QUAL}", Double.toString(ctx.getPhredScaledQual()));
    	}
    	if(js.contains("{FILTER}")){
    		StringBuilder x= new StringBuilder();
    		x.append("[");
    		if(ctx.getFilters().size() > 0){
        		for(String f : ctx.getFilters()){
        			x.append('"' + f + '"' + ", ");
        		}
        		x.append("]");
        		js= js.replace("{FILTER}", x.toString());
    		}
    		else {
    			js= js.replace("{FILTER}", "[null]");
    		}
    	}
		return js;
	}

	/** Replace INFO tags in the jsScript with the actual values found in the variant context object
     * */
    @SuppressWarnings("unchecked")
	private String formatJsScriptWithInfo(String jsScript, VariantContext ctx, VCFHeader vcfHeader){
    	for(String key : ctx.getAttributes().keySet()){
    		if(jsScript.contains('{' + key + '}') || jsScript.contains("{INFO/" + key + '}')){
	    		Object unknValue= ctx.getAttributes().get(key);
	    		String fmtValue;
	    		try{
	    			List<Object> unknList= (List<Object>) unknValue;
	    			StringBuilder listParam= new StringBuilder();
	    			listParam.append("[");
	    			for(Object unk : unknList){
	    				listParam.append(this.formatObjectForJS(key, unk, vcfHeader.getInfoHeaderLine(key).getType()) + ", ");
	    			}
	    			fmtValue= listParam.append("]").toString();
	    		} catch(ClassCastException e){ 
	    			fmtValue= this.formatObjectForJS(key, unknValue, vcfHeader.getInfoHeaderLine(key).getType());
	    		}
	    		jsScript= jsScript.replace("{INFO/" + key + '}', fmtValue);
	    		jsScript= jsScript.replace('{' + key + '}', fmtValue);
    		}
    	}
		return jsScript;
    } 

    /**Similar to formatJsScriptWithInfo but applied to FORMAT.*/
    @SuppressWarnings("unchecked")
	private String formatJsScriptWithFormat(String jsScript, String sampleName, VariantContext ctx, VCFHeader vcfHeader){
    	Iterator<VCFFormatHeaderLine> iter = vcfHeader.getFormatHeaderLines().iterator();
    	while(iter.hasNext()){ 
    		//We iterate through each key in the header and see if there is a match in JS script. 
    		VCFFormatHeaderLine headerLine = iter.next();
    		String key= headerLine.getID();
    		if(jsScript.contains('{' + key + '}') || jsScript.contains("{FMT/" + key + '}')){
	    		Object unkValue= ctx.getGenotype(sampleName).getAnyAttribute(key);
	    		String fmtValue;
	    		try{
	    			List<Object> unknList= (List<Object>) unkValue;
	    			StringBuilder listParam= new StringBuilder();
	    			listParam.append("[");
	    			for(Object unk : unknList){
	    				listParam.append(this.formatObjectForJS(key, unk, headerLine.getType()) + ", ");
	    			}
	    			fmtValue= listParam.append("]").toString();
	    		} catch(ClassCastException e){ 
	    			fmtValue= this.formatObjectForJS(key, unkValue, headerLine.getType());
	    		}
	    		jsScript= jsScript.replace("{FMT/" + key + '}', fmtValue);
	    		jsScript= jsScript.replace('{' + key + '}', fmtValue);
    		}
    	}
		return jsScript;
    } 
    
    /**Return Object unk as a string quoted or not quoted depending on its type 
     * and suitable for a javascript script.
     * */
	private String formatObjectForJS(String key, Object unk, VCFHeaderLineType type) {
		
		if(unk == null){
			return "null"; // Can this actually happen?
		}
		if(type.equals(VCFHeaderLineType.Flag) || 
		   type.equals(VCFHeaderLineType.Integer) || 
		   type.equals(VCFHeaderLineType.Float)){
			return unk.toString();
		} else {
			return '"' + unk.toString() + '"';
		}	
	}

	protected String printToScreen(boolean noFormat, List<IntervalFeature> linf, int terminalWidth, VCFHeader vcfHeader) throws InvalidColourException, InvalidGenomicCoordsException {
    	
		this.makeMatrix(linf, terminalWidth, vcfHeader);
		
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

	protected void setJsScriptFilter(String jsScriptFilter) {
		this.jsScriptFilter= jsScriptFilter;
	}
	
	private String getJsScriptFilter() {
		return this.jsScriptFilter;
	}
	
    private ScriptEngine getEngine() {
		if(this.engine == null){
	    	ScriptEngineManager factory = new ScriptEngineManager();
	    	this.engine = factory.getEngineByName("JavaScript");
		}
		return this.engine;
	}
}
