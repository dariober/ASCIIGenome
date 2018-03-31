package tracks;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Pattern;
import java.io.File;
import java.io.IOException;

import javax.script.ScriptEngine;
import javax.script.ScriptEngineManager;
import javax.script.ScriptException;

import org.apache.commons.io.FilenameUtils;

import com.google.common.base.Joiner;
import com.google.common.base.Splitter;

import exceptions.InvalidColourException;
import exceptions.InvalidGenomicCoordsException;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

class GenotypeMatrix {

	/** Each `List<Genotype>` is a row in the matrix. The String is the sample name*/
	private Map<String, List<FeatureChar>> matrix= new LinkedHashMap<String, List<FeatureChar>>(); 
    
	/**List of regexes to select samples to be displayed */
	private String selectSampleRegex; 
	
	/** Maximum number of samples (rows) to print */
	private Integer nMaxSamples;
	
	private Map<String, String> subSampleRegex= new HashMap<String, String>();

	private String jsScriptFilter; 

	private ScriptEngine engine; // Leave this null. Only if required use the getter to create it.
	                             // The engine maybe called 1000s of times. So if created only once.
	private final Set<String> genotypes= new HashSet<String>();
	
	
    protected GenotypeMatrix() {
    	this.subSampleRegex.put("pattern", "$^"); // Default pattern-replacement to match nothing
    	this.subSampleRegex.put("replacement", "");
    	this.genotypes.add("{HOM}");
    	this.genotypes.add("{HET}");
    	this.genotypes.add("{HOM_REF}");
    	this.genotypes.add("{HOM_VAR}");
    	this.genotypes.add("{HET_NON_REF}");
    	this.genotypes.add("{CALLED}");
    	this.genotypes.add("{NO_CALL}");
    	this.genotypes.add("{MIXED}");
    }
    
    // -------------------------------------------------------------------------

    /** Replace or create the current matrix using the provided input.
     * @throws InvalidGenomicCoordsException 
     * @throws IOException 
     * */
    private void makeMatrix(List<IntervalFeature> variantList, int terminalWidth, VCFHeader vcfHeader) throws InvalidGenomicCoordsException, IOException{

    	this.matrix= new LinkedHashMap<String, List<FeatureChar>>();
    	
    	if(variantList.size() == 0){
    		return;
    	}

    	Map<VariantContext, String> vcfRecordWithScript= new HashMap<VariantContext, String>();
		if(vcfHeader != null && this.getJsScriptFilter() != null &&  ! this.getJsScriptFilter().trim().isEmpty()){
	    	// We assign to each VCF record the JS script formatted with the fields that
			// do not change across samples, so we do the formatting only once.
	    	for(IntervalFeature ctx : variantList){
	    		String js= this.formatJsScriptWithFixedFields(this.getJsScriptFilter(), ctx.getVariantContext());
	    		js= this.formatJsScriptWithInfo(js, ctx.getVariantContext(), vcfHeader);
	    		vcfRecordWithScript.put(ctx.getVariantContext(), js);
	    	}
		}
		List<String> samples= new ArrayList<String>();
		if(vcfHeader != null){
			samples= vcfHeader.getGenotypeSamples(); 
		} else {
			samples= variantList.get(0).getVariantContext().getSampleNamesOrderedByName();
		}
    	int n= 0;
        for(String sampleName : samples){

    		if(n >= this.getnMaxSamples() && this.getnMaxSamples() >= 0){
    			break;
    		}
        	
    		boolean matched= Pattern.compile(this.getSelectSampleRegex()).matcher(sampleName).find();
    		if( ! matched){
    			continue;
			}
    		
    		boolean keep= true;
    		if(vcfHeader != null && this.getJsScriptFilter() != null &&  ! this.getJsScriptFilter().trim().isEmpty()){
    			keep= this.isPassedFilter(vcfRecordWithScript, sampleName, vcfHeader);
        		if( ! keep ){
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
        		if(col < 0){
        			continue;
        		}
                Genotype gt= variant.getVariantContext().getGenotype(sampleName);
                FeatureChar fmtGt= new FeatureChar();
                fmtGt.addFormatGenotype(gt);
                if(genotypeRow.get(col).getText() == '*'){
                	// This cell has mixed genotype
                	continue;
                }
                if( ! (genotypeRow.get(col).getText() == ' ') && ! (genotypeRow.get(col).getText() == fmtGt.getText())){
                	// If the cell is not empty or the genotype is not the same as this one:
                	fmtGt.setText('*');
                	fmtGt.setBgColor(null);
                	fmtGt.setFgColor(null);
                }
                genotypeRow.set(col, fmtGt);
            }
        	this.matrix.put(sampleName, genotypeRow);
        	n++;
        }
    }
    
    /** Parse sampleNames to remove redundant substring(s)
     * */
    private List<String> cleanSampleNames(List<String> sampleNames) {
    	
    	List<String>cleanNames= new ArrayList<String>();
    	for(String x : sampleNames){
        	// Strip dir name if any and see names are still unique 
    		cleanNames.add(new File(x).getName());
    	}
    	Set<String> unique= new HashSet<String>(cleanNames);

    	if(unique.size() == sampleNames.size()){
    		// Try stripping also extension
        	List<String>cleanNames2= new ArrayList<String>();
        	for(String x : cleanNames){
        		cleanNames2.add(FilenameUtils.removeExtension(x));
        	}
        	Set<String> unique2= new HashSet<String>(cleanNames2);
        	if(unique2.size() == sampleNames.size()){
        		return cleanNames2;
        	} else {
        		return cleanNames;
        	}
    	} else {
        	return sampleNames;
    	}
	}

	/**Return true if ANY of the variants in the given sample pass the filters in javascript.
     * @throws InvalidGenomicCoordsException 
     * @throws IOException */
    private boolean isPassedFilter(Map<VariantContext, String> vcfRecordWithScript, String sampleName, VCFHeader vcfHeader) {

//    	Stopwatch sw= Stopwatch.createUnstarted();
    	boolean subGenotype= false;
    	for(String g : this.genotypes){
    		if(this.getJsScriptFilter().contains(g)){
    			subGenotype= true;
    			break;
    		}
    	}
    	
    	// We concatenate all the scripts for each record in a single string that we execute all in one shot.
    	// In this way we avoid calling the JS engine many times. We accumulate the individual scripts in
    	//  a set where duplicates are excluded so we make the final script simpler and faster to execute.
    	// This is important especially for GT where we have few combinations of genotype appearing many times
    	// across markers.
    	Set<String> concatJS= new HashSet<String>();
    	
    	for(VariantContext ctx : vcfRecordWithScript.keySet()){
    		// We format the JS script and apply it to this sample for each record in the 
    		// window. As soon as a record passes the filter, we pass the sample.
    		String js= vcfRecordWithScript.get(ctx);
    		
    		js= this.formatJsScriptWithFormat(js, sampleName, ctx, vcfHeader);
    		if(subGenotype){
    			js= this.formatJsScriptWithGenotype(js, ctx.getGenotype(sampleName));
    		}
    		concatJS.add("(" + js + ")");  		
    	}
    	String concat= Joiner.on(" || ").join(concatJS);
		Object b = null;
//		sw.start();
		try{
			b = this.getEngine().eval(concat);
    		if((boolean)b){
				return true;
			}
		} catch(ClassCastException | ScriptException e){
			String x= concatJS.size() > 0 ? concatJS.iterator().next() : concat; 
			System.err.println(
					  "ERROR in expression. It may be that the expression is not valid Javascript syntax\n"
					+ "or its result is not a boolean (true or false). First evaluated expression was:\n"
					+  x + ". Final result= " + b + "\n");
			this.setJsScriptFilter(null);
		}
		return false;
	}

	private String formatJsScriptWithGenotype(String js, Genotype gt) {
		
		// Comments are copied from htsjdk 
		if(js.contains("{HOM}")){
			js= js.replace("{HOM}", Boolean.toString(gt.isHom()));
		}
		if(js.contains("{HET}")){
			js= js.replace("{HET}", Boolean.toString(gt.isHet()));
		}
		if(js.contains("{HOM_REF}")){
			js= js.replace("{HOM_REF}", Boolean.toString(gt.isHomRef()));
		}
		if(js.contains("{HOM_VAR}")){
			js= js.replace("{HOM_VAR}", Boolean.toString(gt.isHomVar())); // true if all observed alleles are alt; if any alleles are no-calls, return false.
		}
		if(js.contains("{HET_NON_REF}")){
			js= js.replace("{HET_NON_REF}", Boolean.toString(gt.isHetNonRef())); // true if we're het (observed alleles differ) and neither allele is reference; if the ploidy is less than 2 or if any alleles are no-calls, this method will return false.
		}
		if(js.contains("{CALLED}")){
			js= js.replace("{CALLED}", Boolean.toString(gt.isCalled())); // true if this genotype is comprised of any alleles that are not no-calls (even if some are).
		}
		if(js.contains("{NO_CALL}")){
			js= js.replace("{NO_CALL}", Boolean.toString(gt.isNoCall())); // true if this genotype is not actually a genotype but a "no call" (e.g. './.' in VCF); if any alleles are not no-calls (even if some are), this method will return false.
		}
		if(js.contains("{MIXED}")){
			js= js.replace("{MIXED}", Boolean.toString(gt.isMixed())); // true if this genotype is comprised of both calls and no-calls.
		}
		return js;
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
    	Iterator<VCFInfoHeaderLine> iter = vcfHeader.getInfoHeaderLines().iterator();
    	while(iter.hasNext()){ 
    		//We iterate through each key in the header and see if there is a match in JS script. 
    		VCFInfoHeaderLine headerLine = iter.next();
    		String key= headerLine.getID();
    		if(jsScript.contains('{' + key + '}') || jsScript.contains("{INFO/" + key + '}')){
	    		Object unkValue= ctx.getAttributes().get(key);
	    		String fmtValue;
	    		try{
	    			List<Object> unknList= (List<Object>) unkValue;
	    			StringBuilder listParam= new StringBuilder();
	    			listParam.append("[");
	    			for(Object unk : unknList){
	    				listParam.append(this.formatObjectForJS(key, unk, vcfHeader.getInfoHeaderLine(key).getType()) + ", ");
	    			}
	    			fmtValue= listParam.append("]").toString();
	    		} catch(ClassCastException e){ 
	    			fmtValue= this.formatObjectForJS(key, unkValue, vcfHeader.getInfoHeaderLine(key).getType());
	    		} catch(NullPointerException e){
	    			if(headerLine.getType().equals(VCFHeaderLineType.Flag)){ 
	    				// A flag type returns null if the flag is missing, which is odd. Shouldn't it return false?
	    				fmtValue= "false";  
	    			} else {
	    				fmtValue= "null";
	    			}
	    		}
	    		jsScript= jsScript.replace("{INFO/" + key + '}', fmtValue);
	    		jsScript= jsScript.replace('{' + key + '}', fmtValue);
    		}
    	}
		return jsScript;
    } 

    /**Similar to formatJsScriptWithInfo but applied to FORMAT.*/
    private String formatJsScriptWithFormat(String jsScript, String sampleName, VariantContext ctx, VCFHeader vcfHeader){
    	
    	Iterator<VCFFormatHeaderLine> iter = vcfHeader.getFormatHeaderLines().iterator();
    	while(iter.hasNext()){ 
    		//We iterate through each key in the header and see if there is a match in JS script. 
    		VCFFormatHeaderLine headerLine = iter.next();
    		String key= headerLine.getID();
    		if(jsScript.contains('{' + key + '}') || jsScript.contains("{FMT/" + key + '}')){
    			
    			if(key.equals("GT")){ // We put GT as string rather than list.
    				String gt= this.genotypeAsAlleleIndexes(ctx, sampleName); 
    				jsScript= jsScript.replace("{FMT/" + key + '}', gt);
    				jsScript= jsScript.replace('{' + key + '}', gt);
    				continue;
    			} 
    			
	    		Object unkValue= ctx.getGenotype(sampleName).getAnyAttribute(key);
	    		String fmtValue;
	    		if(headerLine.getCount(ctx) == 1){
	    			fmtValue= this.formatObjectForJS(key, unkValue, headerLine.getType()); 
	    		}
	    		else {
	    			List<String> strList= Splitter.on(",").splitToList(unkValue.toString());
	    			StringBuilder listParam= new StringBuilder();
	    			listParam.append("[");
	    			for(String unk : strList){
	    				if(headerLine.getType().equals(VCFHeaderLineType.String) || 
	    				   headerLine.getType().equals(VCFHeaderLineType.Character)){
	    					listParam.append('"' + unk + '"' + ", ");
	    				} else {
	    					listParam.append(unk + ", "); // Float or Int append as it is w/o quoting
	    				}
	    			}
	    			fmtValue= listParam.append("]").toString();
	    		}
	    		jsScript= jsScript.replace("{FMT/" + key + '}', fmtValue);
	    		jsScript= jsScript.replace('{' + key + '}', fmtValue);
    		}
    	}
		return jsScript;
    } 
    
    /** Get the sample genotype in the same format as it appears on the VCF file.
     * I.e. allele indexes separated by '/' or '|' (if phased). E.g. 0/0, 0/1 etc 
     */
	private String genotypeAsAlleleIndexes(VariantContext ctx, String sample) {
		Genotype gt = ctx.getGenotype(sample);
		char sep= gt.isPhased() ? '|' : '/';
		List<String> all= new ArrayList<String>();
		for(Allele a : gt.getAlleles()){
			if(a.isNoCall()){
				all.add(".");
			}
			else {
				int i= ctx.getAlleleIndex(a);
				all.add(Integer.toString(i));
			}
		}
		return '"' + Joiner.on(sep).join(all) + '"';
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

	protected String printToScreen(boolean noFormat, List<IntervalFeature> linf, int terminalWidth, VCFHeader vcfHeader) throws InvalidColourException, InvalidGenomicCoordsException, IOException {
    	
		this.makeMatrix(linf, terminalWidth, vcfHeader);
		
		StringBuilder sb= new StringBuilder();
    	
		List<String> realNames= new ArrayList<String>(this.matrix.keySet());
		List<String> cleanNames = this.cleanSampleNames(realNames);
		
    	for(int j= 0; j < realNames.size(); j++){
    		String sample= realNames.get(j);
    		String printName= cleanNames.get(j).replaceAll(this.subSampleRegex.get("pattern"), this.subSampleRegex.get("replacement"));
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
		if(this.selectSampleRegex == null){
			return ".*";
		}
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
		if(this.nMaxSamples == null){
			return 10;
		}
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
	
	protected String getJsScriptFilter() {
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
