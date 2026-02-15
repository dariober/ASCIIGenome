package tracks;

import colouring.Config;
import colouring.ConfigKey;
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
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Pattern;
import javax.script.ScriptEngine;
import org.apache.commons.io.FilenameUtils;
import samTextViewer.Utils;

class GenotypeMatrix {

  /** Each `List<Genotype>` is a row in the matrix. The String is the sample name */
  private Map<String, List<FeatureChar>> matrix = new LinkedHashMap<String, List<FeatureChar>>();

  /** List of regexes to select samples to be displayed */
  private String selectSampleRegex;

  /** Maximum number of samples (rows) to print */
  private Integer nMaxSamples;

  private Map<String, String> subSampleRegex = new HashMap<String, String>();

  private String pyScriptFilter;

  private ScriptEngine engine; // Leave this null. Only if required use the getter to create it.
  // The engine maybe called 1000s of times. So if created only once.
  private final Set<String> genotypes = new HashSet<String>();

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

  /**
   * Replace or create the current matrix using the provided input.
   *
   * @throws InvalidGenomicCoordsException
   * @throws IOException
   */
  private void makeMatrix(List<VCFFeature> variantList, int terminalWidth, VCFHeader vcfHeader)
      throws IOException {

    this.matrix = new LinkedHashMap<>();

    if (variantList.isEmpty()) {
      return;
    }

    Map<VariantContext, String> vcfRecordWithScript = new HashMap<VariantContext, String>();
    if (vcfHeader != null
        && this.getPyScriptFilter() != null
        && !this.getPyScriptFilter().trim().isEmpty()) {
      // We assign to each VCF record the py script formatted with the fields that
      // do not change across samples, so we do the formatting only once.
      for (VCFFeature ctx : variantList) {
        String py =
            this.formatPyScriptWithFixedFields(this.getPyScriptFilter(), ctx.getVariantContext());
        py = this.formatPyScriptWithInfo(py, ctx.getVariantContext(), vcfHeader);
        vcfRecordWithScript.put(ctx.getVariantContext(), py);
      }
    }
    List<String> samples = new ArrayList<String>();
    if (vcfHeader != null) {
      samples = vcfHeader.getGenotypeSamples();
    } else {
      samples = variantList.get(0).getVariantContext().getSampleNamesOrderedByName();
    }
    int n = 0;
    for (String sampleName : samples) {

      if (n >= this.getnMaxSamples() && this.getnMaxSamples() >= 0) {
        break;
      }

      boolean matched = Pattern.compile(this.getSelectSampleRegex()).matcher(sampleName).find();
      if (!matched) {
        continue;
      }

      boolean keep = true;
      if (vcfHeader != null
          && this.getPyScriptFilter() != null
          && !this.getPyScriptFilter().trim().isEmpty()) {
        keep = this.isPassedFilter(vcfRecordWithScript, sampleName, vcfHeader);
        if (!keep) {
          continue;
        }
      }

      List<FeatureChar> genotypeRow = new ArrayList<FeatureChar>();
      for (int col = 0; col < terminalWidth; col++) {
        FeatureChar na = new FeatureChar();
        na.setText(' ');
        genotypeRow.add(na); // Initialise row. Potentially there is one genotype per screen column.
      }

      for (VCFFeature variant : variantList) {
        int col = variant.getScreenMid();
        if (col < 0) {
          continue;
        }
        Genotype gt = variant.getVariantContext().getGenotype(sampleName);
        FeatureChar fmtGt = new FeatureChar();
        fmtGt.addFormatGenotype(gt);
        if (genotypeRow.get(col).getText() == '*') {
          // This cell has mixed genotype
          continue;
        }
        if (!(genotypeRow.get(col).getText() == ' ')
            && !(genotypeRow.get(col).getText() == fmtGt.getText())) {
          // If the cell is not empty or the genotype is not the same as this one:
          fmtGt.setText('*');
          fmtGt.setBgColour(null);
          fmtGt.setFgColour(null);
        }
        genotypeRow.set(col, fmtGt);
      }
      this.matrix.put(sampleName, genotypeRow);
      n++;
    }
  }

  /** Parse sampleNames to remove redundant substring(s) */
  private List<String> cleanSampleNames(List<String> sampleNames) {

    List<String> cleanNames = new ArrayList<String>();
    for (String x : sampleNames) {
      // Strip dir name if any and see names are still unique
      cleanNames.add(new File(x).getName());
    }
    Set<String> unique = new HashSet<String>(cleanNames);

    if (unique.size() == sampleNames.size()) {
      // Try stripping also extension
      List<String> cleanNames2 = new ArrayList<String>();
      for (String x : cleanNames) {
        cleanNames2.add(FilenameUtils.removeExtension(x));
      }
      Set<String> unique2 = new HashSet<String>(cleanNames2);
      if (unique2.size() == sampleNames.size()) {
        return cleanNames2;
      } else {
        return cleanNames;
      }
    } else {
      return sampleNames;
    }
  }

  /**
   * Return true if ANY of the variants in the given sample pass the filters in python.
   *
   * @throws InvalidGenomicCoordsException
   * @throws IOException
   */
  private boolean isPassedFilter(
      Map<VariantContext, String> vcfRecordWithScript, String sampleName, VCFHeader vcfHeader) {

    //    	Stopwatch sw= Stopwatch.createUnstarted();
    boolean subGenotype = false;
    for (String g : this.genotypes) {
      if (this.getPyScriptFilter().contains(g)) {
        subGenotype = true;
        break;
      }
    }

    // We concatenate all the scripts for each record in a single string that we execute all in one
    // shot.
    // In this way we avoid calling the JS engine many times. We accumulate the individual scripts
    // in
    //  a set where duplicates are excluded so we make the final script simpler and faster to
    // execute.
    // This is important especially for GT where we have few combinations of genotype appearing many
    // times
    // across markers.
    Set<String> concatPy = new HashSet<String>();

    for (VariantContext ctx : vcfRecordWithScript.keySet()) {
      // We format the JS script and apply it to this sample for each record in the
      // window. As soon as a record passes the filter, we pass the sample.
      String py = vcfRecordWithScript.get(ctx);

      py = this.formatPyScriptWithFormat(py, sampleName, ctx, vcfHeader);
      if (subGenotype) {
        py = this.formatPyScriptWithGenotype(py, ctx.getGenotype(sampleName));
      }
      concatPy.add("(" + py + ")");
    }
    String concat = Joiner.on(" or ").join(concatPy);
    try {
      return this.evalScript(concat);
    } catch (IOException e) {
      this.setPyScriptFilter(null);
      return false;
    }
  }

  private boolean evalScript(String script) throws IOException {
    File tmpScriptFile = Utils.createTempFile(".asciigenome.", ".py", true);
    Files.writeString(tmpScriptFile.toPath(), "print(" + script + ")\n");

    ArrayList<String> cmd = new ArrayList<String>();
    cmd.add(Config.get(ConfigKey.python));
    cmd.add(tmpScriptFile.getAbsolutePath());
    ProcessBuilder pb = new ProcessBuilder().command(cmd);
    pb.redirectErrorStream(true);
    Process p = pb.start();
    BufferedReader reader = new BufferedReader(new InputStreamReader(p.getInputStream()));
    List<String> outLines = new ArrayList<String>();
    String line = "";
    while ((line = reader.readLine()) != null) {
      outLines.add(line.strip());
    }
    reader.close();
    try {
      p.waitFor();
      if (p.exitValue() != 0) {
        throw new IOException(Joiner.on("\n").join(outLines));
      }
    } catch (InterruptedException e) {
      e.printStackTrace();
    }
    if (outLines.size() != 1) {
      throw new IOException("Expected 1 line of output. Got:\n" + Joiner.on("\n").join(outLines));
    }
    String output = outLines.get(0);
    if (output.equals("True")) {
      tmpScriptFile.delete();
      return true;
    }
    if (output.equals("False")) {
      tmpScriptFile.delete();
      return false;
    } else {
      throw new IOException("Expected output to be True or False. Got: '" + output + "'");
    }
  }

  private String formatPyScriptWithGenotype(String py, Genotype gt) {

    // Comments are copied from htsjdk
    if (py.contains("{HOM}")) {
      py = py.replace("{HOM}", gt.isHom() ? "True" : "False");
    }
    if (py.contains("{HET}")) {
      py = py.replace("{HET}", gt.isHet() ? "True" : "False");
    }
    if (py.contains("{HOM_REF}")) {
      py = py.replace("{HOM_REF}", gt.isHomRef() ? "True" : "False");
    }
    if (py.contains("{HOM_VAR}")) {
      py =
          py.replace(
              "{HOM_VAR}",
              gt.isHomVar()
                  ? "True"
                  : "False"); // true if all observed alleles are alt; if any alleles are
      // no-calls, return false.
    }
    if (py.contains("{HET_NON_REF}")) {
      py =
          py.replace(
              "{HET_NON_REF}",
              gt.isHetNonRef()
                  ? "True"
                  : "False"); // true if we're het (observed alleles differ) and neither
      // allele is reference; if the ploidy is less than 2 or if
      // any alleles are no-calls, this method will return false.
    }
    if (py.contains("{CALLED}")) {
      py =
          py.replace(
              "{CALLED}",
              gt.isCalled()
                  ? "True"
                  : "False"); // true if this genotype is comprised of any alleles that are
      // not no-calls (even if some are).
    }
    if (py.contains("{NO_CALL}")) {
      py =
          py.replace(
              "{NO_CALL}",
              gt.isNoCall()
                  ? "True"
                  : "False"); // true if this genotype is not actually a genotype but a "no
      // call" (e.g. './.' in VCF); if any alleles are not no-calls
      // (even if some are), this method will return false.
    }
    if (py.contains("{MIXED}")) {
      py =
          py.replace(
              "{MIXED}",
              gt.isMixed()
                  ? "True"
                  : "False"); // true if this genotype is comprised of both calls and no-calls.
    }
    return py;
  }

  /** Replace the JS script string with the actual values in the fixed VCF fields. */
  private String formatPyScriptWithFixedFields(String py, VariantContext ctx) {

    if (py.contains("{CHROM}")) {
      py = py.replace("{CHROM}", '"' + ctx.getContig() + '"');
    }
    if (py.contains("{POS}")) {
      py = py.replace("{POS}", Integer.toString(ctx.getStart()));
    }
    if (py.contains("{ID}")) {
      py = py.replace("{ID}", '"' + ctx.getID() + '"');
    }
    if (py.contains("{REF}")) {
      py = py.replace("{REF}", '"' + ctx.getReference().getBaseString() + '"');
    }
    if (py.contains("{ALT}")) {
      StringBuilder alleles = new StringBuilder();
      alleles.append("[");
      for (Allele a : ctx.getAlternateAlleles()) {
        alleles.append('"' + a.getBaseString() + '"' + ", ");
      }
      alleles.append("]");
      py = py.replace("{ALT}", alleles.toString());
    }
    if (py.contains("{QUAL}")) {
      py = py.replace("{QUAL}", Double.toString(ctx.getPhredScaledQual()));
    }
    if (py.contains("{FILTER}")) {
      StringBuilder x = new StringBuilder();
      x.append("[");
      if (!ctx.getFilters().isEmpty()) {
        for (String f : ctx.getFilters()) {
          x.append('"' + f + '"' + ", ");
        }
        x.append("]");
        py = py.replace("{FILTER}", x.toString());
      } else {
        py = py.replace("{FILTER}", "[null]");
      }
    }
    return py;
  }

  private List<String> ambiguousInfoFmtTags(VCFHeader vcfHeader) {
    Collection<VCFInfoHeaderLine> infoLines = vcfHeader.getInfoHeaderLines();
    List<String> info = new ArrayList<String>();
    for (VCFInfoHeaderLine line : infoLines) {
      info.add(line.getID());
    }
    Collection<VCFFormatHeaderLine> formatLines = vcfHeader.getFormatHeaderLines();
    List<String> found = new ArrayList<String>();
    for (VCFFormatHeaderLine line : formatLines) {
      String key = line.getID();
      if (info.contains(key)) {
        found.add(key);
      }
    }
    return found;
  }

  /**
   * Replace INFO tags in the jsScript with the actual values found in the variant context object
   */
  @SuppressWarnings("unchecked")
  private String formatPyScriptWithInfo(String pyScript, VariantContext ctx, VCFHeader vcfHeader)
      throws IOException {
    for (String key : this.ambiguousInfoFmtTags(vcfHeader)) {
      if (pyScript.contains("{" + key + "}")) {
        throw new IOException(
            "Key '"
                + key
                + "' found in INFO and FORMAT. Please disambiguate using {FMT/"
                + key
                + "} or {INFO/"
                + key
                + "}.");
      }
    }

    Iterator<VCFInfoHeaderLine> iter = vcfHeader.getInfoHeaderLines().iterator();
    while (iter.hasNext()) {
      // We iterate through each key in the header and see if there is a match in python script.
      VCFInfoHeaderLine headerLine = iter.next();
      String key = headerLine.getID();
      if (pyScript.contains('{' + key + '}') || pyScript.contains("{INFO/" + key + '}')) {
        Object unkValue = ctx.getAttributes().get(key);
        String fmtValue;
        try {
          List<Object> unknList = (List<Object>) unkValue;
          StringBuilder listParam = new StringBuilder();
          listParam.append("[");
          for (Object unk : unknList) {
            listParam.append(
                this.formatObjectForPy(unk, vcfHeader.getInfoHeaderLine(key).getType()) + ", ");
          }
          fmtValue = listParam.append("]").toString();
        } catch (ClassCastException e) {
          fmtValue = this.formatObjectForPy(unkValue, vcfHeader.getInfoHeaderLine(key).getType());
        } catch (NullPointerException e) {
          fmtValue = "None";
          if (headerLine.getType().equals(VCFHeaderLineType.Flag)) {
            // A flag type returns null if the flag is missing, which is odd. Shouldn't it return
            // false?
            if (unkValue == null) {
              fmtValue = "False";
            }
          }
        }
        pyScript = pyScript.replace("{INFO/" + key + '}', fmtValue);
        pyScript = pyScript.replace('{' + key + '}', fmtValue);
      }
    }
    return pyScript;
  }

  /** Similar to formatJsScriptWithInfo but applied to FORMAT. */
  private String formatPyScriptWithFormat(
      String pyScript, String sampleName, VariantContext ctx, VCFHeader vcfHeader) {

    Iterator<VCFFormatHeaderLine> iter = vcfHeader.getFormatHeaderLines().iterator();
    while (iter.hasNext()) {
      // We iterate through each key in the header and see if there is a match in JS script.
      VCFFormatHeaderLine headerLine = iter.next();
      String key = headerLine.getID();
      if (pyScript.contains('{' + key + '}') || pyScript.contains("{FMT/" + key + '}')) {

        if (key.equals("GT")) { // We put GT as string rather than list.
          String gt = this.genotypeAsAlleleIndexes(ctx, sampleName);
          pyScript = pyScript.replace("{FMT/" + key + '}', gt);
          pyScript = pyScript.replace('{' + key + '}', gt);
          continue;
        }

        Object unkValue = ctx.getGenotype(sampleName).getAnyAttribute(key);
        String fmtValue;
        if (headerLine.getCount(ctx) == 1) {
          fmtValue = this.formatObjectForPy(unkValue, headerLine.getType());
        } else {
          List<String> strList = Splitter.on(",").splitToList(unkValue.toString());
          StringBuilder listParam = new StringBuilder();
          listParam.append("[");
          for (String unk : strList) {
            if (headerLine.getType().equals(VCFHeaderLineType.String)
                || headerLine.getType().equals(VCFHeaderLineType.Character)) {
              listParam.append('"' + unk + '"' + ", ");
            } else {
              listParam.append(unk + ", "); // Float or Int append as it is w/o quoting
            }
          }
          fmtValue = listParam.append("]").toString();
        }
        pyScript = pyScript.replace("{FMT/" + key + '}', fmtValue);
        pyScript = pyScript.replace('{' + key + '}', fmtValue);
      }
    }
    return pyScript;
  }

  /**
   * Get the sample genotype in the same format as it appears on the VCF file. I.e. allele indexes
   * separated by '/' or '|' (if phased). E.g. 0/0, 0/1 etc
   */
  private String genotypeAsAlleleIndexes(VariantContext ctx, String sample) {
    Genotype gt = ctx.getGenotype(sample);
    char sep = gt.isPhased() ? '|' : '/';
    List<String> all = new ArrayList<String>();
    for (Allele a : gt.getAlleles()) {
      if (a.isNoCall()) {
        all.add(".");
      } else {
        int i = ctx.getAlleleIndex(a);
        all.add(Integer.toString(i));
      }
    }
    return '"' + Joiner.on(sep).join(all) + '"';
  }

  /**
   * Return Object unk as a string quoted or not quoted depending on its type and suitable for a
   * python script.
   */
  private String formatObjectForPy(Object unk, VCFHeaderLineType type) {
    if (unk == null) {
      return "None"; // Can this actually happen?
    }
    if (type.equals(VCFHeaderLineType.Flag)) {
      return (boolean) unk ? "True" : "False";
    }
    if (type.equals(VCFHeaderLineType.Integer) || type.equals(VCFHeaderLineType.Float)) {
      return unk.toString();
    } else {
      return '"' + unk.toString() + '"';
    }
  }

  protected String printToScreen(
      boolean noFormat, List<VCFFeature> linf, int terminalWidth, VCFHeader vcfHeader)
      throws InvalidColourException, InvalidGenomicCoordsException, IOException {

    this.makeMatrix(linf, terminalWidth, vcfHeader);

    StringBuilder sb = new StringBuilder();

    List<String> realNames = new ArrayList<String>(this.matrix.keySet());
    List<String> cleanNames = this.cleanSampleNames(realNames);

    for (int j = 0; j < realNames.size(); j++) {
      String sample = realNames.get(j);
      String printName =
          cleanNames
              .get(j)
              .replaceAll(
                  this.subSampleRegex.get("pattern"), this.subSampleRegex.get("replacement"));
      List<FeatureChar> fmtName = this.formatName(printName);
      List<FeatureChar> row = this.matrix.get(sample);
      for (int i = 0; i < row.size(); i++) {
        if (i < fmtName.size()) {
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

  private List<FeatureChar> formatName(String name) {
    List<FeatureChar> fmt = new ArrayList<FeatureChar>();
    for (int i = 0; i < name.length(); i++) {
      FeatureChar c = new FeatureChar();
      c.setText(name.charAt(i));
      fmt.add(c);
    }
    return fmt;
  }

  private String getSelectSampleRegex() {
    if (this.selectSampleRegex == null) {
      return ".*";
    }
    return selectSampleRegex;
  }

  /** Set list of regexes to capture sample names. If null or 0-length, reset to default ".*" */
  protected void setSelectSampleRegex(String selectSampleRegex) {
    if (selectSampleRegex == null) {
      this.selectSampleRegex = ".*";
    } else {
      this.selectSampleRegex = selectSampleRegex;
    }
  }

  private int getnMaxSamples() {
    if (this.nMaxSamples == null) {
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

  protected void setPyScriptFilter(String pyScriptFilter) {
    this.pyScriptFilter = pyScriptFilter;
  }

  protected String getPyScriptFilter() {
    return this.pyScriptFilter;
  }
}
