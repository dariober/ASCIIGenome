package tracks;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFCodec;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;
import samTextViewer.Utils;

/**
 * Class to hold bed or gtf features. Behaviour should be similar to pybedtools Interval. Feature
 * coords are 1-based. The first ten bases of the chrom have from-to = 1-10
 *
 * @author berald01
 */
public class VCFFeature extends IntervalFeature {
  private VariantContext variantContext;

  /* C o n s t r u c t o r s */

  /**
   * Create an IntervalFeature from a String. Typically this is a line read from file. vcfHeader can
   * be null if trackformat is not VCF.
   */
  public VCFFeature(String line, VCFCodec vcfCodec) {
    this.variantContext = vcfCodec.decode(line);
    this.setChrom(this.variantContext.getContig());
    this.setFrom(this.setFromForVCF());
    this.setTo(this.setToForVCF());
    this.setName(this.variantContext.getID());
    this.setRaw(line);
  }

  /* M e t h o d s */

  /** Replace the existing variant context with a copy having attributes and samples removed */
  protected void stripVariantContext(
      Set<String> infoToRemove, Set<String> formatToRemove, Set<String> samplesToKeep) {

    if (this.variantContext == null) {
      return;
    }

    VariantContext originalVC = this.variantContext;

    VariantContextBuilder vcb = new VariantContextBuilder(originalVC);

    /*
     * -------------------------
     * 1) Strip INFO fields
     * -------------------------
     */
    if (infoToRemove != null && !infoToRemove.isEmpty()) {
      Map<String, Object> newAttributes =
          originalVC.getAttributes().entrySet().stream()
              .filter(e -> !infoToRemove.contains(e.getKey()))
              .collect(Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue));

      vcb.attributes(newAttributes);
    }

    /*
     * -------------------------
     * 2) Filter samples
     * -------------------------
     */
    Collection<Genotype> workingGenotypes;

    if (samplesToKeep == null) {
      // keep all
      workingGenotypes = originalVC.getGenotypes();
    } else if (samplesToKeep.isEmpty()) {
      // keep none
      workingGenotypes = Collections.emptyList();
    } else {
      // keep only specified samples
      workingGenotypes = originalVC.getGenotypes(samplesToKeep);
    }

    /*
     * -------------------------
     * 3) Strip FORMAT fields
     * -------------------------
     */
    if (formatToRemove != null && !formatToRemove.isEmpty()) {

      List<Genotype> newGenotypes = new ArrayList<>(workingGenotypes.size());

      for (Genotype g : workingGenotypes) {

        GenotypeBuilder gb = new GenotypeBuilder(g);

        // Handle common FORMAT fields explicitly
        if (formatToRemove.contains("DP")) gb.noDP();
        if (formatToRemove.contains("AD")) gb.noAD();
        if (formatToRemove.contains("PL")) gb.noPL();
        if (formatToRemove.contains("GQ")) gb.noGQ();

        // Remove custom/extended FORMAT fields
        Map<String, Object> filteredExtended =
            g.getExtendedAttributes().entrySet().stream()
                .filter(e -> !formatToRemove.contains(e.getKey()))
                .collect(Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue));

        gb.attributes(filteredExtended);

        newGenotypes.add(gb.make());
      }

      vcb.genotypes(newGenotypes);

    } else {
      vcb.genotypes(workingGenotypes);
    }

    /*
     * -------------------------
     * 4) Build new immutable VariantContext
     * -------------------------
     */
    this.variantContext = vcb.make();
  }

  /**
   * Map interval to screen coordinates using the provided ruler.
   *
   * @param rulerMap List typically obtained from Ruler_TO_BE_DEPRECTED.mapping of length equal to
   *     the screen. width mapping genome coords to screen coords.
   */
  public void mapToScreen(List<Double> rulerMap) {

    int xfrom = this.getFrom();
    int xto = this.getTo();

    // For INDELS where the first base of the ALT allele is the same as the REF allele,
    // we bump the start position by 1 so that 1bp indels are shown as 1 character rather than 2.
    if (this.variantContext.getReference().getBases().length > 0
        && this.variantContext.getAlleles().size() > 1
        && this.variantContext.getAlleles().get(1).getBases().length > 0
        && this.variantContext.getReference().getBases()[0]
            == this.variantContext.getAlleles().get(1).getBases()[0]) {
      xfrom = xfrom + 1;
      if (xfrom
          > xto) { // This can happen for multiallelic. See docstring this.variantContext.getEnd()
        xto = xfrom;
      }
    }
    /*        |============| <- ruler
     *   ===                  ===  <- Interval(s)
     */
    if ((xfrom < rulerMap.get(0) && xto < rulerMap.get(0))
        || (xfrom > rulerMap.get(rulerMap.size() - 1)) && xto > rulerMap.get(rulerMap.size() - 1)) {
      this.setScreenFrom(-1);
      this.setScreenTo(-1);
      return;
    }

    /*
     * Feature fully contains ruler map?
     *        |============| <- ruler
     *   ===================== <- Interval
     */
    if (xfrom <= rulerMap.get(0) && this.getTo() >= rulerMap.get(rulerMap.size() - 1)) {
      this.setScreenFrom(0);
      this.setScreenTo(rulerMap.size() - 1);
      return;
    }

    // Feature is all or partially contained
    this.setScreenFrom(Utils.getIndexOfclosestValue(xfrom, rulerMap));
    this.setScreenTo(Utils.getIndexOfclosestValue(xto, rulerMap));
    /*        |============|      <- ruler
     *   ========   ===    =====  <- Interval(s)
     */
    if (this.getScreenFrom() == -1) {
      this.setScreenFrom(0);
    }
    if (this.getScreenTo() == -1) {
      this.setScreenTo(rulerMap.size() - 1);
    }
    if (this.getScreenTo() == -1) {
      throw new RuntimeException("Unexpected mapping of features to ruler.");
    }
  }

  private char getCharForVCFIdeogram() {

    if (this.variantContext.getAlleles().size() > 2) { // Multiallelic
      return '|';
    } else if (this.variantContext.isSNP()) {
      Allele alt = this.variantContext.getAlleles().get(1);
      return alt.getBaseString().charAt(0);
    } else if (this.variantContext.isSimpleInsertion()) {
      return 'I';
    } else if (this.variantContext.isSimpleDeletion()) {
      return 'D';
    } else if (this.variantContext.isMNP()) {
      return 'M';
    } else if (this.variantContext.isComplexIndel()) {
      return 'X';
    } else {
      return '|';
    }
  }

  @Override
  protected void makeIdeogram(boolean addName) {
    int ideogramLength = this.getScreenTo() - this.getScreenFrom() + 1;
    this.setIdeogram(new ArrayList<>(ideogramLength), false);

    for (int i = 0; i < ideogramLength; i++) {
      FeatureChar c = new FeatureChar();
      char ideogramChar;
      ideogramChar = this.getCharForVCFIdeogram();
      c.addFormatVCF(ideogramChar);
      c.setText(ideogramChar);
      this.getIdeogram(false, false).add(c);
    }

    if (addName) {
      this.addNameToIdeogram();
    }
  }

  protected VariantContext getVariantContext() {
    return this.variantContext;
  }

  private int setFromForVCF() {
    int from = this.variantContext.getStart();
    this.setToForVCF();
    return from;
  }

  private int setToForVCF() {
    if (this.variantContext.getAlleles().size() > 2) { // Multiallelic
      return this.variantContext.getEnd();
    } else if (this.variantContext.isSNP()) {
      return this.variantContext.getStart();
    } else if (this.variantContext.isSimpleInsertion()) {
      int alt_len = this.variantContext.getAlleles().get(1).length();
      return this.variantContext.getStart() + alt_len - 1;
    } else if (this.variantContext.isSimpleDeletion() || this.variantContext.isMNP()) {
      int ref_len = this.variantContext.getAlleles().get(0).length();
      return this.variantContext.getStart() + ref_len - 1;
    } else if (this.variantContext.isComplexIndel()) {
      int ref_len = this.variantContext.getAlleles().get(0).length();
      int alt_len = this.variantContext.getAlleles().get(1).length();
      if (ref_len > alt_len) { // Similar to a deletion
        return this.variantContext.getStart() + ref_len - 1;
      } else {
        // Similar to an insertion
        return this.variantContext.getStart() + alt_len - 1;
      }
    } else {
      return this.variantContext.getEnd();
    }
  }
}
