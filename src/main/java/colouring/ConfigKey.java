package colouring;

import com.google.common.base.Joiner;
import exceptions.InvalidConfigException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import samTextViewer.Utils;

// See http://stackoverflow.com/questions/15989316/how-to-add-a-description-for-each-entry-of-enum
// for how this class was prepared.
public enum ConfigKey {
  // Keep all these fields lower case.
  background("Background colour"),
  foreground("Foreground colour"),
  seq_a("Colour for nucleotide A"),
  seq_c("Colour for nucleotide C"),
  seq_g("Colour for nucleotide G"),
  seq_t("Colour for nucleotide T"),
  seq_other("Colour for any other nucleotide"),
  shade_low_mapq("Colour for shading reads with low MAPQ"),
  low_mapq("Shade reads below this MAPQ"),
  methylated_foreground("Foreground colour for methylated C"),
  unmethylated_foreground("Foreground colour for unmethylated C"),
  methylated_background("Background colour for methylated C"),
  unmethylated_background("Background colour for unmethylated C"),
  title_colour("Default Colour for titles"),
  feature_background_positive_strand("Colour for features on forward strand"),
  feature_background_negative_strand("Colour for features on reverse strand"),
  feature_background_no_strand("Colour for features without strand information"),
  footer("Colour for footer line"),
  chrom_ideogram("Colour for chromosome ideogram"),
  ruler("Colour for ruler"),
  max_reads_in_stack("Max number of reads to accumulate when showing read tracks"),
  shade_baseq("Shade read base when quality is below this threshold"),
  shade_structural_variant(
      "Background colour for reads suggesting structural variation or 'false' for no shading"),
  highlight_mid_char("Highlight mid-character in read tracks?"),
  nucs_as_letters("Show read nucleotides as letters at single base resolution?"),
  show_soft_clip("NOT IN USE YET - Show soft clipped bases in read tracks?"),
  stop_codon("Colour for stop codon"),
  start_codon("Colour for start codon"),
  codon("Colour for codons other than start and stop");

  private String value;

  ConfigKey(String value) {
    this.value = value;
  }

  /** U P D A T E M E Configuration parameters that are NOT COLOURS */
  public static Set<ConfigKey> colourKeys() {
    Set<ConfigKey> colourKeys = new HashSet<ConfigKey>();
    colourKeys.add(ConfigKey.background);
    colourKeys.add(ConfigKey.foreground);
    colourKeys.add(ConfigKey.seq_a);
    colourKeys.add(ConfigKey.seq_c);
    colourKeys.add(ConfigKey.seq_g);
    colourKeys.add(ConfigKey.seq_t);
    colourKeys.add(ConfigKey.seq_other);
    colourKeys.add(ConfigKey.shade_low_mapq);
    colourKeys.add(ConfigKey.methylated_foreground);
    colourKeys.add(ConfigKey.unmethylated_foreground);
    colourKeys.add(ConfigKey.methylated_background);
    colourKeys.add(ConfigKey.unmethylated_background);
    colourKeys.add(ConfigKey.title_colour);
    colourKeys.add(ConfigKey.feature_background_positive_strand);
    colourKeys.add(ConfigKey.feature_background_negative_strand);
    colourKeys.add(ConfigKey.feature_background_no_strand);
    colourKeys.add(ConfigKey.footer);
    colourKeys.add(ConfigKey.chrom_ideogram);
    colourKeys.add(ConfigKey.ruler);
    colourKeys.add(ConfigKey.shade_structural_variant);
    colourKeys.add(ConfigKey.start_codon);
    colourKeys.add(ConfigKey.stop_codon);
    colourKeys.add(ConfigKey.codon);
    return colourKeys;
  }

  /** U P D A T E M E */
  public static Set<ConfigKey> booleanKeys() {
    Set<ConfigKey> booleanKeys = new HashSet<ConfigKey>();
    booleanKeys.add(ConfigKey.highlight_mid_char);
    booleanKeys.add(ConfigKey.nucs_as_letters);
    booleanKeys.add(ConfigKey.show_soft_clip);
    return booleanKeys;
  }

  public static Set<ConfigKey> integerKeys() {
    Set<ConfigKey> integerKeys = new HashSet<ConfigKey>();
    integerKeys.add(ConfigKey.max_reads_in_stack);
    integerKeys.add(ConfigKey.shade_baseq);
    integerKeys.add(ConfigKey.low_mapq);
    return integerKeys;
  }

  public static ConfigKey getEnum(String value) {
    if (value == null) throw new IllegalArgumentException();
    for (ConfigKey v : values()) if (value.equalsIgnoreCase(v.getDescription())) return v;
    throw new IllegalArgumentException();
  }

  public String getDescription() {
    return value;
  }

  public static Set<String> getValues() {
    HashSet<String> values = new HashSet<String>();

    for (ConfigKey c : ConfigKey.values()) {
      values.add(c.name());
    }
    return values;
  }

  public static ConfigKey getConfigKeyFromShort(String configkey)
      throws InvalidConfigException, IOException {
    List<String> optionKeys = new ArrayList<String>(ConfigKey.getValues());
    ConfigKey key;
    if (optionKeys.contains(configkey)) {
      key = ConfigKey.valueOf(configkey);
    } else {
      List<String> candidates = Utils.suggestCommand(configkey, optionKeys);
      if (candidates.size() == 0) {
        System.err.println(
            Utils.padEndMultiLine(
                "No config key found for " + configkey, Utils.getTerminalWidth()));
        throw new InvalidConfigException();
      }
      if (candidates.size() > 1) {
        System.err.println(
            Utils.padEndMultiLine(
                "Did you mean " + Joiner.on(" or ").join(candidates) + "?",
                Utils.getTerminalWidth()));
        throw new InvalidConfigException();
      }
      String candidate = candidates.get(0);
      boolean found = false;
      for (String x : optionKeys) {
        if (x.toLowerCase().contains(configkey.toLowerCase())) {
          found = true;
          break;
        }
      }
      if (!found) {
        System.err.println(
            Utils.padEndMultiLine(
                "No config key found for " + configkey, Utils.getTerminalWidth()));
        throw new InvalidConfigException();
      }
      System.err.println(Utils.padEndMultiLine("Setting " + candidate, Utils.getTerminalWidth()));
      key = ConfigKey.valueOf(candidate);
    }
    return key;
  }
}
