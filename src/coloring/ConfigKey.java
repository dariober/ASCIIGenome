package coloring;

import java.util.HashSet;
import java.util.Set;

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
	shade_low_mapq("Colour for shading reads wit low MAPQ"),
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
	shade_structural_variant("Background colour for reads suggesting structural variation");
	
	private String value;

	ConfigKey(String value) {
        this.value = value;
    }

	/** 			U P D A T E   M E 
	 * Configuration parameters that are NOT COLOURS
	 * */
	public static Set<ConfigKey> nonColorKeys(){
		Set<ConfigKey> nonColorKeys= new HashSet<ConfigKey>(); 
		nonColorKeys.add(ConfigKey.max_reads_in_stack);
		nonColorKeys.add(ConfigKey.shade_baseq);
		return nonColorKeys;
	}
	
    public static ConfigKey getEnum(String value) {
        if(value == null)
            throw new IllegalArgumentException();
        for(ConfigKey v : values())
            if(value.equalsIgnoreCase(v.getDescription())) return v;
        throw new IllegalArgumentException();
    }

    public String getDescription() {
        return value;
    }
	
	public static Set<String> getValues(){
		HashSet<String> values = new HashSet<String>();

		  for (ConfigKey c : ConfigKey.values()) {
		      values.add(c.name());
		  }
		  return values;
	}
	
}
