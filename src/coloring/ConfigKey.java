package coloring;

import java.util.HashSet;
import java.util.Set;

public enum ConfigKey {
	// Keep all these fields lower case.
	background,
	foreground,
	read_negative_strand,
	read_positive_strand,
	seq_a,
	seq_c,
	seq_g,
	seq_t,
	seq_other,
	shade_low_mapq,
	methylated_foreground,
	unmethylated_foreground,
	methylated_background,
	unmethylated_background,
	read_default_background,
	read_default_foreground,
	title_colour,
	feature_foreground,
	feature_background_positive_strand,
	feature_background_negative_strand,
	feature_background_no_strand,
	footer,
	chrom_ideogram,
	ruler;
	
	public static Set<String> getValues(){
		HashSet<String> values = new HashSet<String>();

		  for (ConfigKey c : ConfigKey.values()) {
		      values.add(c.name());
		  }
		  return values;
	}
	
}
