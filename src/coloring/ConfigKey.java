package coloring;

import java.util.HashSet;
import java.util.Set;

public enum ConfigKey {
	// Keep all these fields lower case.
	background,
	foreground,
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
	title_colour,
	feature_background_positive_strand,
	feature_background_negative_strand,
	feature_background_no_strand,
	footer,
	chrom_ideogram,
	ruler,
	max_reads_in_stack;
	
	public static Set<String> getValues(){
		HashSet<String> values = new HashSet<String>();

		  for (ConfigKey c : ConfigKey.values()) {
		      values.add(c.name());
		  }
		  return values;
	}
	
}
