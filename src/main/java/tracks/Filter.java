package tracks;

public enum Filter {
  DEFAULT_HIDE_REGEX("^$"),
  DEFAULT_SHOW_REGEX(".*"),
  DEFAULT_AWK(""),
  DEFAULT_f_FLAG("0"),
  DEFAULT_F_FLAG("4"),
  DEFAULT_MAPQ("0"),
  DEFAULT_VARIANT_CHROM("");

  private String value;

  private Filter(String value) {
    this.value = value;
  }

  public String getValue() {
    return value;
  }
}
