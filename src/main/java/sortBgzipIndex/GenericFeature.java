package sortBgzipIndex;

import htsjdk.tribble.Feature;

class GenericFeature implements Feature {
  private final String contig;
  private final int start;
  private final int end;

  public GenericFeature(String contig, int start, int end) {
    this.contig = contig;
    this.start = start;
    this.end = end;
  }

  public GenericFeature(String contig, int start) {
    this.contig = contig;
    this.start = start;
    this.end = -1;
  }

  @Override
  public String getContig() {
    return contig;
  }

  @Override
  public int getStart() {
    return start;
  }

  @Override
  public int getEnd() {
    return end;
  }
}
