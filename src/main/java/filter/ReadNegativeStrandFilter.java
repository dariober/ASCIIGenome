package filter;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.filter.SamRecordFilter;

public class ReadNegativeStrandFilter implements SamRecordFilter {

  private boolean include = false;

  public ReadNegativeStrandFilter(final boolean include) {
    this.include = include;
  }

  @Override
  public boolean filterOut(final SAMRecord record) {
    if (include) {
      if (record.getReadNegativeStrandFlag()) {
        return false;
      }
    } else {
      // exclude
      if (!record.getReadNegativeStrandFlag()) {
        return false;
      }
    }

    return true;
  }

  @Override
  public boolean filterOut(final SAMRecord first, final SAMRecord second) {
    throw new UnsupportedOperationException("Paired *Filter not implemented!");
  }
}
