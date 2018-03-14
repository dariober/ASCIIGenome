package filter;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.filter.SamRecordFilter;

public class NotPrimaryAlignmentFilter implements SamRecordFilter {
    
	private boolean include= false;
	
	/**
     * @param record the SAMRecord to evaluate
     * @return true if the SAMRecord matches the filter, otherwise false
     */
    public NotPrimaryAlignmentFilter(final boolean include) {
        this.include= include;
    }

    /**
     * Determines whether a SAMRecord matches this filter
     *
     * @param record the SAMRecord to evaluate
     *
     * @return true if the SAMRecord matches the filter, otherwise false
     */
    public boolean filterOut(final SAMRecord record) {
        if (include) {
            if (record.getNotPrimaryAlignmentFlag()) {
                return false;
            }
        } else {
            // exclude
            if (!record.getNotPrimaryAlignmentFlag()) {
                return false;
            }
        }

        return true;
    }

    /**
     * Determines whether a pair of SAMRecord matches this filter
     *
     * @param first  the first SAMRecord to evaluate
     * @param second the second SAMRecord to evaluate
     *
     * @return true if the SAMRecords matches the filter, otherwise false
     */
    public boolean filterOut(final SAMRecord first, final SAMRecord second) {
        throw new UnsupportedOperationException("Paired *Filter not implemented!");
    }
}
