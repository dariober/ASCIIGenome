package filter;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.filter.SamRecordFilter;

/**
 * Extension of samflags to filter reads coming from top strand. Useful for
 * methylation calling. Assign to this filter the next bit integer (4096) 
 * @author berald01
 *
 */
public class ReadFromTopStrandFilter implements SamRecordFilter {
    
	private boolean include= false;
	
	/**
     * @param record the SAMRecord to evaluate
     * @return true if the SAMRecord matches the filter, otherwise false
     */
    public ReadFromTopStrandFilter(final boolean include) {
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
        
    	boolean isTopStrand= (
    			(!record.getReadNegativeStrandFlag() && !record.getReadPairedFlag()) ||  // +ve unpaired
		        (!record.getReadNegativeStrandFlag() && record.getReadPairedFlag() && record.getFirstOfPairFlag()) ||  // +ve 1st in pair
		        (record.getReadNegativeStrandFlag() && record.getReadPairedFlag() && record.getSecondOfPairFlag()));   // -ve 2nd in pair

//    	boolean isBottomStrand= (
//   			(record.getReadNegativeStrandFlag() && !record.getReadPairedFlag()) || 
//			    (record.getReadNegativeStrandFlag() && record.getFirstOfPairFlag()) || 
//					!record.getReadNegativeStrandFlag() && record.getSecondOfPairFlag()
//   			);
    	
    	if (include) {
            if ( isTopStrand ) {     
                return false;
            }
        } else {
            // exclude
            if ( !isTopStrand ) {
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
