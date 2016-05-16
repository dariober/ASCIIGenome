package filter;

import static org.junit.Assert.*;

import java.io.File;
import java.util.List;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.filter.AlignedFilter;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.SamLocusIterator;
import htsjdk.samtools.util.SamLocusIterator.LocusInfo;
import htsjdk.samtools.util.SamLocusIterator.RecordAndOffset;

import org.junit.Test;

public class FilterTest {

	final SAMRecord rec= new SAMRecord(null);
	
	@Test
	public void exampleFromSamtools(){
		// How filters behave:
		// filterOut(): 
		// Return true if the SAMRecord matches the filter, otherwise false
		
		// Example: Keep aligned reads
		rec.setFlags(0);

		AlignedFilter f= new AlignedFilter(true); // Set true because we require aligned
		assertFalse(f.filterOut(rec)); // filterOut() returns false because with flag 0 it's not filtered out.

		f= new AlignedFilter(false);
		assertTrue(f.filterOut(rec)); // Exclude if aligned

		
		// Keep non aligned reads (flag 4)
		rec.setFlags(4);
		f= new AlignedFilter(false); // Keep if non aligned
		assertFalse(f.filterOut(rec)); // Returns false because read is kept.
		
		f= new AlignedFilter(true); // Keep if aligned
		assertTrue(f.filterOut(rec)); // Exclude read because 4 is non aligned
	}
	
	@Test
	public void testFilters(){

		rec.setFlags(1);
		ReadPairedFilter f= new ReadPairedFilter(true); // Keep if paired
		assertFalse(f.filterOut(rec));

		f= new ReadPairedFilter(false); // Exclude if paired
		assertTrue(f.filterOut(rec));

		rec.setFlags(1+2);
		assertFalse(new ProperPairFilter(true).filterOut(rec)); // Keep
		assertTrue(new ProperPairFilter(false).filterOut(rec)); // Exclude

		rec.setFlags(1);
		assertTrue(new ProperPairFilter(true).filterOut(rec)); // Exclude becasue paired but proper

		rec.setFlags(1+8);
		assertFalse(new MateUnmappedFilter(true).filterOut(rec)); // Keep if mate unmapped
		assertTrue(new MateUnmappedFilter(false).filterOut(rec)); // Exclude

		rec.setFlags(1);
		assertTrue(new MateUnmappedFilter(true).filterOut(rec)); // Exclude because mate mapped

		rec.setFlags(16);
		assertFalse(new ReadNegativeStrandFilter(true).filterOut(rec)); // Keep if -ve: Kept
		assertTrue(new ReadNegativeStrandFilter(false).filterOut(rec)); // Exclude if -ve: Excluded

		rec.setFlags(0);
		assertFalse(new ReadNegativeStrandFilter(false).filterOut(rec)); // Keep if +ve: Kept
		rec.setFlags(16);
		assertTrue(new ReadNegativeStrandFilter(false).filterOut(rec)); // Kept if +ve: Excluded

		rec.setFlags(1 + 32);
		assertFalse(new MateNegativeStrandFilter(true).filterOut(rec)); // Keep if set: Kept
		assertTrue(new MateNegativeStrandFilter(false).filterOut(rec)); // Exclude if set: Excluded

		rec.setFlags(1 + 64);
		assertFalse(new FirstOfPairFilter(true).filterOut(rec)); // Keep if set: Kept
		assertTrue(new FirstOfPairFilter(false).filterOut(rec)); // Exclude if set: Excluded

		rec.setFlags(1 + 128);
		assertFalse(new SecondOfPairFilter(true).filterOut(rec)); // Keep if set: Kept
		assertTrue(new SecondOfPairFilter(false).filterOut(rec)); // Exclude if set: Excluded

		rec.setFlags(256);
		assertFalse(new NotPrimaryAlignmentFilter(true).filterOut(rec)); // Keep if set: Kept
		assertTrue(new NotPrimaryAlignmentFilter(false).filterOut(rec)); // Exclude if set: Excluded

		rec.setFlags(512);
		assertFalse(new ReadFailsVendorQualityCheckFilter(true).filterOut(rec)); // Keep if set: Kept
		assertTrue(new ReadFailsVendorQualityCheckFilter(false).filterOut(rec)); // Exclude if set: Excluded

		rec.setFlags(1024);
		assertFalse(new DuplicateReadFilter(true).filterOut(rec)); // Keep if set: Kept
		assertTrue(new DuplicateReadFilter(false).filterOut(rec)); // Exclude if set: Excluded

		rec.setFlags(2048);
		assertFalse(new SupplementaryAlignmentFilter(true).filterOut(rec)); // Keep if set: Kept
		assertTrue(new SupplementaryAlignmentFilter(false).filterOut(rec)); // Exclude if set: Excluded		
	}
	
	@Test
	public void filterForTopBottomStrand(){
		rec.setFlags(0);
		assertFalse(new ReadFromTopStrandFilter(true).filterOut(rec)); // Keep if set: Kept
		assertTrue(new ReadFromTopStrandFilter(false).filterOut(rec)); // Exclude if set: Excluded
		rec.setFlags(16);
		assertTrue(new ReadFromTopStrandFilter(true).filterOut(rec)); // Keep if set: Excluded
		assertFalse(new ReadFromTopStrandFilter(false).filterOut(rec)); // Exclude if set: Kept
		
		rec.setFlags(17); // paired -ve strand
		// Keep if from bottom strand
		assertFalse(new ReadFromTopStrandFilter(false).filterOut(rec));
		// Keep if from top strand
		assertTrue(new ReadFromTopStrandFilter(true).filterOut(rec));
	}
	
	@Test
	public void canFilterFromIntFlag(){
		
		int f_incl= 131;
		int F_excl= 72;
		String chrom= "chrY";
		int from= 1;
		int to= 100;
		
		SamReaderFactory srf=SamReaderFactory.make();
		SamReader samReader= srf.open(new File("test_data/mjb050_oxBS.bam"));
		SAMFileHeader fh= samReader.getFileHeader();
		IntervalList il= new IntervalList(fh);
		Interval interval= new Interval(chrom, from, to);
		il.add(interval);
		
		List<SamRecordFilter> filters= FlagToFilter.flagToFilterList(f_incl, F_excl);
	
		SamLocusIterator samLocIter= new SamLocusIterator(samReader, il, true);
		samLocIter.setSamFilters(filters);
		while(samLocIter.hasNext()){
			LocusInfo locus= samLocIter.next();
			if(locus.getRecordAndPositions().size() > 0){
				//System.out.println(locus.getPosition() + " " + locus.getRecordAndPositions().size() + " " +
				//		locus.getRecordAndPositions().get(0).getRecord().getFlags());
			}
		}
	}
	
	@Test
	public void STUBcanCallBS(){
		
		System.out.println("\u203E");
		
		int f_incl= 0;
		int F_excl= 0;
		String chrom= "chrY";
		int from= 1;
		int to= 1;
		
		SamReaderFactory srf=SamReaderFactory.make();
		SamReader samReader= srf.open(new File("test_data/mjb050_oxBS.bam"));
		
		SAMFileHeader fh= samReader.getFileHeader();
		IntervalList il= new IntervalList(fh);
		Interval interval= new Interval(chrom, from, to);
		il.add(interval);
		
		List<SamRecordFilter> filters= FlagToFilter.flagToFilterList(f_incl, F_excl);
	
		SamLocusIterator samLocIter= new SamLocusIterator(samReader, il, true);
		samLocIter.setSamFilters(filters);
		while(samLocIter.hasNext()){
			LocusInfo locus= samLocIter.next();
			int M= 0;
			int U= 0;
			int mism= 0;
			for(RecordAndOffset recOff : locus.getRecordAndPositions()){
				int pos= locus.getPosition();
				// Code to get ref sequence at pos
				
				// If ref sequence is C, count read bases if:
				char refbase= Character.toUpperCase('\0');
				char rb= Character.toUpperCase((char)recOff.getReadBase());
				
				boolean isTopStrand= (new ReadFromTopStrandFilter(true)).filterOut(recOff.getRecord());
				
				if(refbase == 'C'){

					if( isTopStrand	){ // -ve 2nd pair
						if(rb == 'C'){
							M++;
						} else if(rb == 'T'){
							U++;
						} else {
							mism++;
						}
					}  					
				} else if (refbase == 'G'){

					if(	!isTopStrand ){
						if(rb == 'G'){
							M++;
						} else if(rb == 'A'){
							U++;
						} else {
							mism++;
						}
							// System.out.println(locus.getPosition() + " ");					
						}  
				} else {
					// Not a C on ref
				}

			}
			
//			if(locus.getRecordAndPositions().size() > 0){
//				System.out.println(locus.getPosition() + " " + locus.getRecordAndPositions().size() + " " +
//						locus.getRecordAndPositions().get(0).getRecord().getFlags());
//			}
		}
	}
}

/*
// READ_PAIRED_FLAG = 0x1;
// PROPER_PAIR_FLAG = 0x2;
// READ_UNMAPPED_FLAG = 0x4; <- Implemented via AlignedFilter(false)
// MATE_UNMAPPED_FLAG = 0x8;
// READ_STRAND_FLAG = 0x10;
// MATE_STRAND_FLAG = 0x20;
// FIRST_OF_PAIR_FLAG = 0x40;
// SECOND_OF_PAIR_FLAG = 0x80;
// NOT_PRIMARY_ALIGNMENT_FLAG = 0x100;
// READ_FAILS_VENDOR_QUALITY_CHECK_FLAG = 0x200;
// DUPLICATE_READ_FLAG = 0x400;
// SUPPLEMENTARY_ALIGNMENT_FLAG = 0x800;
*/

/*
 		String sam= "test_data/ds051.actb.bam";
		String chrom= "chr7";
		int from= 5566781;
		int to= from+10000;
		int windowSize= 50;
		int mapq= 50;

		SamReader samReader= ReadWriteBAMUtils.reader(sam, ValidationStringency.SILENT);
		SAMFileHeader fh= samReader.getFileHeader();
		
		IntervalList il= new IntervalList(fh);
		Interval interval= new Interval(chrom, from, to);
		il.add(interval);
		
		List<SamRecordFilter> filters= new ArrayList<SamRecordFilter>();
		//filters.add(new FailsVendorReadQualityFilter());
		//filters.add(new NotPrimaryAlignmentFilter());
		//filters.add(new DuplicateReadFilter());
		//filters.add(new AlignedFilter(true));
		//filters.add(new MappingQualityFilter(mapq));
		
	
		//ReadPairedFilter f= new ReadPairedFilter(false);
		//System.out.println(f.filterOut(rec));
		
		SamLocusIterator samLocIter= new SamLocusIterator(samReader, il, true);
		samLocIter.setSamFilters(filters);
		while(samLocIter.hasNext()){
			int locus= samLocIter.next().getRecordAndPositions().size();			
			System.out.println(locus);
		}
		//		
 */