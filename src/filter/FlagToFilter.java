package filter;

import htsjdk.samtools.filter.AlignedFilter;
import htsjdk.samtools.filter.SamRecordFilter;

import java.util.ArrayList;
import java.util.List;


public class FlagToFilter {
	
	public static List<SamRecordFilter> flagToFilterList(int f_incl, int F_excl){
		
		List<SamRecordFilter> list= new ArrayList<SamRecordFilter>();
		
		if((f_incl & 1) == 1){
			list.add(new ReadPairedFilter(true));
		}
		if((F_excl & 1) == 1){
			list.add(new ReadPairedFilter(false));
		}
		
		if((f_incl & 2) == 2){
			list.add(new ProperPairFilter(true));
		}
		if((F_excl & 2) == 2){
			list.add(new ProperPairFilter(false));
		}
		
		if((f_incl & 4) == 4){ // sic! tru/false is the other way round.
			list.add(new AlignedFilter(false));
		}
		if((F_excl & 4) == 4){
			list.add(new AlignedFilter(true));
		}

		if((f_incl & 8) == 8){
			list.add(new MateUnmappedFilter(true));
		}
		if((F_excl & 8) == 8){
			list.add(new MateUnmappedFilter(false));
		}

		if((f_incl & 16) == 16){
			list.add(new ReadNegativeStrandFilter(true));
		}
		if((F_excl & 16) == 16){
			list.add(new ReadNegativeStrandFilter(false));
		}
		
		if((f_incl & 32) == 32){
			list.add(new MateNegativeStrandFilter(true));
		}
		if((F_excl & 32) == 32){
			list.add(new MateNegativeStrandFilter(false));
		}		

		if((f_incl & 64) == 64){
			list.add(new FirstOfPairFilter(true));
		}
		if((F_excl & 64) == 64){
			list.add(new FirstOfPairFilter(false));
		}		
		
		if((f_incl & 128) == 128){
			list.add(new SecondOfPairFilter(true));
		}
		if((F_excl & 128) == 128){
			list.add(new SecondOfPairFilter(false));
		}		

		if((f_incl & 256) == 256){
			list.add(new NotPrimaryAlignmentFilter(true));
		}
		if((F_excl & 256) == 256){
			list.add(new NotPrimaryAlignmentFilter(false));
		}		

		if((f_incl & 512) == 512){
			list.add(new ReadFailsVendorQualityCheckFilter(true));
		}
		if((F_excl & 512) == 512){
			list.add(new ReadFailsVendorQualityCheckFilter(false));
		}		

		if((f_incl & 1024) == 1024){
			list.add(new DuplicateReadFilter(true));
		}
		if((F_excl & 1024) == 1024){
			list.add(new DuplicateReadFilter(false));
		}		

		if((f_incl & 2048) == 2048){
			list.add(new SupplementaryAlignmentFilter(true));
		}
		if((F_excl & 2048) == 2048){
			list.add(new SupplementaryAlignmentFilter(false));
		}		

		if((f_incl & 4096) == 4096){
			list.add(new ReadFromTopStrandFilter(true));
		}
		if((F_excl & 4096) == 4096){
			list.add(new ReadFromTopStrandFilter(false));
		}		

		
		return list;
	} 

	//READ_PAIRED_FLAG = 0x1;
	//PROPER_PAIR_FLAG = 0x2;
	//READ_UNMAPPED_FLAG = 0x4; <- Implemented via AlignedFilter(false)
	//MATE_UNMAPPED_FLAG = 0x8;
	//READ_STRAND_FLAG = 0x10;
	//MATE_STRAND_FLAG = 0x20;
	//FIRST_OF_PAIR_FLAG = 0x40;
	//SECOND_OF_PAIR_FLAG = 0x80;
	//NOT_PRIMARY_ALIGNMENT_FLAG = 0x100;
	//READ_FAILS_VENDOR_QUALITY_CHECK_FLAG = 0x200;
	//DUPLICATE_READ_FLAG = 0x400;
	//SUPPLEMENTARY_ALIGNMENT_FLAG = 0x800;

	
}
