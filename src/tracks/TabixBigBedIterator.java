package tracks;

import java.io.IOException;

import org.broad.igv.bbfile.BBFileReader;
import org.broad.igv.bbfile.BedFeature;

import htsjdk.tribble.readers.TabixReader;
import htsjdk.tribble.readers.TabixReader.Iterator;
import org.broad.igv.bbfile.BigBedIterator;

public class TabixBigBedIterator {

	Iterator tabixIterator;
	BigBedIterator bigBedIterator;
	
	protected TabixBigBedIterator(TabixReader reader, String chrom, int start, int end){
		this.tabixIterator= reader.query(chrom, start, end);
	}

	protected TabixBigBedIterator(BBFileReader reader, String chrom, int start, int end){
		this.bigBedIterator= reader.getBigBedIterator(chrom, start, chrom, end, false);
	}
	
	protected String next() throws IOException{
		
		if(this.tabixIterator != null){
			return this.tabixIterator.next();

		} else if(this.bigBedIterator != null) {
			if( ! this.bigBedIterator.hasNext()){
				return null;
			}			
			BedFeature x= this.bigBedIterator.next();
			if(x == null){
				return null;
			}
			StringBuilder sb= new StringBuilder();
			sb.append(x.getChromosome());
			sb.append("\t");
			sb.append(x.getStartBase());
			sb.append("\t");
			sb.append(x.getEndBase());
			for(String field : x.getRestOfFields()){
				sb.append("\t");
				sb.append(field);
			}
			return sb.toString();
		} else {
			throw new RuntimeException();
		}
	}
	
}
