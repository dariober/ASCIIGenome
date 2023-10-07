package tracks;

import java.util.HashSet;
import java.util.Set;

import org.broad.igv.bbfile.BBFileReader;
import htsjdk.tribble.readers.TabixReader;

/** Adapter to make tabixReader and bigBed behave in the same way.
 * */
public class TabixBigBedReader {

    private TabixReader tabixReader; 
    private BBFileReader bigBedReader;
    
    protected TabixBigBedReader(TabixReader tabixReader){
        this.tabixReader= tabixReader;
    };
    
    protected TabixBigBedReader(BBFileReader bigBedReader){
        this.bigBedReader = bigBedReader;
    };
    
    protected TabixBigBedIterator query(String chrom, int start, int end){

        if(this.tabixReader != null){
            return new TabixBigBedIterator(this.tabixReader, chrom, start, end);
        
        } else if(this.bigBedReader != null){
            return new TabixBigBedIterator(this.bigBedReader, chrom, start, end);
        
        } else {
            throw new RuntimeException();
        }
    }

    public Set<String> getChromosomes() {

        if(this.tabixReader != null && this.bigBedReader == null){
        return this.tabixReader.getChromosomes();		
    
        } else if(this.tabixReader == null && this.bigBedReader != null){
            return new HashSet<String>(this.bigBedReader.getChromosomeNames());
            
        } else{
            System.err.println("Tabix and bigBed readers cannto be both set. One must be null.");
            throw new RuntimeException();
        }


    }
    
}
