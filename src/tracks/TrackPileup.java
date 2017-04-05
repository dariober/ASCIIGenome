package tracks;

import java.io.IOException;
import java.net.MalformedURLException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.lang3.StringUtils;

import exceptions.InvalidGenomicCoordsException;
import exceptions.InvalidRecordException;
import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.filter.AggregateFilter;
import samTextViewer.GenomicCoords;
import samTextViewer.Utils;

/** Collect depth and coverage info over a region of a chromosome
 * */
class TrackPileup extends TrackWiggles {

	/** Key: Position of the locus. Value: information at this position
	 * */
	private Map<Integer, Locus> loci= new LinkedHashMap<Integer, Locus>();
	private long nRecsInWindow;
	
	/*        C O N S T R U C T O R         */

	/** Initialize pileup with a list of empty loci at the region given by chrom:from-to
	 * At the start, the loci map is empty as no information is provided yet. 
	 * Positions without info (i.e. 0 coverage) are not present in the map at all 
	 * (to save memory and cpu). 
	 * @throws IOException 
	 * @throws SQLException 
	 * @throws InvalidRecordException 
	 * @throws InvalidGenomicCoordsException 
	 * @throws ClassNotFoundException 
	 * */
	protected TrackPileup(String bam, GenomicCoords gc) throws IOException, ClassNotFoundException, InvalidGenomicCoordsException, InvalidRecordException, SQLException {

		if(!Utils.bamHasIndex(bam)){
			System.err.println("\nAlignment file " + bam + " has no index.\n");
			throw new RuntimeException();
		}
		this.setFilename(bam);
		this.setWorkFilename(bam);
		this.setGc(gc);	
	}

	/*       M E T H O D S        */
	
	@Override
	public void update() throws MalformedURLException{
		
		SamReader samReader= Utils.getSamReader(this.getWorkFilename());
		
		this.nRecsInWindow= Utils.countReadsInWindow(this.getWorkFilename(), this.getGc(), this.getSamRecordFilter());
		
		Iterator<SAMRecord> sam= samReader.query(this.getGc().getChrom(), this.getGc().getFrom(), this.getGc().getTo(), false);
		AggregateFilter aggregateFilter= new AggregateFilter(this.getSamRecordFilter());

		while(sam.hasNext()){
			
			SAMRecord rec= sam.next();

			Boolean passAwk= true; // Utils.passAwkFilter(rec.getSAMString(), this.getAwk());
			
			if( !rec.getReadUnmappedFlag() && !aggregateFilter.filterOut(rec) && passAwk){
				this.add(rec);
			}
		}

	}
	
	/** Update the map of loci with the information in this record. 
	 * */
	protected void add(SAMRecord samRecord){
		
		// Is this read forward or reverse? First or second in pair?
		boolean isFirstOFPair= samRecord.getFirstOfPairFlag();
		boolean isReverse= samRecord.getReadNegativeStrandFlag();
		
		List<AlignmentBlock> alnBlocks = samRecord.getAlignmentBlocks();
		if(alnBlocks.size() == 0){
			// Nothing to be done. This may happen with e.g. fully clipped reads
			return; 
		}
		int readPos= alnBlocks.get(0).getReadStart() - 1; // Leftmost position where the read starts aligning.
		for(AlignmentBlock block : alnBlocks){
			
			if(block.getReferenceStart() > this.getGc().getTo()){
				// This block is completely to the right of the user's coordinates. 
		        // No need to process this record anymore.
				return; 
			}
			if((block.getReferenceStart() + block.getLength() - 1) < this.getGc().getFrom()){
				// This block is completely to the left of the user's coordinates.
				// So no useful information here.
				continue;
			}
			
			int refPos= block.getReferenceStart() - 1; // -1 because we increment immediately in the loop
			for(int i= 0; i < block.getLength(); i++){
				// Iterate through each base of this block
				refPos++;
				readPos++;
				if(refPos < this.getGc().getFrom()){
					continue; // Skip this position. Part of this block is outside user's coordinates. 
				}
				if(refPos > this.getGc().getTo()){
					return; // We have passed the user's endpoint, no need to process this record anymore
				}
				
				// Start collecting info	
				if( ! this.loci.containsKey(refPos)){
					// Add this position to the map.
					this.loci.put(refPos, new Locus(this.getGc().getChrom(), refPos));
				}

				// What read base do we have at this position?
				char base= (char) samRecord.getReadBases()[readPos-1];				
				this.loci.get(refPos).add(base, isReverse, isFirstOFPair);
			}
		}
		// Now we need to increment counts corresponding to deletions in the reference
		List<int[]>deletedBlocks= this.getRefPositionOfDeletedBlocks(samRecord);
		for(int[] block : deletedBlocks){
			for(int refPos= block[0]; refPos <= block[1]; refPos++){
				if( ! this.loci.containsKey(refPos)){
					// Add this position to the map.
					this.loci.put(refPos, new Locus(this.getGc().getChrom(), refPos));
				}
				this.loci.get(refPos).add('D', isReverse, isFirstOFPair);
			}
		}
	}

	/**
	 * MEMO: Deletion does not consume read bases. It consumes reference bases:
	 * ref  NNNNNNNN
	 * read NNN---NN
	 * */
	private List<int[]> getRefPositionOfDeletedBlocks(SAMRecord samRecord){
		
		int ndel= StringUtils.countMatches(samRecord.getCigarString().toUpperCase(), 'D');
		List<int[]> deletedBlocks= new ArrayList<int[]>(ndel); // int array contains start and end of deletion
		if( ndel == 0){
			return deletedBlocks; // No deletetions, nothing to be done.
		}
		
		// 1S 2M 3D 4M 5D 6M 7S
		// Walk along the cigar string. When an operator consumes read bases, advance the position tracker
		// When you hit a deletion record the start and end of the deleted block in ref coordinates.
		List<CigarElement> cigarOps = samRecord.getCigar().getCigarElements();
		int readTracker= 0;
		int[] deletedBlock= new int[2]; // Length 2 becouse it contains start and end of deleted block
		for(CigarElement op : cigarOps){
			if(op.getOperator().equals(CigarOperator.DELETION)){
				deletedBlock[0]= samRecord.getReferencePositionAtReadPosition(readTracker) + 1;
				deletedBlock[1]= samRecord.getReferencePositionAtReadPosition(readTracker + 1) - 1;
				deletedBlocks.add(deletedBlock);
				if(deletedBlocks.size() == ndel){
					// We have accumulated all the deleted blocks. No need to process this read any more.
					return deletedBlocks; 
				}
			}
			if(op.getOperator().consumesReadBases()){
				readTracker += op.getLength();
			}
		}
		return deletedBlocks;
	}

	protected Map<Integer, Integer> getDepth(){
		Map<Integer, Integer> depth= new LinkedHashMap<Integer, Integer>();
		for(int pos : this.loci.keySet()){
			Locus loc= this.loci.get(pos);
			depth.put(pos, loc.getDepth());
		}
		return depth;
	}
	
	protected Map<Integer, Locus> getLoci(){
		return this.loci;
	}
	
}
