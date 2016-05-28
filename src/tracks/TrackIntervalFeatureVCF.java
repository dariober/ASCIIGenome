package tracks;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import exceptions.InvalidGenomicCoordsException;
import samTextViewer.GenomicCoords;

public class TrackIntervalFeatureVCF extends TrackIntervalFeature {

	private List<IntervalFeature> intervalFeatureList= new ArrayList<IntervalFeature>();
	
	public TrackIntervalFeatureVCF(String filename, GenomicCoords gc)
			throws IOException, InvalidGenomicCoordsException {
		super(filename, gc);
		
		this.setGc(gc);
		this.setFilename(filename);
		this.intervalFeatureSet= new IntervalFeatureSet(filename);
		this.update();
		
	}

	/* M e t h o d s */
	@Override
	public void update() throws IOException, InvalidGenomicCoordsException{
		this.intervalFeatureList = this.intervalFeatureSet.getFeaturesInInterval(
				this.getGc().getChrom(), this.getGc().getFrom(), this.getGc().getTo());
						
		for(IntervalFeature ift : intervalFeatureList){
			ift= (IntervalFeatureVCF)ift;			
			ift.mapToScreen(this.getGc().getMapping());
		}
				
	}
	
}
