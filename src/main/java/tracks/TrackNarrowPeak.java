package tracks;

import java.io.IOException;
import java.sql.SQLException;

import exceptions.InvalidGenomicCoordsException;
import exceptions.InvalidRecordException;
import samTextViewer.GenomicCoords;
import samTextViewer.Utils;

public class TrackNarrowPeak extends TrackWiggles {

	public TrackNarrowPeak(String filename, GenomicCoords gc) throws IOException, InvalidGenomicCoordsException,
			ClassNotFoundException, InvalidRecordException, SQLException {

		this.setTrackFormat(TrackFormat.NARROWPEAK);
		this.setFilename(filename);
		this.setWorkFilename(filename);
		this.bdgDataColIdx= 7;
		
		if(! Utils.hasTabixIndex(filename)){
			String tabixBdg= this.tabixBedgraphToTmpFile(filename);
			this.setWorkFilename(tabixBdg);
		}
		this.setGc(gc);
	}

	@Override
	public void update() throws IOException, InvalidRecordException, InvalidGenomicCoordsException, ClassNotFoundException, SQLException {

		if(this.bdgDataColIdx < 5){
			System.err.println("Invalid index for narrowPeak column of data value. Expected >= 5. Got " + this.bdgDataColIdx);
			this.bdgDataColIdx= 7;
			throw new InvalidRecordException();
		}
		this.bedGraphToScores(this.getWorkFilename());
	}	
}
