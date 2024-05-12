package session;

import samTextViewer.GenomicCoords;
import tracks.TrackSet;

public class Session {
  private GenomicCoords genomicCoords;
  private TrackSet trackSet;

  public Session(String yamlFile) {
    // Read yaml file
    // Get genome file if set and parse it into genomicCoords
    // For each track populate trackSet

  }

  public TrackSet getTrackSet() {
    return this.trackSet;
  }

  public GenomicCoords getGenomicCoords() {
    return this.genomicCoords;
  }
}
