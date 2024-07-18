package session;

import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import tracks.Track;
import tracks.TrackSet;

public class Session {
  private String sessionName;
  private String lastRead;
  private SessionGenome genome;
  private Map<String, SessionTrack> tracks;

  public Session() {}

  public Session(SessionGenome sessionGenome, Map<String, SessionTrack> tracks) {
    LocalDateTime date = LocalDateTime.now();
    this.setLastRead(date.format(DateTimeFormatter.ISO_LOCAL_DATE_TIME));
    this.genome = sessionGenome;
    this.tracks = tracks;
  }

  public TrackSet toTrackSet() {
    List<Track> trackList = new ArrayList<>();
    for (String trackName : tracks.keySet()) {
      try {
        SessionTrack st = tracks.get(trackName);
        Track tr = st.toTrack(trackName, this.getGenome().toGenomicCoords());
        trackList.add(tr);
      } catch (Exception e) {
        //
      }
    }
    return new TrackSet(trackList);
  }

  public String getSessionName() {
    return sessionName;
  }

  public void setSessionName(String sessionName) {
    this.sessionName = sessionName;
  }

  public void setLastRead(String lastRead) {
    this.lastRead = lastRead;
  }

  public String getLastRead() {
    return lastRead;
  }

  public SessionGenome getGenome() {
    return this.genome;
  }

  public void setGenome(SessionGenome genome) {
    this.genome = genome;
  }

  public Map<String, SessionTrack> getTracks() {
    return this.tracks;
  }

  public void setTracks(Map<String, SessionTrack> track) {
    this.tracks = track;
  }

  @Override
  public boolean equals(Object o) {
    if (o == this) {
      return true;
    }
    if (!(o instanceof Session)) {
      return false;
    }
    Session s = (Session) o;
    return this.getSessionName().equals(s.getSessionName());
  }

  @Override
  public String toString() {
    return "Session{"
        + "genome="
        + genome
        + ", lastRead='"
        + lastRead
        + '\''
        + ", sessionName='"
        + sessionName
        + '\''
        + ", tracks="
        + tracks
        + '}';
  }
}
