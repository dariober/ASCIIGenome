package session;

import com.esotericsoftware.yamlbeans.YamlException;
import com.esotericsoftware.yamlbeans.YamlReader;
import exceptions.InvalidGenomicCoordsException;
import exceptions.InvalidRecordException;
import exceptions.InvalidTrackTypeException;
import samTextViewer.GenomicCoords;
import samTextViewer.Utils;
import tracks.Track;
import tracks.TrackPileup;
import tracks.TrackReads;
import tracks.TrackSet;

import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.sql.SQLException;
import java.util.*;
import exceptions.SessionException;

public class Session {
  private ArrayList<Map<String, Object>> sessionList;
  private Map<String, Object> session;

  public Session(String sessionYamlFile, String session) throws FileNotFoundException, YamlException, SessionException {
    // Read yaml file
    // Get genome file if set and parse it into genomicCoords
    // For each track populate trackSet
    YamlReader reader = new YamlReader(new FileReader(sessionYamlFile));
    this.sessionList = (ArrayList<Map<String, Object>>) reader.read();
    this.sortListSessionsLastRead(this.sessionList);
    for (Map<String, Object> ss : this.sessionList) {
      if (((String) ss.get("name")).equals(session)) {
        this.session = ss;
        break;
      }
    }
    if (this.session == null) {
      throw new SessionException("Session '" + session + "' not found in file " + sessionYamlFile);
    }
  }

  public GenomicCoords getGenome() throws IOException, InvalidGenomicCoordsException {
    Map<String, Object> genome = (Map<String, Object>) this.session.get("genome");
    return new GenomicCoords((String) genome.get("region"), Utils.getTerminalWidth(),null, (String) genome.get("fasta"));
  }

  private String getRegion() {
    Map<String, Object> genome = (Map<String, Object>) this.session.get("genome");
    return (String) genome.get("region");
  }

  private void sortListSessionsLastRead(List<Map<String, Object>> ss) {
    ss.sort((o1, o2) -> ((String) o2.get("lastRead")).compareTo((String) o1.get("lastRead")));
  }

  public TrackSet getTrackSet() throws SQLException, InvalidGenomicCoordsException, IOException, InvalidRecordException, ClassNotFoundException, InvalidTrackTypeException {
    List<Map<String, Object>> tracks = (List<Map<String, Object>>) this.session.get("tracks");
    List<String> tl = new ArrayList<String>();
    TrackSet trackSet = new TrackSet(tl, null);
    GenomicCoords gc = new GenomicCoords(this.getRegion(), Utils.getTerminalWidth(), null, null);
    for (Map<String, Object> map : tracks) {
      Track tr;
      String type = (String) map.get("type");
      if (type.equals("TrackPileup")) {
        tr = new TrackPileup((String) map.get("source"), gc);
      } else if (type.equals("TrackReads")) {
        tr = new TrackReads((String) map.get("source"), gc);
      } else {
        throw new InvalidTrackTypeException("Invalid type: " + type);
      }
      trackSet.addTrack(tr, (String) map.get("name"));
    }
    return trackSet;
  }
}
