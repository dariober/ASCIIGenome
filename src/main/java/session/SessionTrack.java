package session;

import exceptions.InvalidGenomicCoordsException;
import exceptions.InvalidRecordException;
import exceptions.InvalidTrackTypeException;
import java.io.IOException;
import java.lang.reflect.Constructor;
import java.lang.reflect.InvocationTargetException;
import java.sql.SQLException;
import samTextViewer.GenomicCoords;
import tracks.Track;

public class SessionTrack {
  public String filename;
  public String type;
  public String colour;
  public String awk;

  public SessionTrack() {}

  /* Handle Track to make it suitable for serialization */
  public SessionTrack(Track tr) {
    this.filename = tr.getFilename();
    this.awk = tr.getAwk();
    this.colour = tr.getTitleColour();
    this.type = tr.getClass().getName();
  }

  public Track toTrack(String trackName, GenomicCoords gc) throws InvalidTrackTypeException {
    try {
      Constructor<?> c = Class.forName(this.type).getConstructor(String.class, GenomicCoords.class);
      Track tr = (Track) c.newInstance(this.filename, gc);
      tr.setTrackTag(trackName);
      tr.setTitleColour(this.colour);
      if (this.awk != null) tr.setAwk(this.awk);
      return tr;
    } catch (NoSuchMethodException
        | InstantiationException
        | IllegalAccessException
        | InvocationTargetException
        | ClassNotFoundException
        | IOException
        | InvalidGenomicCoordsException
        | InvalidRecordException
        | SQLException e) {
      throw new InvalidTrackTypeException(
          "Unable to make track '"
              + trackName
              + "' of type '"
              + this.type
              + "'. Got error message:\n"
              + e.toString());
    }
  }

  @Override
  public String toString() {
    return "SessionTrack{"
        + "awk='"
        + awk
        + '\''
        + ", filename='"
        + filename
        + '\''
        + ", type='"
        + type
        + '\''
        + ", colour='"
        + colour
        + '\''
        + '}';
  }
}
