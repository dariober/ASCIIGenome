package utils;

import java.util.HashMap;
import java.util.Map;
import tracks.TrackFormat;

public class CsvPresets {
  public static Map<TrackFormat, CsvFormat> map = new HashMap<>();

  static {
    map.put(TrackFormat.BED, new CsvFormat(0, 1, 2, -1, true, -1, '#', '\t'));
    map.put(TrackFormat.BEDGRAPH, new CsvFormat(0, 1, 2, 3, true, -1, '#', '\t'));
    map.put(TrackFormat.GFF, new CsvFormat(0, 3, 4, -1, false, -1, '#', '\t'));
    map.put(TrackFormat.GTF, new CsvFormat(0, 3, 4, -1, false, -1, '#', '\t'));
  }

  public static CsvFormat get(TrackFormat format) {
    return map.get(format);
  }
}
