package tracks;

import aligner.Sassy;
import colouring.Config;
import colouring.ConfigKey;
import exceptions.InvalidGenomicCoordsException;
import exceptions.InvalidRecordException;
import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.List;
import samTextViewer.GenomicCoords;
import samTextViewer.Utils;

public class TrackSassySearch extends TrackReads {

  private List<String> sassyCliOpts = new ArrayList<>();
  private final Path tempWorkDir;

  public TrackSassySearch(List<String> sassyCliOpts, GenomicCoords gc)
      throws IOException,
          InvalidGenomicCoordsException {
    this.downsampleUsingTemplateName = false;
    this.MAX_REGION_SIZE = Integer.MAX_VALUE;
    this.sassyCliOpts = sassyCliOpts;
    this.tempWorkDir = Utils.createTempDir(".asciigenome.sassy.", true);
    String bamFile = Paths.get(this.tempWorkDir.toString(), "aln.bam").toString();
    this.setWorkFilename(bamFile);
    this.setFilename(bamFile);
    this.setTrackFormat(TrackFormat.BAM);
    this.gc = gc;
    this.update();
  }

  @Override
  public void update() throws InvalidGenomicCoordsException, IOException {
    this.sassySearch();
    super.update();
  }

  private void sassySearch() throws IOException {
    Sassy sassy = new Sassy(this.getGc(), this.tempWorkDir);
    sassy.setExecPath(Path.of(Config.get(ConfigKey.sassy)));
    try {
      sassy.search(this.getSassyCliOpts());
    } catch (InterruptedException e) {
      throw new RuntimeException(e.getMessage());
    }
    String samFile = this.getWorkFilename().replaceAll(".bam$", ".sam");
    sassy.writeSAMFile(Path.of(samFile));
    Utils.sortAndIndexSamOrBam(samFile, this.getWorkFilename(), true, this.getGc().getFastaFile());
    Path.of(samFile).toFile().delete();
  }

  protected List<String> getSassyCliOpts() {
    return this.sassyCliOpts;
  }

  protected void setSassyCliOpts(List<String> sassyCliOpts)
      throws InvalidGenomicCoordsException, IOException {
    this.sassyCliOpts = sassyCliOpts;
    this.sassySearch();
    this.update();
  }
}
