package aligner;

import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import org.apache.commons.io.FileUtils;
import samTextViewer.Utils;

public class Vmatch {
  // private String query;
  // private String reference;
  Path tempWorkDir;

  public Vmatch() throws IOException {
    this.tempWorkDir = Utils.createTempDir(".asciigenome_vmatch", true);
  }

  public void mkvtree(String reference) throws IOException {
    ProcessBuilder builder = new ProcessBuilder();

    String referenceFile = String.valueOf(
        Paths.get(String.valueOf(this.getTempWorkDir()), "reference.fa"));
    FileUtils.writeStringToFile(new File(referenceFile), ">reference\n" + reference);

    String cmd = "mkvtree -pl -tis -ois -dna -db " + referenceFile;
    builder.command("bash", "-c", cmd);
  }

  public Path getTempWorkDir() {
    return this.tempWorkDir;
  }


}
