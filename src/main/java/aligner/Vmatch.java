package aligner;

import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Stream;

import com.google.common.base.Joiner;
import com.google.common.base.Splitter;
import faidx.Faidx;
import faidx.UnindexableFastaFileException;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import org.apache.commons.io.FileUtils;
import samTextViewer.Utils;
import systemCommand.SystemCommand;

public class Vmatch {
  // private String query;
  private final Path referenceFasta;
  private final Path referenceIndex;
  Path vmatchExecDir = Paths.get("");
  private final Path workDir;
  private String header;
  private List<VmatchAlignmentRecord> vmatchAlignmentRecords = new ArrayList<>();

  // CONDA_SUBDIR=osx-64 conda install vmatch
  public Vmatch(Path referenceFasta, Path workDir) throws IOException {
    this.referenceFasta = referenceFasta;
    this.workDir = workDir;
    this.referenceIndex = Paths.get(String.valueOf(this.workDir), String.valueOf(referenceFasta.getFileName()));
  }

  public void mkvtree() throws IOException, InterruptedException {


    String mkvtree = String.valueOf(Paths.get(String.valueOf(this.getVmatchExecDir()), "mkvtree"));
    List<String> cmd = new ArrayList<>();
    cmd.add(mkvtree);
    cmd.add("-dna");
    cmd.add("-pl");
    cmd.add("-allout");
    cmd.add("-db");
    cmd.add(this.referenceFasta.toString());
    cmd.add("-indexname");
    cmd.add(String.valueOf(this.referenceIndex));
    Utils.execSystemCommand(new String[]{}, cmd);
  }

  public void vmatch(Path queryFasta, int minMatchLen, int editDistance, int best) throws IOException, InterruptedException {
    String vmatch = String.valueOf(Paths.get(String.valueOf(this.getVmatchExecDir()), "vmatch"));
    String opts = " -sort ea -d -p -showdesc 0 -s -e " + editDistance + " -l " + minMatchLen;
    opts += best > 0 ? " -best " + best : "";
    List<String> cmd = Splitter.on(" ").splitToList(vmatch + opts + " -q " + queryFasta + " " + this.referenceIndex);
    List<String> out = Utils.execSystemCommand(new String[]{}, cmd);
    this.vmatchAlignmentRecords = this.processVmatchOutput(out.stream(), this.getQuerySequence(queryFasta.toFile()));
  }

  private Map<String,  String> getQuerySequence(File queryFasta) {
    ReferenceSequenceFile ref =
            ReferenceSequenceFileFactory.getReferenceSequenceFile(queryFasta);
    ReferenceSequence seq;
    Map<String, String> querySeq = new HashMap<>();
    while ((seq = ref.nextSequence()) != null) {
      querySeq.put(seq.getName(), seq.getBaseString());
    }
    return querySeq;
  }

  private List<VmatchAlignmentRecord> processVmatchOutput(Stream<String> vmatchOutput, Map<String, String> querySequence) {
    List<String> rawLines = new ArrayList<>();
    List<VmatchAlignmentRecord> vmatchAlignmentRecord = new ArrayList<>();
    Iterator<String> iter = vmatchOutput.iterator();
    while (iter.hasNext()) {
      String line = iter.next();
      if (line.startsWith("#")) {
        this.header = line;
        continue;
      }
      if (!rawLines.isEmpty() && rawLines.get(rawLines.size() - 1).isEmpty() && line.isEmpty()) {
        rawLines.remove(rawLines.size() - 1);
        String queryName = Splitter.on(" ").omitEmptyStrings().splitToList(rawLines.get(0)).get(5);
        String fullQuerySequence = querySequence.get(queryName);
        vmatchAlignmentRecord.add(new VmatchAlignmentRecord(rawLines, fullQuerySequence));
        rawLines.clear();
      } else {
        rawLines.add(line);
      }
    }
    return vmatchAlignmentRecord;
  }

  public List<VmatchAlignmentRecord> getVmatchAlignmentRecords(){
    return this.vmatchAlignmentRecords;
  }

  public Path getWorkDir() {
    return this.workDir;
  }

  public void setVmatchExecDir(Path vmatchExecDir) {
    this.vmatchExecDir = vmatchExecDir;
  }

  public Path getVmatchExecDir() {
    return vmatchExecDir;
  }
}
