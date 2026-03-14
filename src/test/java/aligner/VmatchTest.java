package aligner;

import com.google.common.base.Joiner;
import htsjdk.samtools.SAMRecord;
import org.apache.commons.io.FileUtils;
import org.biojava.nbio.alignment.Alignments.PairwiseSequenceAlignerType;
import org.biojava.nbio.core.alignment.template.Profile;
import org.biojava.nbio.core.alignment.template.SequencePair;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.compound.NucleotideCompound;
import org.junit.Test;
import samTextViewer.Utils;

import java.io.File;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;

public class VmatchTest {

  @Test
  public void testAlign() throws Exception {
    Path tempWorkDir = Utils.createTempDir("tmp.", false);

//    this.referenceFasta = String.valueOf(
//            Paths.get(String.valueOf(this.getTempWorkDir()), "reference.fa"));
//    FileUtils.writeStringToFile(new File(referenceFasta), ">reference\n" + reference);


    Vmatch vm = new Vmatch(Path.of("test_data/vmatch/chr7.fa"), tempWorkDir);
    vm.setVmatchExecDir(Path.of("/Users/dario.beraldi/miniforge3/envs/tritume/bin"));
    vm.mkvtree();

    vm.vmatch(Path.of("test_data/vmatch/query.fa"), 4, 3, 5);

    System.out.println(vm.getVmatchAlignmentRecords().stream().map(x -> x.getVmatchAlignmentRecordString()).toList());
    System.out.println(vm.getVmatchAlignmentRecords().stream().map(x -> x.getSAMRecord().getSAMString()).toList());
  }
}
