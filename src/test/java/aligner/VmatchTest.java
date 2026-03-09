package aligner;

import htsjdk.samtools.SAMRecord;
import org.biojava.nbio.alignment.Alignments.PairwiseSequenceAlignerType;
import org.biojava.nbio.core.alignment.template.Profile;
import org.biojava.nbio.core.alignment.template.SequencePair;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.compound.NucleotideCompound;
import org.junit.Test;

public class VmatchTest {

  @Test
  public void testAlign() throws Exception {
    Vmatch vm = new Vmatch();
    System.out.println(vm.getTempWorkDir());
  }
}
