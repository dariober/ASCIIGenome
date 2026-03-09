package aligner;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.SAMRecord;
import java.util.Arrays;
import org.biojava.nbio.alignment.Alignments;
import org.biojava.nbio.alignment.Alignments.PairwiseSequenceAlignerType;
import org.biojava.nbio.alignment.Alignments.PairwiseSequenceScorerType;
import org.biojava.nbio.alignment.SimpleGapPenalty;
import org.biojava.nbio.alignment.template.PairwiseSequenceAligner;
import org.biojava.nbio.core.alignment.matrices.SubstitutionMatrixHelper;
import org.biojava.nbio.core.alignment.template.AlignedSequence;
import org.biojava.nbio.core.alignment.template.Profile;
import org.biojava.nbio.core.alignment.template.SequencePair;
import org.biojava.nbio.core.alignment.template.SubstitutionMatrix;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.compound.NucleotideCompound;

public class Aligner {

  private final SubstitutionMatrix<NucleotideCompound> matrix;
  private final SimpleGapPenalty gapPenalty;
  private SequencePair<DNASequence, NucleotideCompound> pair;
  private byte[] readBases;

  public Aligner() {
    // Simple nucleotide scoring matrix
    this.matrix = SubstitutionMatrixHelper.getNuc4_4();

    // Gap penalties
    gapPenalty = new SimpleGapPenalty();
    gapPenalty.setOpenPenalty((short) 5);
    gapPenalty.setExtensionPenalty((short) 2);
  }

  public void align (
      String queryDNA,
      String referenceDNA,
      PairwiseSequenceAlignerType pairwiseSequenceAlignerType
  ) throws Exception {

    DNASequence query = new DNASequence(queryDNA);
    DNASequence reference = new DNASequence(referenceDNA);
    this.readBases = queryDNA.getBytes();

    PairwiseSequenceAligner<DNASequence, NucleotideCompound> aligner =
        Alignments.getPairwiseAligner(
            query,
            reference,
            pairwiseSequenceAlignerType, // Smith-Waterman
            gapPenalty,
            matrix
        );
    this.pair = aligner.getPair();

    double xx = Alignments.getAllPairsScores(
        Arrays.asList(query, reference), PairwiseSequenceScorerType.LOCAL, new SimpleGapPenalty(),
        matrix)[0];
  }

  private String toCigar() {

    AlignedSequence<DNASequence, NucleotideCompound> q = this.getPair().getQuery();

    int alnLen = this.getPair().getLength();
    int queryLen = q.getOriginalSequence().getLength();

    int firstQueryIdx = -1;
    int lastQueryIdx = -1;

    for (int i = 1; i <= alnLen; i++) {
      int seqIndex = q.getSequenceIndexAt(i);

      if (seqIndex != -1) {
        if (firstQueryIdx == -1) firstQueryIdx = seqIndex;
        lastQueryIdx = seqIndex;
      }
    }

    int leftSoft = firstQueryIdx - 1;
    int rightSoft = queryLen - lastQueryIdx;

    StringBuilder cigar = new StringBuilder();

    if (leftSoft > 0) cigar.append(leftSoft).append('S');

    char prev = 0;
    int count = 0;

    for (int i = 1; i <= alnLen; i++) {

      NucleotideCompound qc = this.getPair().getCompoundAt(1, i);
      NucleotideCompound rc = this.getPair().getCompoundAt(2, i);

      char op;

      if (qc == null || qc.toString().equals("-")) op = 'D';
      else if (rc == null || rc.toString().equals("-")) op = 'I';
      else op = 'M';

      if (prev != op && count > 0) {
        cigar.append(count).append(prev);
        count = 0;
      }

      prev = op;
      count++;
    }

    if (count > 0) cigar.append(count).append(prev);

    if (rightSoft > 0) cigar.append(rightSoft).append('S');

    return cigar.toString();
  }

  public SequencePair<DNASequence, NucleotideCompound> getPair() {
    return this.pair;
  }

  private int computeNM() {
    int nm = 0;
    int length = this.getPair().getLength();

    for (int i = 1; i <= length; i++) { // BioJava is 1-based
      NucleotideCompound q = this.getPair().getCompoundAt(1, i);
      NucleotideCompound r = this.getPair().getCompoundAt(2, i);

      if (q == null || q.toString().equals("-")) {
        // Deletion in query
        nm += 1;
      } else if (r == null || r.toString().equals("-")) {
        // Insertion in query
        nm += 1;
      } else if (!q.equals(r)) {
        // Mismatch
        nm += 1;
      }
      // else match → 0
    }

    return nm;
  }

  public SAMRecord toSAMRecord() {
    SAMRecord samRecord = new SAMRecord(null);
    samRecord.setReadBases(this.readBases);
    samRecord.setAlignmentStart(this.getPair().getIndexInTargetForQueryAt(1));
    samRecord.setCigarString(this.toCigar());
    samRecord.setAttribute("NM", this.computeNM());
    return samRecord;
  }
}