package aligner;

import com.google.common.base.Joiner;
import com.google.common.base.Splitter;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.SequenceUtil;

import java.util.List;

public class VmatchAlignmentRecord {
  private String vmatchAlignmentRecordString;
  private String query = "";
  private String reference = "";
  private SAMRecord samRecord;

  public VmatchAlignmentRecord(List<String> rawLines, String fullQuerySequence) {
    if (rawLines.isEmpty()) {
      return;
    }
    this.vmatchAlignmentRecordString = Joiner.on("\n").join(rawLines);
    String resultLine = rawLines.remove(0);

    for (String line : rawLines) {
      line = line.trim();
      if (line.startsWith("Sbjct:")) {
        this.reference += line.replaceAll("Sbjct: ", "").replaceAll(" .*", "");
      } else if (line.startsWith("Query:")) {
        this.query += line.replaceAll("Query: ", "").replaceAll(" .*", "");
        }
    }
    if (this.reference.isEmpty() || this.query.isEmpty() || this.query.length() != this.reference.length()) {
      throw new RuntimeException("Invalid alignment");
    }
    this.samRecord = this.makeSAMRecord(resultLine, this.reference, this.query, fullQuerySequence);
  }

  private SAMRecord makeSAMRecord(String result, String reference, String query, String fullQuerySequence) {
    List<String> lst = Splitter.on(" ").omitEmptyStrings().splitToList(result);
    SAMRecord samRecord = new SAMRecord(null);
    samRecord.setReadName(lst.get(5));
    samRecord.setReadBases(fullQuerySequence.getBytes());
    samRecord.setReferenceName(lst.get(1));
    samRecord.setAlignmentStart(Integer.parseInt(lst.get(2)) + 1);
    samRecord.setAttribute("NM", Integer.parseInt(lst.get(7)));
    samRecord.setAttribute("ev", Float.parseFloat(lst.get(8)));
    samRecord.setAttribute("AS", Integer.parseInt(lst.get(9)));
    samRecord.setAttribute("pi", Float.parseFloat(lst.get(10)));
    // samRecord.setAttribute("al", Integer.parseInt(lst.get(0)));

    if (lst.get(3).equals("D")) {
      samRecord.setReadNegativeStrandFlag(false);
    } else if (lst.get(3).equals("P")) {
      samRecord.setReadNegativeStrandFlag(true);
    } else {
      throw new RuntimeException("Unable to determine strand");
    }
    String cigar = "";
    int leftClip = Integer.parseInt(lst.get(6));
    if ( leftClip != 0) {
      cigar += leftClip + "S";
    }
    cigar += this.generateCigar(query, reference);
    int rightClip = fullQuerySequence.length() - query.replaceAll("-", "").length() - leftClip;
    if (rightClip != 0 ) {
      cigar += rightClip + "S";
    }
    samRecord.setCigarString(cigar);
    return samRecord;
  }

  private String generateCigar(String query, String subject) {
    if (query.length() != subject.length()) {
      throw new IllegalArgumentException("Aligned strings must have the same length");
    }

    StringBuilder cigar = new StringBuilder();
    char prevOp = ' ';
    int count = 0;

    for (int i = 0; i < query.length(); i++) {
      char q = query.charAt(i);
      char s = subject.charAt(i);

      char op;

      if (q == '-') {
        op = 'D'; // deletion from query
      } else if (s == '-') {
        op = 'I'; // insertion in query
      } else {
        op = 'M'; // match or mismatch
      }

      if (i == 0) {
        prevOp = op;
        count = 1;
      } else if (op == prevOp) {
        count++;
      } else {
        cigar.append(count).append(prevOp);
        prevOp = op;
        count = 1;
      }
    }

    cigar.append(count).append(prevOp);
    return cigar.toString();
  }

  public SAMRecord getSAMRecord() {
    return this.samRecord;
  }

  public String getVmatchAlignmentRecordString() {
    return this.vmatchAlignmentRecordString;
  }

}
