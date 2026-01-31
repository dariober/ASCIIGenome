package tracks;

import colouring.Config;
import colouring.ConfigKey;
import exceptions.InvalidColourException;
import exceptions.InvalidGenomicCoordsException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/** Model a sequenced fragment typically represented by a pair of reads. */
class SamSequenceFragment {

  private static final char CONNECTING_CHAR = '~';
  private TextRead leftRead = null;
  private TextRead rightRead = null;
  private boolean isSingleton;

  /** Singleton if read has properly paired flag not set. */

  // C O N S T R U C T O R S

  protected SamSequenceFragment(TextRead tr) {
    this.leftRead = tr;
    // Only one read is present in this fragment.
    // We need to decide whether it is a singleton. It may not be a singleton
    // if the mate of this read exist but it wasn't found.
    if (tr.getSamRecord().getReadPairedFlag()
        && tr.getSamRecord().getProperPairFlag()
        && tr.getSamRecord().getMateAlignmentStart() > 0) {
      // Note that the position of the mate must be known.
      this.isSingleton = false;
    } else {
      this.isSingleton = true;
    }
  }

  protected SamSequenceFragment(TextRead tr, TextRead mate) {
    if (tr.getAlignmentStart() <= mate.getAlignmentStart()) {
      this.leftRead = tr;
      this.rightRead = mate;
    } else {
      this.leftRead = mate;
      this.rightRead = tr;
    }
    this.isSingleton = false;
  }

  // M E T H O D S

  protected int getTextStart() {
    if (this.getRightRead() == null) {
      return this.getLeftRead().getTextStart();
    } else {
      return Math.min(this.getRightRead().getTextStart(), this.getLeftRead().getTextStart());
    }
  }

  protected int getTextEnd() {
    if (this.getRightRead() == null) {
      return this.getLeftRead().getTextEnd();
    } else {
      return Math.max(this.getRightRead().getTextEnd(), this.getLeftRead().getTextEnd());
    }
  }

  TextRead getLeftRead() {
    return this.leftRead;
  }

  protected TextRead getRightRead() {
    return this.rightRead;
  }

  protected String getPrintableFragment(boolean bs, boolean noFormat)
      throws IOException, InvalidGenomicCoordsException, InvalidColourException {
    if (this.getRightRead() == null || this.isSingleton) {
      return this.getLeftRead().getPrintableTextRead(bs, noFormat, false);
    } else {
      // We need a string that goes from start to end and fills up the bit
      // in the middle or make pairs overlap:
      // AAAAAAAAA-----TTTTTTTT
      // Overlap:
      // AAAAAAAAA
      //        tttttttttt
      // becomes:
      // AAAAAAAtttttttttt
      //
      // Prepare chars:
      List<FeatureChar> lst = new ArrayList<FeatureChar>();
      int len = this.getTextEnd() - this.getTextStart() + 1; // MEMO: Start is 1-based
      for (int i = 0; i < len; i++) {
        FeatureChar x = new FeatureChar();
        x.setText(CONNECTING_CHAR);
        lst.add(x);
      }
      int i = 0;
      for (FeatureChar x : this.getLeftRead().getTextReadAsFeatureChars(bs)) {
        lst.set(i, x);
        i++;
      }
      // Fill up with the right read. NB: Reads might be fully contained one into the other
      // Where does the right read starts from?
      int from = this.getRightRead().getTextStart() - this.getLeftRead().getTextStart();
      for (FeatureChar rightNt : this.getRightRead().getTextReadAsFeatureChars(bs)) {
        FeatureChar leftNt = lst.get(from);
        if (leftNt.getText() == CONNECTING_CHAR) {
          lst.set(from, rightNt);
        } else if (!leftNt.isFaint() && !rightNt.isFaint()) {
          if (Character.toUpperCase(leftNt.getText()) == Character.toUpperCase(rightNt.getText())) {
            lst.set(from, rightNt);
          } else {
            FeatureChar n = (FeatureChar) rightNt.clone();
            n.setText(Character.isUpperCase(n.getText()) ? 'N' : 'n');
            lst.set(from, n);
          }
        } else if (!leftNt.isFaint() && rightNt.isFaint()) {
          lst.set(from, leftNt);
        } else if (leftNt.isFaint() && !rightNt.isFaint()) {
          lst.set(from, rightNt);
        } else if (leftNt.isFaint() && rightNt.isFaint()) {
          // Both faint, just pick one.
          lst.set(from, rightNt);
        }
        from++;
      }
      StringBuilder sb = new StringBuilder();
      for (FeatureChar x : lst) {
        sb.append(x.format(noFormat));
      }
      return sb.toString();
    }
  }

  protected void setSingleton(boolean singleton) {
    this.isSingleton = singleton;
  }
}
