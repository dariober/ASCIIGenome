package tracks;

import org.apache.commons.lang3.NotImplementedException;

/** Info about a screen position for wiggle-like data. Somewhat counterpart to ScreenLocusInfo */
public class ScreenWiggleLocusInfo {

  private int cntGenomicLoci = 0; // Count of genomic loci mapped to this screen position
  private double sumScore =
      0; // Sum of scores accumulated from wiggle sites mapped to this screen locus
  private double min = Double.MAX_VALUE;
  private double max = Double.MIN_VALUE;

  /* C o n s t r u c t o r */
  public ScreenWiggleLocusInfo() {}

  /* M e t h o d s */
  /** Increment attributes by given score */
  public void increment(double score, DataTransformation dataTransformation) {
    //    if (score <= 0 && (dataTransformation == DataTransformation.LOG10 || dataTransformation ==
    // DataTransformation.MINUS_LOG10)) {
    //      throw  new RuntimeException("Invalid transformation for value: " + score);
    //    }
    if (dataTransformation == DataTransformation.IDENTITY) {
      //
    } else if (dataTransformation == DataTransformation.LOG10) {
      score = Math.log10(score);
    } else if (dataTransformation == DataTransformation.MINUS_LOG10) {
      score = -Math.log10(score);
    } else {
      throw new NotImplementedException(
          "Transformation " + dataTransformation + " not implemented yet");
    }
    cntGenomicLoci++;
    sumScore += score;
    this.min = Math.min(this.min, score);
    this.max = Math.max(this.max, score);
  }

  public String toString() {
    String str = "cntGenomicLoci: " + this.cntGenomicLoci + "; sumScores: " + this.sumScore;
    return str;
  }

  /*   G e t t e r s   */
  protected double getMeanScore() {
    return this.sumScore / this.cntGenomicLoci;
  }

  public double getMax() {
    if (this.cntGenomicLoci > 0) {
      return max;
    } else {
      return Double.NaN;
    }
  }

  public double getMin() {
    if (this.cntGenomicLoci > 0) {
      return min;
    } else {
      return Double.NaN;
    }
  }
}
