package samTextViewer;

import colouring.Config;
import colouring.ConfigKey;
import exceptions.InvalidColourException;
import exceptions.InvalidGenomicCoordsException;
import java.util.*;
import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.exceptions.TranslationException;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
import org.biojava.nbio.core.sequence.io.IUPACParser;
import org.biojava.nbio.core.sequence.io.IUPACParser.IUPACTable;
import org.biojava.nbio.core.sequence.template.Sequence;
import org.biojava.nbio.core.sequence.transcription.Frame;
import org.biojava.nbio.core.sequence.transcription.TranscriptionEngine;
import org.biojava.nbio.core.sequence.transcription.TranscriptionEngine.Builder;
import tracks.FeatureChar;

public class GenomicSequence {

  private final byte[] sequence;
  private boolean noFormat = false;
  private String geneticCode = "UNIVERSAL";
  private final Map<Frame, Sequence<AminoAcidCompound>> sixFrameTranslation = new HashMap<>();
  private List<Frame> frames = new ArrayList<>();
  private PrintCodon printCodon = PrintCodon.ALL;
  private final Integer forwardOffset;
  private final Integer reverseOffset;

  protected GenomicSequence(byte[] sequence, Integer genomicPositionStart, Integer chromSize)
      throws InvalidGenomicCoordsException {
    this.sequence = sequence;
    if (sequence == null) {
      this.forwardOffset = null;
      this.reverseOffset = null;
    } else {
      if (chromSize - (genomicPositionStart + sequence.length - 1) < 0) {
        throw new RuntimeException("Invalid sequence length or start position or chromosome size");
      }
      this.forwardOffset = (genomicPositionStart - 1) % 3;
      this.reverseOffset = (chromSize - (genomicPositionStart - 1) - sequence.length) % 3;
      this.translate();
    }
  }

  protected GenomicSequence() throws InvalidGenomicCoordsException {
    this.sequence = null;
    this.forwardOffset = null;
    this.reverseOffset = null;
  }

  private void translate() throws InvalidGenomicCoordsException {
    if (this.sequence == null) {
      return;
    }
    IUPACTable table = IUPACParser.getInstance().getTable(this.geneticCode.toUpperCase());
    if (table == null) {
      throw new InvalidGenomicCoordsException(
          "Invalid translation table " + this.geneticCode.toUpperCase());
    }
    TranscriptionEngine transcriptionEngine =
        new Builder().table(table).initMet(true).trimStop(false).build();

    HashMap<Frame, Sequence<AminoAcidCompound>> sixFrame = new HashMap<>();
    DNASequence dna;
    try {
      dna = new DNASequence(new String(this.sequence));
    } catch (CompoundNotFoundException e) {
      throw new InvalidGenomicCoordsException(e.getMessage());
    }
    for (Frame frame : Frame.getAllFrames()) {
      Map<Frame, Sequence<AminoAcidCompound>> tr = Map.of();
      try {
        tr = transcriptionEngine.multipleFrameTranslation(dna, frame);
      } catch (TranslationException e) {
        // System.err.println(frame);
      }
      sixFrame.put(frame, tr.get(frame));
    }
    this.sixFrameTranslation.putAll(this.recodeFrames(sixFrame));
  }

  private Map<Frame, Sequence<AminoAcidCompound>> recodeFrames(
      Map<Frame, Sequence<AminoAcidCompound>> tr) {
    Map<Frame, Sequence<AminoAcidCompound>> tr2 = new HashMap<>(Map.of());
    tr.forEach(
        (frame, sequence) -> {
          // MEMO: frame.ordinal() is between 0 and 5
          int oldFrame = frame.ordinal() <= 2 ? frame.ordinal() + 1 : frame.ordinal() - 2;

          if (frame.ordinal() <= 2) {
            int newForwardFrame = (oldFrame + this.forwardOffset) % 3;
            newForwardFrame = newForwardFrame == 0 ? 3 : newForwardFrame;
            if (newForwardFrame == 1) {
              tr2.put(Frame.ONE, sequence);
            } else if (newForwardFrame == 2) {
              tr2.put(Frame.TWO, sequence);
            } else if (newForwardFrame == 3) {
              tr2.put(Frame.THREE, sequence);
            } else {
              throw new RuntimeException("Unexpected frame");
            }
          } else {
            int newReverseFrame = (oldFrame + this.reverseOffset) % 3;
            newReverseFrame = newReverseFrame == 0 ? 3 : newReverseFrame;
            if (newReverseFrame == 1) {
              tr2.put(Frame.REVERSED_ONE, sequence);
            }
            if (newReverseFrame == 2) {
              tr2.put(Frame.REVERSED_TWO, sequence);
            }
            if (newReverseFrame == 3) {
              tr2.put(Frame.REVERSED_THREE, sequence);
            }
          }
        });
    return tr2;
  }

  public String getPrintableSequence()
      throws InvalidColourException, InvalidGenomicCoordsException {
    if (this.sequence == null) {
      return "";
    }
    StringBuilder faSeqStr = new StringBuilder();
    for (int i = Frame.getForwardFrames().length - 1; i >= 0; i--) {
      Frame x = Frame.getForwardFrames()[i];
      if (this.frames.contains(x)) {
        faSeqStr.append(this.proteinToString(x, this.printCodon)).append('\n');
      }
    }

    if (this.noFormat) {
      faSeqStr.append(new String(this.sequence));
    } else {
      for (byte c : this.sequence) {
        // For colour scheme see http://www.umass.edu/molvis/tutorials/dna/atgc.htm
        char base = (char) c;
        String prefix = "\033[48;5;" + Config.get256Colour(ConfigKey.background) + ";38;5;";
        if (base == 'A' || base == 'a') {
          faSeqStr
              .append(prefix)
              .append(Config.get256Colour(ConfigKey.seq_a))
              .append("m")
              .append(base);
        } else if (base == 'C' || base == 'c') {
          faSeqStr
              .append(prefix)
              .append(Config.get256Colour(ConfigKey.seq_c))
              .append("m")
              .append(base);
        } else if (base == 'G' || base == 'g') {
          faSeqStr
              .append(prefix)
              .append(Config.get256Colour(ConfigKey.seq_g))
              .append("m")
              .append(base);
        } else if (base == 'T' || base == 't') {
          faSeqStr
              .append(prefix)
              .append(Config.get256Colour(ConfigKey.seq_t))
              .append("m")
              .append(base);
        } else {
          faSeqStr
              .append(prefix)
              .append(Config.get256Colour(ConfigKey.seq_other))
              .append("m")
              .append(base);
        }
      }
    }
    faSeqStr.append('\n');

    for (Frame x : Frame.getReverseFrames()) {
      if (this.frames.contains(x)) {
        faSeqStr.append(this.proteinToString(x, this.printCodon)).append('\n');
      }
    }
    return faSeqStr.toString();
  }

  protected LinkedHashMap<Integer, String> geneticCodeNames() {
    LinkedHashMap<Integer, String> tables = new LinkedHashMap<Integer, String>();
    IUPACParser.getInstance()
        .getTables()
        .forEach(
            (x) -> {
              tables.put(x.getId(), x.getName());
            });
    return tables;
  }

  //  private String startStopCodonsToString(Sequence<AminoAcidCompound> protein) {
  //
  //  }
  //
  private String proteinToString(Frame frame, PrintCodon printCodon)
      throws InvalidColourException, InvalidGenomicCoordsException {
    Sequence<AminoAcidCompound> protein = this.sixFrameTranslation.get(frame);
    ArrayList<FeatureChar> fmtSeq = new ArrayList<>();

    int sidePadding =
        frame.name().startsWith("REVERSED")
            ? (3 - this.reverseOffset) % 3
            : (3 - this.forwardOffset) % 3;
    if (frame.equals(Frame.TWO) || frame.equals(Frame.REVERSED_TWO)) {
      sidePadding = (sidePadding + 1) % 3;
    }
    if (frame.equals(Frame.THREE) || frame.equals(Frame.REVERSED_THREE)) {
      sidePadding = (sidePadding + 2) % 3;
    }

    int width = 0;
    for (int i = 0; i < sidePadding; i++) {
      FeatureChar c = new FeatureChar();
      c.setText(' ');
      c.setFgColour(Config.get(ConfigKey.codon));
      fmtSeq.add(c);
      width += 1;
    }

    if (protein != null) {
      for (AminoAcidCompound aa : protein.getAsList()) {
        char[] aaStr = {' ', aa.getShortName().charAt(0), ' '};
        for (char s : aaStr) {
          FeatureChar c = new FeatureChar();
          if (printCodon.equals(PrintCodon.ALL)) {
            c.setText(s);
          } else if ((printCodon.equals(PrintCodon.START)) && s == 'M') {
            c.setText(s);
          } else if (printCodon.equals(PrintCodon.STOP) && s == '*') {
            c.setText(s);
          } else if (printCodon.equals(PrintCodon.START_AND_STOP) && (s == '*' || s == 'M')) {
            c.setText(s);
          } else {
            c.setText(' ');
          }

          if (aa.getShortName().equals("*")
              && ((printCodon.equals(PrintCodon.STOP)
                  || printCodon.equals(PrintCodon.START_AND_STOP)
                  || printCodon.equals(PrintCodon.ALL)))) {
            c.setFgColour(Config.get(ConfigKey.stop_codon));
            c.setInvertFgBgColour(true);
          } else if (aa.getShortName().equalsIgnoreCase("M")
              && ((printCodon.equals(PrintCodon.START)
                  || printCodon.equals(PrintCodon.START_AND_STOP)
                  || printCodon.equals(PrintCodon.ALL)))) {
            c.setFgColour(Config.get(ConfigKey.start_codon));
            c.setInvertFgBgColour(true);
          } else {
            c.setFgColour(Config.get(ConfigKey.codon));
          }
          fmtSeq.add(c);
        }
        width += 3;
      }
    }

    if (frame.name().startsWith("REVERSED")) {
      Collections.reverse(fmtSeq);
      int npad = this.sequence.length - width;
      while (npad > 0) {
        FeatureChar c = new FeatureChar();
        c.setText(' ');
        c.setFgColour(Config.get(ConfigKey.codon));
        fmtSeq.add(0, c);
        npad--;
      }
      fmtSeq.set(fmtSeq.size() - 1, this.formatFrameChar(frame));
    } else {
      fmtSeq.set(0, this.formatFrameChar(frame));
    }
    StringBuilder sb = new StringBuilder();
    for (FeatureChar x : fmtSeq) {
      sb.append(x.format(this.noFormat));
    }
    return new String(sb);
  }

  private FeatureChar formatFrameChar(Frame frame) {
    FeatureChar c = new FeatureChar();
    if (frame.ordinal() == 0 || frame.ordinal() == 1 || frame.ordinal() == 2) {
      c.setText(String.valueOf(frame.ordinal() + 1).charAt(0));
    } else if (frame.ordinal() == 3 || frame.ordinal() == 4 || frame.ordinal() == 5) {
      c.setText(String.valueOf(frame.ordinal() - 2).charAt(0));
    } else {
      throw new RuntimeException();
    }
    c.setInvertFgBgColour(true);
    return c;
  }

  public boolean isNoFormat() {
    return noFormat;
  }

  public void setNoFormat(boolean noFormat) {
    this.noFormat = noFormat;
  }

  public String getGeneticCode() {
    return geneticCode;
  }

  public void setGeneticCode(String geneticCode) throws InvalidGenomicCoordsException {
    this.geneticCode = geneticCode;
    this.translate();
  }

  public ArrayList<Frame> getFrames() {
    return (ArrayList<Frame>) frames;
  }

  public void setFrames(ArrayList<Frame> frames) {
    this.frames = frames;
  }

  public void setFrames(Frame frame) {
    ArrayList<Frame> frames = new ArrayList<>();
    frames.add(frame);
    this.frames = frames;
  }

  public void setFrames(Frame[] frames) {
    this.frames = new ArrayList<>(Arrays.asList(frames));
  }

  public PrintCodon getPrintCodon() {
    return printCodon;
  }

  public void setPrintCodon(PrintCodon printCodon) {
    this.printCodon = printCodon;
  }

  public byte[] getSequence() {
    return sequence;
  }

  public String toString() {
    String sequence = "n/a";
    if (this.getSequence() != null) {
      sequence =
          new String(this.getSequence()).substring(0, Math.min(10, this.getSequence().length));
    }
    String frames = String.valueOf(this.getFrames());
    return "Frames: "
        + frames
        + "\nSequence (start of): "
        + sequence
        + "\nreverseOffset"
        + this.reverseOffset
        + "\nforwardOffset"
        + this.forwardOffset
        + "\nGenetic code: "
        + this.getGeneticCode();
  }
}
