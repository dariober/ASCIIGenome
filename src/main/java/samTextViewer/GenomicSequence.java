package samTextViewer;

import colouring.Config;
import colouring.ConfigKey;
import exceptions.InvalidColourException;
import exceptions.InvalidCommandLineException;
import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.compound.AmbiguityDNACompoundSet;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
import org.biojava.nbio.core.sequence.io.IUPACParser;
import org.biojava.nbio.core.sequence.template.Sequence;
import org.biojava.nbio.core.sequence.transcription.Frame;
import org.biojava.nbio.core.sequence.transcription.TranscriptionEngine;

import java.util.*;

public class GenomicSequence {
    private final byte[] sequence;
    private boolean noFormat = false;
    private String geneticCode = "UNIVERSAL";
    private Map<Frame, Sequence<AminoAcidCompound>> sixFrameTranslation = new HashMap<>();
    private ArrayList<Frame> frames = new ArrayList<>();

    public GenomicSequence(byte[] sequence) throws CompoundNotFoundException {
        this.sequence = sequence;
        IUPACParser.IUPACTable table = IUPACParser.getInstance().getTable(this.geneticCode);
        TranscriptionEngine transcriptionEngine = new TranscriptionEngine.Builder()
                .table(table)
                .initMet(true)
                .trimStop(false)
                .build();
        DNASequence dna = new DNASequence(new String(this.sequence), AmbiguityDNACompoundSet.getDNACompoundSet());
        for (Frame frame : Frame.getAllFrames()) {
            Map<Frame, Sequence<AminoAcidCompound>> tr = new LinkedHashMap<>();
            try {
                tr = transcriptionEngine.multipleFrameTranslation(dna, frame);
            } catch(Exception e) {
                //
            }
            this.sixFrameTranslation.put(frame, tr.get(frame));
        }
    }

    public String getPrintableSequence() throws InvalidColourException, CompoundNotFoundException, InvalidCommandLineException {
        if (this.sequence == null) {
            return "";
        }
        StringBuilder faSeqStr = new StringBuilder();

        for (Frame x : Frame.getForwardFrames()) {
            if (this.frames.contains(x)) {
                faSeqStr.append(this.proteinToString(x)).append('\n');
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
                    faSeqStr.append(prefix).append(Config.get256Colour(ConfigKey.seq_a)).append("m").append(base);
                } else if (base == 'C' || base == 'c') {
                    faSeqStr.append(prefix).append(Config.get256Colour(ConfigKey.seq_c)).append("m").append(base);
                } else if (base == 'G' || base == 'g') {
                    faSeqStr.append(prefix).append(Config.get256Colour(ConfigKey.seq_g)).append("m").append(base);
                } else if (base == 'T' || base == 't') {
                    faSeqStr.append(prefix).append(Config.get256Colour(ConfigKey.seq_t)).append("m").append(base);
                } else {
                    faSeqStr.append(prefix).append(Config.get256Colour(ConfigKey.seq_other)).append("m").append(base);
                }
            }
        }
        faSeqStr.append('\n');

        for (Frame x : Frame.getReverseFrames()) {
            if (this.frames.contains(x)) {
                faSeqStr.append(this.proteinToString(x)).append('\n');
            }
        }
        return faSeqStr.toString();
    }

    protected LinkedHashMap<Integer, String> geneticCodeNames() {
        LinkedHashMap<Integer, String> tables = new LinkedHashMap<Integer, String>();
        IUPACParser.getInstance().getTables().forEach((x) -> {
            tables.put(x.getId(), x.getName());
        });
        return tables;
    }

    private String proteinToString(Frame frame) throws InvalidColourException {
        Sequence<AminoAcidCompound> protein = this.sixFrameTranslation.get(frame);
        StringBuilder out = new StringBuilder();
        int width = 0;
        if (frame.equals(Frame.TWO)) {
            out.append("_");
            width += 1;
        }
        if (frame.equals(Frame.THREE)) {
            out.append("__");
            width += 2;
        }

        if (protein != null) {
            for (AminoAcidCompound aa : protein.getAsList()) {
                String aaStr = "_" + aa.getShortName() + "_";
                if (!this.noFormat && aa.getShortName().equals("*")) {
                    String prefix = "\033[48;5;" + Config.get256Colour(ConfigKey.background) + ";38;5;";
                    aaStr = prefix + Config.get256Colour(ConfigKey.stop_codon) + "m" + aaStr + "\033[0m";
                    System.out.println("HERE " + aaStr);
                }
                out.append(aaStr);
                width += 3;
            }
        }

        if (frame.equals(Frame.REVERSED_TWO)) {
            out.insert(0, "_");
            width += 1;
        }
        if (frame.equals(Frame.REVERSED_THREE)) {
            out.insert(0,"__");
            width += 2;
        }

        if (frame.name().startsWith("REVERSED")) {
            String str = new String(out.reverse()).toLowerCase();
            int npad = this.sequence.length - width;
            if (npad > 0) {
                str = "_".repeat(npad) + str;
            }
            return str;
        }
        return new String(out).toUpperCase();
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

    public void setGeneticCode(String geneticCode) {
        this.geneticCode = geneticCode;
    }

    public ArrayList<Frame> getFrames() {
        return frames;
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
}
