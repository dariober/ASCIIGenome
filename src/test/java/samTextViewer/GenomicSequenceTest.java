package samTextViewer;

import colouring.Config;
import exceptions.InvalidColourException;
import exceptions.InvalidCommandLineException;
import exceptions.InvalidConfigException;
import exceptions.InvalidGenomicCoordsException;

import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.RNASequence;
import org.biojava.nbio.core.sequence.compound.AmbiguityDNACompoundSet;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
import org.biojava.nbio.core.sequence.io.IUPACParser;
import org.biojava.nbio.core.sequence.template.Sequence;
import org.biojava.nbio.core.sequence.transcription.Frame;
import org.biojava.nbio.core.sequence.transcription.TranscriptionEngine;
import org.biojava.nbio.core.sequence.transcription.Table;

import org.junit.Before;
import org.junit.Test;

import static org.junit.Assert.assertEquals;

public class GenomicSequenceTest {

    @Before
    public void initConfig() throws IOException, InvalidConfigException {
        new Config(null);
    }

    @Test
    public void canGetGeneticCodeNames()
            throws CompoundNotFoundException, InvalidGenomicCoordsException {
        String dna = "";
        GenomicSequence gs = new GenomicSequence(dna.getBytes());
        assertEquals(17, gs.geneticCodeNames().size());
        assertEquals("UNIVERSAL", gs.geneticCodeNames().get(1));
    }

    @Test
    public void canTranslateSequence()
            throws InvalidColourException, CompoundNotFoundException, InvalidCommandLineException, InvalidGenomicCoordsException {
        String dna = "ATGCTGTAG";
        GenomicSequence gs = new GenomicSequence(dna.getBytes());
        gs.setGeneticCode("UNIVERSAL");
        gs.setNoFormat(false);
        gs.setPrintCodon(PrintCodon.STOP);

        gs.setFrames(Frame.getAllFrames());
        System.out.print(gs.getPrintableSequence() + "---\n");
        gs.setFrames(Frame.getForwardFrames());
        System.out.print(gs.getPrintableSequence() + "---\n");
        gs.setFrames(Frame.getReverseFrames());
        System.out.print(gs.getPrintableSequence() + "---\n");

        gs.setFrames(new ArrayList<>());
        System.out.print(gs.getPrintableSequence() + "---\n");
    }
}