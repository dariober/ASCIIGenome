package samTextViewer;

import static org.junit.Assert.*;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.io.FileUtils;
import org.junit.Test;

import com.itextpdf.text.DocumentException;

import exceptions.BamIndexNotFoundException;
import exceptions.InvalidColourException;
import exceptions.InvalidCommandLineException;
import exceptions.InvalidConfigException;
import exceptions.InvalidGenomicCoordsException;
import exceptions.InvalidRecordException;
import faidx.UnindexableFastaFileException;

public class MainTest {

	@Test
	public void canSuggestCommand() throws ClassNotFoundException, IOException, InvalidGenomicCoordsException, InvalidCommandLineException, InvalidRecordException, BamIndexNotFoundException, SQLException, DocumentException, UnindexableFastaFileException, InvalidColourException, InvalidConfigException {
		String[] args= new String[] {"-ni", "-nf", "--exec", "prnt"};
		List<String> out = this.runMain(args);
		assertTrue(out.get(1).contains("Maybe you mean print?"));
	}

	@Test
	public void canSetConfig() throws ClassNotFoundException, IOException, InvalidGenomicCoordsException, InvalidCommandLineException, InvalidRecordException, BamIndexNotFoundException, SQLException, DocumentException, UnindexableFastaFileException, InvalidColourException, InvalidConfigException {
		String[] args= new String[] {"-ni", "-nf", "--exec", "setConfig nucs f", "test_data/ds051.short.bam"};
		List<String> out = this.runMain(args);
		assertTrue(out.get(0).contains(">>>>>>>>>>>>>>>>>>>"));
	}
	
	@Test
	public void canFlipBooleanConfig() throws ClassNotFoundException, IOException, InvalidGenomicCoordsException, InvalidCommandLineException, InvalidRecordException, BamIndexNotFoundException, SQLException, DocumentException, UnindexableFastaFileException, InvalidColourException, InvalidConfigException {
		String[] args= new String[] {"-ni", "-nf", "--exec", "setConfig nucs", "test_data/ds051.short.bam"};
		List<String> out = this.runMain(args);
		assertTrue(out.get(0).contains(">>>>>>>>>>>>>>>>>>>"));
	}

	@Test
	public void doNotSetInvalidBooleanConfig() throws ClassNotFoundException, IOException, InvalidGenomicCoordsException, InvalidCommandLineException, InvalidRecordException, BamIndexNotFoundException, SQLException, DocumentException, UnindexableFastaFileException, InvalidColourException, InvalidConfigException {
		String[] args= new String[] {"-ni", "-nf", "--exec", "setConfig nucs 999", "test_data/ds051.short.bam"};
		List<String> out = this.runMain(args);
		System.err.println(out);
		assertTrue(out.get(1).contains("Unable to set"));
	}
	
	@Test
	public void doNotSetInvalidColourConfig() throws ClassNotFoundException, IOException, InvalidGenomicCoordsException, InvalidCommandLineException, InvalidRecordException, BamIndexNotFoundException, SQLException, DocumentException, UnindexableFastaFileException, InvalidColourException, InvalidConfigException {
		String[] args= new String[] {"-ni", "-nf", "--exec", "setConfig seq_a 999", "test_data/ds051.short.bam"};
		List<String> out = this.runMain(args);
		System.err.println(out);
		assertTrue(out.get(1).contains("Unable to set"));
	}
	
	@Test
	public void doNotSetInvalidIntegerConfig() throws ClassNotFoundException, IOException, InvalidGenomicCoordsException, InvalidCommandLineException, InvalidRecordException, BamIndexNotFoundException, SQLException, DocumentException, UnindexableFastaFileException, InvalidColourException, InvalidConfigException {
		String[] args= new String[] {"-ni", "-nf", "--exec", "setConfig shade_baseq foo", "test_data/ds051.short.bam"};
		List<String> out = this.runMain(args);
		System.err.println(out);
		assertTrue(out.get(1).contains("Unable to set"));
	}
	
	/* H E L P E R S */
	
	/** Execute main with the given array of arguments and return a list of length 2 containing 1) stdout and 2) stderr.
	 * */
	private List<String> runMain(String[] args) throws ClassNotFoundException, IOException, InvalidGenomicCoordsException, InvalidCommandLineException, InvalidRecordException, BamIndexNotFoundException, SQLException, DocumentException, UnindexableFastaFileException, InvalidColourException, InvalidConfigException{

		PrintStream stdout= System.out;
		ByteArrayOutputStream baosOut= new ByteArrayOutputStream();
		System.setOut(new PrintStream(baosOut));
		
		PrintStream stderr= System.err;
		ByteArrayOutputStream baosErr= new ByteArrayOutputStream();
		System.setErr(new PrintStream(baosErr));

		Main.main(args);

		String out= baosOut.toString();
	    System.setOut(stdout);

		String err= baosErr.toString();
	    System.setErr(stderr);

	    List<String> outErr= new ArrayList<String>();
	    outErr.add(out);
	    outErr.add(err);
	    return outErr;
	}
}
