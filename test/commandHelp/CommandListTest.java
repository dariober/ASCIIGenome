package commandHelp;

import static org.junit.Assert.*;

import org.junit.Test;

import exceptions.InvalidCommandLineException;

public class CommandListTest {

	@Test
	public void printMarkdownDocs() throws InvalidCommandLineException {
		System.out.println(CommandList.markdownDocs());
	}
	
	//@Test
	public void canPrintFullHelp() throws InvalidCommandLineException {
		System.out.println("FULL HELP");
		System.out.println(CommandList.fullHelp());
	}

	//@Test
	public void canPrintBriefHelp() throws InvalidCommandLineException {
		System.out.println("BRIEF HELP");
		System.out.println(CommandList.briefHelp());
	}
	
	//@Test
	public void canPrintDocstringForCommand() throws InvalidCommandLineException {
		System.out.println("DOCSTRING");
		System.out.println(CommandList.getHelpForCommand("ylim"));
	}

}
