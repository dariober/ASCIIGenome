package commandHelp;

import static org.junit.Assert.*;

import java.io.File;
import java.io.IOException;

import org.junit.Test;

import exceptions.InvalidCommandLineException;

public class CommandListTest {

	@Test
	public void updateMarkdownFile() throws InvalidCommandLineException, IOException{
		CommandList.updateCommandHelpMdFile(new File("commandHelp.md"));
	}
	
	@Test
	public void canPrintFullHelp() throws InvalidCommandLineException {
		System.out.println("FULL HELP");
		System.out.println(CommandList.fullHelp());
	}

	@Test
	public void canPrintBriefHelp() throws InvalidCommandLineException {
		System.out.println("BRIEF HELP");
		System.out.println(CommandList.briefHelp());
	}
	
	@Test
	public void canPrintDocstringForCommand() throws InvalidCommandLineException {
		System.out.println("DOCSTRING");
		System.out.println(CommandList.getHelpForCommand("ylim"));
	}

}
