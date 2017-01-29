package commandHelp;

import java.io.File;
import java.io.IOException;

import org.junit.Test;

import exceptions.InvalidCommandLineException;

public class CommandListTest {

	@Test
	public void updateReStructuredFile() throws InvalidCommandLineException, IOException{
		CommandList.updateCommandHelpMdFile(new File("docs/commandHelp.rst"));
	}
	
	@Test
	public void canPrintBriefHelp() throws InvalidCommandLineException {
		System.out.println("BRIEF HELP");
		System.out.println(CommandList.briefHelp());
	}
	
	@Test
	public void canPrintDocstringForCommand() throws InvalidCommandLineException {
		System.out.println("DOCSTRING");
		System.out.println(CommandList.getHelpForCommand("history"));
	}

}
