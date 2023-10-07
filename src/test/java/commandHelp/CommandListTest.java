package commandHelp;

import java.io.File;
import java.io.IOException;

import org.junit.Test;

import coloring.Config;
import exceptions.InvalidColourException;
import exceptions.InvalidCommandLineException;
import exceptions.InvalidConfigException;

public class CommandListTest {

        @Test
    public void updateReStructuredFile() throws InvalidCommandLineException, IOException, InvalidColourException{
        CommandList.updateCommandHelpMdFile(new File("docs/commandHelp.rst"));
    }
    
    @Test
    public void canPrintBriefHelp() throws InvalidCommandLineException, InvalidColourException {
        System.out.println("BRIEF HELP");
        System.out.println(CommandList.briefHelp());
    }
    
    @Test
    public void canPrintDocstringForCommand() throws InvalidCommandLineException, IOException, InvalidColourException, InvalidConfigException {
        new Config(null);
        System.out.println("DOCSTRING");
        System.out.println(CommandList.getHelpForCommand("featureColor"));
    }
    
}
