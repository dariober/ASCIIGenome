package commandHelp;

import static org.junit.Assert.*;

import java.util.List;

import org.junit.Test;

import com.google.common.base.Splitter;

import exceptions.InvalidCommandLineException;
import samTextViewer.Utils;

public class CommandHelpTest {

    @Test
    public void testSplitter() throws InvalidCommandLineException {
        Iterable<String> cmdInputList= Splitter.on("&&").trimResults().omitEmptyStrings().split("zi && zo && ylim ");
        for(String x : cmdInputList){
            List<String> cmdSubInput= Utils.tokenize(x, " ");
            System.out.println(cmdSubInput);
        }
        
        cmdInputList= Splitter.on("&&").trimResults().omitEmptyStrings().split("");
        for(String x : cmdInputList){
            List<String> cmdSubInput= Utils.tokenize(x, " ");
            System.out.println(cmdSubInput);
        }
        
    }
    
    @Test
    public void canPrintHelpString() throws InvalidCommandLineException {

        CommandHelp cmd= new CommandHelp();
        cmd.setName("ylim");
        cmd.setArgs("<min> <max> [track_re = .*]...");
        cmd.setBriefDescription("This is the very short description. TYpically one or two lines.");
        cmd.setAdditionalDescription("Set limits of y axis for all track IDs containing regex. (This~strecth~of~words~in~brackets~should~stay~together.) "
                + "New line here -> \n "
                + "ylim 0 na. If regex is omitted all tracks will be captured. Default: 'ylim na na .*' "
                + "Set limits of y axis for all track IDs containing regex. Use na to autoscale to min and/or max. E.g. ylim 0 na. If regex is omitted all tracks will be captured. Default: 'ylim na na .*'");
        
        System.out.print(cmd.printCommandHelp());
        
        System.out.print(cmd.printBriefHelp());
        
    }

}
