package commandHelp;

import java.util.LinkedHashMap;
import java.util.Map;

import net.sourceforge.argparse4j.ArgumentParsers;
import net.sourceforge.argparse4j.impl.Arguments;
import net.sourceforge.argparse4j.inf.ArgumentParser;

public class Argparse {

	public static Map<Command, ArgumentParser> parser(){
		
		Map<Command, ArgumentParser> argMap= new LinkedHashMap<Command, ArgumentParser>();
		// ----------------------------------
		ArgumentParser parser= ArgumentParsers
				.newArgumentParser("print")
				.defaultHelp(true)
				.version("Print lines for the tracks matched by `track_regex`.")
				.description("Print stuff");
		
		parser.addArgument("--nlines", "-n").help("Print up to this many lines, default 10. No limit if < 0.").type(Integer.class).setDefault(10);
		parser.addArgument("--full", "-full").help("Print full linesa dn wrap them if wider than the screen.").action(Arguments.storeTrue());
		parser.addArgument("--clip", "-clip").help("Clip lines longer than the screen width. This is the default.")
			.action(Arguments.storeTrue());
		parser.addArgument("--off", "-off").help("Turn off printing.").action(Arguments.storeTrue());
		parser.addArgument("--invert", "-v").help("Invert selection: apply changes to the tracks not selected by list of track_regex.").action(Arguments.storeTrue());
		parser.addArgument("track_regex").nargs("*").setDefault( new String[] {".*"} );
		parser.addArgument("--brief").action(Arguments.version());
		// ----------------------------------
		return argMap;
	}
	
}
