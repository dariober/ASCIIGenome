package commandHelp;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import exceptions.InvalidCommandLineException;
import jline.console.ConsoleReader;
import jline.console.completer.StringsCompleter;

public class CommandList {
	
	private static String SEE_ALSO= "\nSee also and feel free to report issues to https://github.com/dariober/ASCIIGenome\n\n";
	
	public static ConsoleReader initConsole() throws IOException{
		
		ConsoleReader console= new ConsoleReader(); 
		try {
			for(CommandHelp x : CommandList.commandHelpList()){
				if(x.getName().length() > 2){
					console.addCompleter(new StringsCompleter(x.getName()));
				}
			}
		} catch (InvalidCommandLineException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return console;
	}
	
	public static String fullHelp() throws InvalidCommandLineException{
		String help= "\n      N a v i g a t i o n \n\n";
		for(CommandHelp x : CommandList.getCommandsForSection(Section.NAVIGATION)){
			help += (x.printCommandHelp() + "\n");
		}
		help += "\n      F i n d \n\n";
		for(CommandHelp x : CommandList.getCommandsForSection(Section.FIND)){
			help += (x.printCommandHelp() + "\n");
		}
		help += "\n      D i s p l a y \n\n";
		for(CommandHelp x : CommandList.getCommandsForSection(Section.DISPLAY)){
			help += (x.printCommandHelp() + "\n");
		}
		
		help += "\n      A l i g n m e n t s \n\n";
		for(CommandHelp x : CommandList.getCommandsForSection(Section.ALIGNMENTS)){
			help += (x.printCommandHelp() + "\n");
		}

		for(CommandHelp x : CommandList.getCommandsForSection(Section.GENERAL)){
			help += (x.printCommandHelp() + "\n");
		}
		
		help += "\n      G e n e r a l \n\n";
		for(CommandHelp x : CommandList.getCommandsForSection(Section.GENERAL)){
			help += (x.printCommandHelp());
		}
		help += SEE_ALSO;
		return help;
	}
	
	public static String briefHelp() throws InvalidCommandLineException{
		String help= "\n      N a v i g a t i o n \n\n";
		for(CommandHelp x : CommandList.getCommandsForSection(Section.NAVIGATION)){
			help += (x.printBriefHelp());
		}
		help += "\n      F i n d \n\n";
		for(CommandHelp x : CommandList.getCommandsForSection(Section.FIND)){
			help += (x.printBriefHelp());
		}
		help += "\n      D i s p l a y \n\n";
		for(CommandHelp x : CommandList.getCommandsForSection(Section.DISPLAY)){
			help += (x.printBriefHelp());
		}
		
		help += "\n      A l i g n m e n t s \n\n";
		for(CommandHelp x : CommandList.getCommandsForSection(Section.ALIGNMENTS)){
			help += (x.printBriefHelp());
		}
		help += "\n      G e n e r a l \n\n";
		for(CommandHelp x : CommandList.getCommandsForSection(Section.GENERAL)){
			help += (x.printBriefHelp());
		}
		help += SEE_ALSO;
		return help;
	}

	
	private final static List<CommandHelp> commandHelpList() throws InvalidCommandLineException{
		List<CommandHelp> cmdList= new ArrayList<CommandHelp>();
		CommandHelp cmd= new CommandHelp();		

		cmd.setName("f"); cmd.setArgs(""); cmd.inSection= Section.NAVIGATION;
		cmd.setBriefDescription("Move forward by 1/10 of a window");  
		cmdList.add(cmd);
		
		cmd= new CommandHelp(); 
		cmd.setName("b"); cmd.setArgs(""); cmd.inSection= Section.NAVIGATION; 
		cmd.setBriefDescription("Move backward by 1/10 of a window"); 
		cmdList.add(cmd);
		
		cmd= new CommandHelp();
		cmd.setName("ff"); cmd.setArgs(""); cmd.inSection= Section.NAVIGATION;
		cmd.setBriefDescription("Move forward by 1/2 of a window");  
		cmdList.add(cmd);
		
		cmd= new CommandHelp();
		cmd.setName("bb"); cmd.setArgs(""); cmd.inSection= Section.NAVIGATION; 
		cmd.setBriefDescription("Move backward by 1/2 of a window"); 
		cmdList.add(cmd);
		
		cmd= new CommandHelp();
		cmd.setName("goto"); cmd.setArgs("<chrom>:[from]-[to]"); cmd.inSection= Section.NAVIGATION; 
		cmd.setBriefDescription("Go to region chrom:from-to or to chrom:from or to start of chrom. "); 
		cmd.setAdditionalDescription("The character ':' is a shortcut for `goto`. Examples:\n"
				+ "goto chr8:1-1000 \n"
				+ "goto chr8:1 \n"
				+ "goto chr8 \n"
				+ "## Or the same\n"
				+ ":chr8:1-1000 \n"
				+ ":chr8:1 \n"
				+ ":chr8");
		cmdList.add(cmd);
				
		cmd= new CommandHelp();
		cmd.setName("zi"); cmd.setArgs("[INT = 1]"); cmd.inSection= Section.NAVIGATION; 
		cmd.setBriefDescription("Zoom in INT times. Each zoom halves the window size."); 
		cmd.setAdditionalDescription("To zoom quickly use INT=~5 or 10 e.g. `zi~10`");
		cmdList.add(cmd);
		
		cmd= new CommandHelp();
		cmd.setName("zo"); cmd.setArgs("[INT = 1]"); cmd.inSection= Section.NAVIGATION; 
		cmd.setBriefDescription("Zoom out INT times. Each zoom doubles the window size.");
		cmd.setAdditionalDescription("To zoom quickly use INT=~5 or 10 e.g. `zo 10`");
		cmdList.add(cmd);
		
		cmd= new CommandHelp();
		cmd.setName("INT"); cmd.setArgs("[INT]"); cmd.inSection= Section.NAVIGATION; 
		cmd.setBriefDescription(""
				+ "Go to position `INT` or to region `INT INT` on current chromosome.");
		cmd.setAdditionalDescription(""
				+ "Allowed is the hyphenated format  separating the two positions. "
				+ "If a list of integers is given, the first and last are taken as *from* and *to*. "
				+ "This is handy to copy and paste intervals from the ruler above the prompt. "
				+ "\nExamples:\n"
				+ "10~~~~~~~~~~~~~~~~~~~## Will jump to position 10 \n"
				+ "10 1000~~~~~~~~~~~~~~## Go to region 10-1000 \n"
				+ "10-1000~~~~~~~~~~~~~~## Same as above\n"
				+ "10 250 500 750 1000~~## Same as above again");
		cmdList.add(cmd);
		
		cmd= new CommandHelp();
		cmd.setName("+"); cmd.setArgs("<INT> [k|m]"); cmd.inSection= Section.NAVIGATION; 
		cmd.setBriefDescription("Move forward by INT bases. Suffixes k (kilo) and M (mega) are expanded to x1000 and x1,000,000. "
				+ "Examples: `-2m` or `+10k` or `+10.5k`"); 
		cmdList.add(cmd);
		
		cmd= new CommandHelp();
		cmd.setName("-"); cmd.setArgs("<INT> [k|m]"); cmd.inSection= Section.NAVIGATION; 
		cmd.setBriefDescription("Move backwards by INT bases. Suffixes k (kilo) and M (mega) are expanded to x1000 and x1,000,000.\n"
				+ "Examples: `-100` or `-10k` or `-10.5m`"); 
		cmdList.add(cmd);
		
		cmd= new CommandHelp();
		cmd.setName("p"); cmd.setArgs(""); cmd.inSection= Section.NAVIGATION; 
		cmd.setBriefDescription("Go to the previous visited position.");
		cmd.setAdditionalDescription("Similar to the back and forward arrows of an Internet browser.");
		cmdList.add(cmd);

		cmd= new CommandHelp();
		cmd.setName("n"); cmd.setArgs(""); cmd.inSection= Section.NAVIGATION; 
		cmd.setBriefDescription("Go to the next visited position.");
		cmd.setAdditionalDescription("Similar to the back and forward arrows of an Internet browser.");
		cmdList.add(cmd);
		

		cmd= new CommandHelp();
		cmd.setName("next"); cmd.setArgs("[track_id]"); cmd.inSection= Section.NAVIGATION; 
		cmd.setBriefDescription("Move to the next feature on track_id on *current* chromosome. "); 
		cmd.setAdditionalDescription(""
				+ "`next` centers the window on the found feature and zooms out. "
				+ "This is useful for quickly browsing through annotation files of genes or ChIP-Seq "
				+ "peaks in combination with read coverage tracks (bigwig, tdf, etc.). "
				+ "`next_start` instead sets the window right at the start of the feature.\n "
				+ "\n"
				+ "The `next` command does exactly that, it moves to the next feature. "
				+ "If there are no more features after the current position it doesn't rewind to the beginning "
				+ "(use `1` for that) and it doesn't move to another chromosome, "
				+ "use `goto chrom` for that.\n "
				+ "\n"
				+ "If `track_id` is omitted, the first annotation track is used. "
				+ "If trackId is not a feature track (bed, gtf, etc) a more or less ugly warning is issued."); 

		cmdList.add(cmd);

		cmd= new CommandHelp();
		cmd.setName("next_start"); cmd.setArgs("[track_id]"); cmd.inSection= Section.NAVIGATION; 
		cmd.setBriefDescription("Move to the next feature on track_id on *current* chromosome. "); 
		cmd.setAdditionalDescription(""
				+ "`next` centers the window on the found feature and zooms out. "
				+ "This is useful for quickly browsing through annotation files of genes or ChIP-Seq "
				+ "peaks in combination with read coverage tracks (bigwig, tdf, etc.). "
				+ "`next_start` instead sets the window right at the start of the feature.\n "
				+ "\n"
				+ "The `next` command does exactly that, it moves to the next feature. "
				+ "If there are no more features after the current position it doesn't rewind to the beginning "
				+ "(use `1` for that) and it doesn't move to another chromosome, "
				+ "use `goto chrom` for that.\n "
				+ "\n"
				+ "If `track_id` is omitted, the first annotation track is used. "
				+ "If track_id is not a feature track (bed, gtf, etc) a more or less ugly warning is issued."); 

		cmdList.add(cmd);
		
		cmd= new CommandHelp();
		cmd.setName("find_first"); cmd.setArgs("<regex> [track_id]"); cmd.inSection= Section.FIND; 
		cmd.setBriefDescription("Find the first record in track_id containing regex."); 
		cmd.setAdditionalDescription(""
				+ "The search starts from the *end* of the current window "
				+ "(so the current window is not searched) and moves forward on the current chromosome. "
				+ "At the end  of the current chromosome move to the next chromosomes and then restart at "
				+ " the start of the initial one. The search stops at the first match found."); 
		cmdList.add(cmd);

		cmd= new CommandHelp();
		cmd.setName("find_all"); cmd.setArgs("<regex> [track_id]"); cmd.inSection= Section.FIND; 
		cmd.setBriefDescription("Find the region containing *all* the records on chromosome containing regex. "); 
		cmd.setAdditionalDescription(""
				+ "The search starts at the current chromosome before moving to the other ones. "
				+ "It stops at the first chromosome returning one or more hits. "
				+ "Useful to get all gtf records of a gene.\n"
				+ "\n"
				+ "E.g. `find_all ACTB genes.gtf` will find the entire ACTB gene "
				+ "(provided the regex is specific enough of course)."); 

		cmdList.add(cmd);

		cmd= new CommandHelp();
		cmd.setName("seqRegex"); cmd.setArgs("<regex>"); cmd.inSection= Section.FIND; 
		cmd.setBriefDescription("Find <regex> in reference sequence and show matches as and additional track. ");
		cmd.setAdditionalDescription(""
				+ "Useful to show restriction enzyme sites, "
				+ "transcription factor motifs, etc. The tag of this track is "
				+ "`seqRegex` and it is not displayed. To adjust its height use "
				+ "`trackHeight~10~seqRegex`. If regex is omitted the matching is disabled Matching is case sensitive, to ignore case use the regex syntax "
				+ "`(?i)`. Example\n"
				+ "seqReg~(?i)ACTG\n"
				+ "This command is ignored if the reference fasta file is missing.");
		cmdList.add(cmd);
		
		cmd= new CommandHelp();
		cmd.setName("visible"); cmd.setArgs("[show_regex = .*] [hide_regex = ''] [track_regex = .*]..."); cmd.inSection= Section.DISPLAY; 
		cmd.setBriefDescription("Display features matching show_regex, hide those matching hide_regex. Apply to tracks matched by track_regex.");
		cmd.setAdditionalDescription(""
				+ "\nThis command is useful to filter the annotation in GTF or BED files, for example: "
				+ "`visible~RNA~mRNA~gtf` \n "
				+ "Will show the rows containing 'RNA' but will hide those containing 'mRNA', applies "
				+ "to tracks whose name matches 'gtf'."
				+ ""
				+ "\nWith no arguments reset to default: `visible~.*~^$~.*` which means show everything, hide nothing, apply to all tracks. "
				);
		cmdList.add(cmd);		

		cmd= new CommandHelp();
		cmd.setName("squash"); cmd.setArgs("[track_regex = .*]..."); cmd.inSection= Section.DISPLAY; 
		cmd.setBriefDescription("Toggle the squashing of features with the same coordinates. ");
		cmd.setAdditionalDescription("If set, features with the same start, end, and strand are squashed in a single one. "
				+ "The displayed feature is the first one found in the group of features with the same coordinates. "
				+ "Useful to compact GTF where e.g. CDS and exons have the same coordinates. "
				+ "Applies only to annotation tracks captured by track_regex");
		cmdList.add(cmd);		

		cmd= new CommandHelp();
		cmd.setName("merge"); cmd.setArgs("[track_regex = .*]..."); cmd.inSection= Section.DISPLAY; 
		cmd.setBriefDescription("Toggle the merging of overlapping features.");
		cmd.setAdditionalDescription("If set, features with overalapping coordinates are merged in a single one. "
				+ "Merged features will not have strand and name information. Note that merging is done without considering strand information. "
				+ "Applies only to annotation tracks captured by the list of track_regex");
		cmdList.add(cmd);		
		
		cmd= new CommandHelp();
		cmd.setName("gffNameAttr"); cmd.setArgs("[attribute_name = NULL] [track_regex = .*]..."); cmd.inSection= Section.DISPLAY; 
		cmd.setBriefDescription("For GTF/GFF tracks, choose the attribute to get the feature name from.");
		cmd.setAdditionalDescription("Use attribute NULL to reset to default choice of attribute. "
				+ "Applies to all GFF/GTF tracks captured by the list of `track_regex`. Example\n"
				+ "gffNameAttr gene_name genes.gtf .*gff");
		cmdList.add(cmd);

		cmd= new CommandHelp();
		cmd.setName("trackHeight"); cmd.setArgs("INT [track_regex = .*]..."); cmd.inSection= Section.DISPLAY; 
		cmd.setBriefDescription("Set track height to INT lines of text for all tracks matching regexes.\n");
		cmd.setAdditionalDescription("Example: `trackHeight 5 aln.*bam gtf`");
		cmdList.add(cmd);

		cmd= new CommandHelp();
		cmd.setName("ylim"); cmd.setArgs("<min> <max> [track_regex = .*]..."); cmd.inSection= Section.DISPLAY; 
		cmd.setBriefDescription("Set the y-axis limit for all tracks matched by regexes.");
		cmd.setAdditionalDescription("Use `na` to autoscale to min and/or max. "
				+ "This command applies only to tracks displaying quantitative data on y-axis (e.g. bigwig, tdf), "
				+ "the other tracks are unaffected.\n"
				+ "Examples:\n "
				+ "ylim 0 50~~~~~~## Set min= 0 and max= 50 in all tracks.\n"
				+ "ylim 0 na~~~~~~## Set min to 0 and autoscale the max. Apply to all tracks\n"
				+ "ylim na na tdf~## Autoscale min and max. Apply to all tracks matching 'tdf'");
		cmdList.add(cmd);

		cmd= new CommandHelp();
		cmd.setName("colorTrack"); cmd.setArgs("<color> [track_regex = .*]..."); cmd.inSection= Section.DISPLAY; 
		cmd.setBriefDescription("Set colour for tracks matched by regex. ");
		cmd.setAdditionalDescription(""
				+ "Available colours: red, green, yellow, blue, magenta, cyan, grey, "
				+ "light_red, light_green, light_yellow, light_blue, light_magenta, light_cyan, light_grey, "
				+ "white, black, default. The 'default' colour reset to the system default colour. "
				+ "Colouring is rendered with ANSI codes 8/16. Example:\n"
				+ "colorTrack~light_blue~ts.*gtf~ts.*bam");
		cmdList.add(cmd);

		cmd= new CommandHelp();
		cmd.setName("dataCol"); cmd.setArgs("[index = 4] [track_regex = .*]..."); cmd.inSection= Section.DISPLAY; 
		cmd.setBriefDescription("Select data column for bedgraph tracks containing regex. ");
		cmd.setAdditionalDescription("index: 1-based column index. This command applies only to "
				+ "tracks of type bedgraph.\n For example, use column 5 on tracks containing #1 and #3:\n "
				+ "dataCol 5 #1 #3");
		cmdList.add(cmd);

		cmd= new CommandHelp();
		cmd.setName("print"); cmd.setArgs("[track_regex = .*]..."); cmd.inSection= Section.DISPLAY; 
		cmd.setBriefDescription("Toggle the printing of lines in the tracks matched by track_regex. Long lines clipped. ");
		cmd.setAdditionalDescription("Useful to show exactly what features are present in the current window. "
				+ "Features are filtered in/out according to the `visible` command. Applies only to annotation tracks");
		cmdList.add(cmd);

		cmd= new CommandHelp();
		cmd.setName("printFull"); cmd.setArgs("[track_regex = .*]..."); cmd.inSection= Section.DISPLAY; 
		cmd.setBriefDescription("Toggle the printing of lines in the tracks matched by track_regex. Long lines wrapped. ");
		cmd.setAdditionalDescription("Useful to show exactly what features are present in the current window. "
				+ "Features are filtered in/out according to the `visible` command. Applies only to annotation tracks");
		cmdList.add(cmd);

		cmd= new CommandHelp();
		cmd.setName("showGenome"); cmd.setArgs(""); cmd.inSection= Section.DISPLAY; 
		cmd.setBriefDescription("Print the genome dictionary with a representation of chromosome sizes.");
		cmd.setAdditionalDescription("");
		cmdList.add(cmd);

		cmd= new CommandHelp();
		cmd.setName("addTracks"); cmd.setArgs("[file or URL]..."); cmd.inSection= Section.DISPLAY; 
		cmd.setBriefDescription("Add tracks from local or remote files.");
		cmd.setAdditionalDescription("");
		cmdList.add(cmd);
		
		cmd= new CommandHelp();
		cmd.setName("orderTracks"); cmd.setArgs("[track_regex]..."); cmd.inSection= Section.DISPLAY; 
		cmd.setBriefDescription("Reorder tracks according to the list of regexes.");
		cmd.setAdditionalDescription("Not all the tracks need to be listed, the missing ones "
				+ "follow the listed ones in unchanged order.\n"
				+ "For example, given the track list: `[hela.bam#1, hela.bed#2, hek.bam#3, hek.bed#4]`:\n\n"
				+ "orderTracks #2 #1   -> [hela.bed#2, hela.bam#1, hek.bam#3, hek.bed#4]\n"
				+ "orderTracks bam bed -> [hela.bam#1, hek.bam#3, hela.bed#2, hek.bed#4]");
		cmdList.add(cmd);


		cmd= new CommandHelp();
		cmd.setName("history"); cmd.setArgs(""); cmd.inSection= Section.DISPLAY; 
		cmd.setBriefDescription("Show the list of visited positions.");
		cmd.setAdditionalDescription("");
		cmdList.add(cmd);
		
		cmd= new CommandHelp();
		cmd.setName("rpm"); cmd.setArgs("[track_regex = .*]"); cmd.inSection= Section.ALIGNMENTS; 
		cmd.setBriefDescription("Toggle read coverage from raw count to reads per million.");
		cmd.setAdditionalDescription("");
		cmdList.add(cmd);		

		cmd= new CommandHelp();
		cmd.setName("-f"); cmd.setArgs("INT [track_regex = .*]..."); cmd.inSection= Section.ALIGNMENTS; 
		cmd.setBriefDescription("Include reads with INT bits set in tracks matched by regexes.");
		cmd.setAdditionalDescription("Same as in `samtools view`. Note that the flag 4096 "
				+ "can be used to filter in or out reads on top strand, this is useful in bisulfite mode.");
		cmdList.add(cmd);		

		cmd= new CommandHelp();
		cmd.setName("-F"); cmd.setArgs("INT [track_regex = .*]..."); cmd.inSection= Section.ALIGNMENTS; 
		cmd.setBriefDescription("Exclude reads with INT bits set in tracks matched by regexes.");
		cmd.setAdditionalDescription("Same as in `samtools view`. Note that the flag 4096 "
				+ "can be used to filter in or out reads on top strand, this is useful in bisulfite mode.");
		cmdList.add(cmd);

		cmd= new CommandHelp();
		cmd.setName("mapq"); cmd.setArgs("INT [track_regex = .*]..."); cmd.inSection= Section.ALIGNMENTS; 
		cmd.setBriefDescription("Include reads with mapq >= INT in tracks matched by regexes");
		cmd.setAdditionalDescription("For example: mapq 30 aln1.bam aln2.bam");
		cmdList.add(cmd);
		
		cmd= new CommandHelp();
		cmd.setName("BSseq"); cmd.setArgs("[track_regex = .*]..."); cmd.inSection= Section.ALIGNMENTS; 
		cmd.setBriefDescription("Toggle bisulfite mode for read tracks matched by regex.");
		cmd.setAdditionalDescription("In bisulfite mode, the characters M and m mark methylated bases "
				+ "(i.e. unconverted C to T) and U and u are used for unmethylated bases "
				+ "(i.e. C converted to T). Upper case is used for reads on  forward strand, small case for reverse. "
				+ "Ignored without reference fasta sequence.");
		cmdList.add(cmd);		

		cmd= new CommandHelp();
		cmd.setName("save"); cmd.setArgs("[filename = chrom_start_end.txt']"); cmd.inSection= Section.GENERAL; 
		cmd.setBriefDescription("Save current screenshot to file in either text or png format.");
		cmd.setAdditionalDescription("Default filename is generated from the current coordinates and the default format is txt. "
				+ "With filename .png save as png using current coordinates as filename. "
				+ "Use extension .png to save as png format. Note that colours are not retained.");
		cmdList.add(cmd);

		cmd= new CommandHelp();
		cmd.setName("q"); cmd.setArgs(""); cmd.inSection= Section.GENERAL; 
		cmd.setBriefDescription("Quit");
		cmd.setAdditionalDescription("");
		cmdList.add(cmd);

		cmd= new CommandHelp();
		cmd.setName("h"); cmd.setArgs(""); cmd.inSection= Section.GENERAL; 
		cmd.setBriefDescription("Show this help. For help on specific commands use `<command name> -h`");
		cmd.setAdditionalDescription("");
		cmdList.add(cmd);
		
		// Make sure ther are no undocumented cmds
		List<String> documented= new ArrayList<String>();
		for(CommandHelp x : cmdList){
			if(documented.contains(x.getName())){
				System.err.println(x.getName() + " already documented!");
				throw new InvalidCommandLineException();
			}
			documented.add(x.getName());
		}
		for(String x : CommandList.cmds()){
			if(!documented.contains(x)){
				System.err.println("Undocumented command: " + x);
				// throw new InvalidCommandLineException();
			}
		}
		
		return cmdList;
			
		}

	protected static List<CommandHelp> getCommandsForSection(Section section) throws InvalidCommandLineException{
		List<CommandHelp> cmdList= new ArrayList<CommandHelp>();
		for(CommandHelp x : commandHelpList()){
			if(x.inSection.equals(section)){
				cmdList.add(x);
			}
		}
		return cmdList;
	}
	
	
	/* Known commnds */
	protected static final List<String> cmds(){
		List<String> paramList= new ArrayList<String>();
		paramList.add("q");
		paramList.add("h");
		paramList.add("f");
		paramList.add("b");
		paramList.add("ff");
		paramList.add("bb");
		paramList.add("zi");
		paramList.add("zo");
		paramList.add("goto");
		paramList.add("INT");
		paramList.add("+");
		paramList.add("-");
		paramList.add("p");
		paramList.add("n");
		paramList.add("next");
		paramList.add("next_start");
		paramList.add("find_first");
		paramList.add("find_all");
		paramList.add("seqRegex");
		paramList.add("visible");
		paramList.add("gffNameAttr");
		paramList.add("squash");
		paramList.add("merge");
		paramList.add("trackHeight");
		paramList.add("colorTrack");
		paramList.add("ylim");
		paramList.add("dataCol");
		paramList.add("print");
		paramList.add("printFull");
		paramList.add("showGenome");
		paramList.add("addTracks");
		paramList.add("orderTracks");
		paramList.add("history");
		paramList.add("rpm");
		paramList.add("-f");
		paramList.add("-F");
		paramList.add("mapq");
		paramList.add("BSseq");
		paramList.add("save");
		//paramList.add("savef");
	
		return paramList;
	}

	public static String getHelpForCommand(String commandName) {
		try {
			for(CommandHelp x : CommandList.commandHelpList()){
				if(x.getName().equals(commandName)){
					return x.printCommandHelp();
				}
			}
		} catch (InvalidCommandLineException e1) {
			e1.printStackTrace();
		}
		return "";
	}

}
