package samTextViewer;

import java.util.LinkedHashMap;

import net.sourceforge.argparse4j.ArgumentParsers;
import net.sourceforge.argparse4j.impl.Arguments;
import net.sourceforge.argparse4j.inf.ArgumentParser;
import net.sourceforge.argparse4j.inf.ArgumentParserException;
import net.sourceforge.argparse4j.inf.Namespace;

public class ArgParse {
	
	public static String PROG_NAME= "ASCIIGenome";
	public static String VERSION= "1.11.0";
	public static String WEB_ADDRESS= "https://github.com/dariober/ASCIIGenome";
	public static String WEB_RTD= "http://asciigenome.readthedocs.io/";
	
	public static LinkedHashMap<String, String> docstrings= new LinkedHashMap<String, String>(); 
	
	
	/* Parse command line args */
	public static Namespace argParse(String[] args){
		ArgumentParser parser= ArgumentParsers
				.newArgumentParser(PROG_NAME)
				.defaultHelp(true)
				.version("${prog} " + VERSION)
				.description("DESCRIPTION\n"
+ "Genome browser at the command line.\n"
+ "\nFor details see " + WEB_RTD);	
		parser.addArgument("input")
			.type(String.class)
			.required(false)
			.nargs("*")
			.help("Input files to be displayed: bam, bed, gtf, bigwig, bedgraph, etc");

		parser.addArgument("--batchFile", "-b")
			.type(String.class)
			.required(false)
			.help("Bed or gff file of regions to process in batch. Use - to read from stdin. "
					+ "ASCIIGenome will iterate through the regions in this file");

		parser.addArgument("--region", "-r")
			.type(String.class)
			.required(false)
			.help("Go to region. Format 1-based as 'chrom:start-end' or 'chrom:start' or 'chrom'. E.g. chr1:1-1000");

//		parser.addArgument("--genome", "-g")
//			.type(String.class)
//			.help("A genome file or a tag identifying a genome build (e.g. hg19), or bam file with suitable header");
		
		parser.addArgument("--fasta", "-fa")
			.type(String.class)
			.help("Optional reference fasta file.\n"
					+ "If given, must be indexed, e.g. with `samtools faidx ref.fa`");

		parser.addArgument("--exec", "-x")
			.type(String.class)
			.help("Commands to be executed at the prompt. Either a single string, e.g. 'goto chr1 && next && seqRegex ACTG' or a file with one command per line.");

		parser.addArgument("--noFormat", "-nf")
			.action(Arguments.storeTrue())
			.help("Do not format output with non ascii chars (colour, bold, etc.)");

		parser.addArgument("--nonInteractive", "-ni")
			.action(Arguments.storeTrue())
			.help("Non interactive mode: Exit after having processed cmd line args.");
		
		parser.addArgument("--config", "-c")
			.type(String.class)
			.required(false)
			.setDefault("metal")
			.help("Source of configuration settings. "
					+ "It can be a local file or a tag matching a built-in configuration: "
					+ "'black_on_white', 'white_on_black', 'metal'. "
					+ "If null, first try to read configuration from file '~/.asciigenome_config'. "
					+ "If this file is missing use a built-in setting. "
					+ "For examples of configuration files see https://github.com/dariober/ASCIIGenome/blob/master/resources/config/");
		
		parser.addArgument("--debug")
			.type(Integer.class)
			.choices(0, 1, 2)
			.setDefault(0)
			.help("Set debugging mode. 0: off; 1: print exception stack traces; 2: print stack traces and exit.");
		
		parser.addArgument("--version", "-v").action(Arguments.version());
		
		Namespace opts= null;
		try{
			opts= parser.parseArgs(args);
		}
		catch(ArgumentParserException e) {
			parser.handleError(e);
			System.exit(1);
		}		
		return(opts);
	}
	
	public static String getDocstrings(){
		String docstring= "";
		for(String h : docstrings.keySet()){
			docstring += h + "\n    " + docstrings.get(h) + "\n";
		}
		return docstring;
	}
	
}
