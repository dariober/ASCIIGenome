package samTextViewer;

import java.util.LinkedHashMap;

import net.sourceforge.argparse4j.ArgumentParsers;
import net.sourceforge.argparse4j.impl.Arguments;
import net.sourceforge.argparse4j.inf.ArgumentParser;
import net.sourceforge.argparse4j.inf.ArgumentParserException;
import net.sourceforge.argparse4j.inf.Namespace;

public class ArgParse {
	
	public static String PROG_NAME= "ASCIIGenome";
	public static String VERSION= "0.4.0";
	public static String WEB_ADDRESS= "https://github.com/dariober/ASCIIGenome";
	
	public static LinkedHashMap<String, String> docstrings= new LinkedHashMap<String, String>(); 
	
	
	/* Parse command line args */
	public static Namespace argParse(String[] args){
		ArgumentParser parser= ArgumentParsers
				.newArgumentParser(PROG_NAME)
				.defaultHelp(true)
				.version("${prog} " + VERSION)
				.description("DESCRIPTION\n"
+ "Text viewer for genome alignment and annotation files.\n"
+ "For details see " + WEB_ADDRESS + "\n"
+ "");	
		parser.addArgument("input")
			.type(String.class)
			.required(false)
			.nargs("*")
			.help("Input files to be displayed: bam, bed, gtf, bigwig, bedgraph, etc");

		parser.addArgument("--batchFile", "-b")
			.type(String.class)
			.required(false)
			.setDefault("")
			.help("Bed or gff file of regions to process in batch. "
					+ "ASCIIGenome will iterate through the regions in this file");

		parser.addArgument("--region", "-r")
			.type(String.class)
			.required(false)
			.setDefault("")
			.help("Go to region. Format 1-based as 'chrom:start-end' or 'chrom:start' or 'chrom'. E.g. chr1:1-1000");

		parser.addArgument("--genome", "-g")
			.type(String.class)
			.setDefault("")
			.help("A genome file or a tag identifying a genome build (e.g. hg19), or bam file with suitable header");
		
		parser.addArgument("--fasta", "-fa")
			.type(String.class)
			.help("Optional reference fasta file.\n"
					+ "If given, must be indexed, e.g. with `samtools faidx ref.fa`");

		parser.addArgument("--exec", "-x")
			.type(String.class)
			.setDefault("")
			.help("Commands to be executed at the prompt. Must be a single string. E.g. 'goto chr1 && next && seqRegex ACTG'");

		
		//parser.addArgument("--maxReadsStack", "-M")
		//	.type(Integer.class)
		//	.setDefault(2000)
		//	.help("Maximum number of reads to accumulate before printing. If more than this many reads map to the window\n"
		//			+ "randomy select them");
		
		parser.addArgument("--noFormat", "-nf")
			.action(Arguments.storeTrue())
			.help("Do not format output with non ascii chars (colour, bold, etc.)");

		parser.addArgument("--nonInteractive", "-ni")
			.action(Arguments.storeTrue())
			.help("Non interactive mode: Exit after having processed cmd line args.");
		
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
