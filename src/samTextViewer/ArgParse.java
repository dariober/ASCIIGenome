package samTextViewer;

import java.util.LinkedHashMap;

import net.sourceforge.argparse4j.ArgumentParsers;
import net.sourceforge.argparse4j.impl.Arguments;
import net.sourceforge.argparse4j.inf.ArgumentParser;
import net.sourceforge.argparse4j.inf.ArgumentParserException;
import net.sourceforge.argparse4j.inf.Namespace;

public class ArgParse {
	
	public static String PROG_NAME= "SamTextViewer";
	public static String VERSION= "0.1.0";
	
	public static LinkedHashMap<String, String> docstrings= new LinkedHashMap<String, String>(); 
	
	
	/* Parse command line args */
	public static Namespace argParse(String[] args){
		ArgumentParser parser= ArgumentParsers
				.newArgumentParser(PROG_NAME)
				.defaultHelp(true)
				.version("${prog} " + VERSION)
				.description("DESCRIPTION\n"
+ "Text viewer for genome alignment and annotation files.\n"
+ "For details see https://github.com/dariober/Java-cafe/tree/master/SamTextViewer\n"
+ "Example\n"
+ "java /.../SamTextViewer.jar -r chr18:1000-2000 reads.bam ann.gtf.gz"
+ "");	
		parser.addArgument("insam")
			.type(String.class)
			.required(false)
			.nargs("*")
			.help("Input files. bam/cram must be sorted and indexed. Large bed/gtf files should be indexed with tabix.");
		
		parser.addArgument("--region", "-r")
			.type(String.class)
			.required(false)
			.setDefault("")
			.help("Go to region. Format 1-based as 'chrom:start-end' or 'chrom:start' or 'chrom'. E.g. chr1:1-1000");

		parser.addArgument("--genome", "-g")
			.type(String.class)
			.setDefault("")
			.help("A genome file or a tag identifying a genome build (e.g. hg19), or bam file with suitable header");
		
//		parser.addArgument("--windowSize", "-w")
//			.type(Integer.class)
//			.setDefault(-1)
//			.help("Window size to display. Default to terminal width x0.95");
		
		parser.addArgument("--fasta", "-fa")
			.type(String.class)
			.help("Optional reference fasta reference file.\n"
					+ "If given, must be indexed, e.g. with `samtools faidx ref.fa`");

		docstrings.put("-f", "Required sam flags. Use 4096 for reads on top strand");
		parser.addArgument("--f", "-f")
			.type(Integer.class)
			.setDefault(0)
			.help(docstrings.get("-f"));
		
		docstrings.put("-F", "Filtering sam flags. Use 4096 for reads on top strand");
		parser.addArgument("--F", "-F")
			.type(Integer.class)
			.setDefault(0)
			.help(docstrings.get("-f"));

		docstrings.put("-q", "Minumum mapping quality for a read to be considered");
		parser.addArgument("--mapq", "-q")
			.type(Integer.class)
			.setDefault(0)
			.help(docstrings.get("-q"));
		
		docstrings.put("-m", "Maximum number of lines to print for read tracks.");
		parser.addArgument("--maxLines", "-m")
			.type(Integer.class)
			.setDefault(10)
			.help(docstrings.get("-m"));

		docstrings.put("-rpm", "Toggle on/off the normalization of Reads Per Million for bam input. Default off");
		parser.addArgument("--rpm", "-rpm")
			.action(Arguments.storeTrue())
			.help(docstrings.get("-rpm"));

		
		// docstrings.put("-d", "Maximum number of lines to print for coverage tracks");
		//parser.addArgument("--maxDepthLines", "-d")
		//	.type(Integer.class)
		//	.setDefault(10)
		//	.help("Track height: Maximum number of lines to print for each track");

		docstrings.put("-ml", "Maximum number of lines to print for each methylation track");
		parser.addArgument("--maxMethylLines", "-ml")
			.type(Integer.class)
			.setDefault(10)
			.help(docstrings.get("-ml"));

		parser.addArgument("--maxReadsStack", "-M")
			.type(Integer.class)
			.setDefault(2000)
			.help("Maximum number of reads to accumulate before printing. If more than this many reads map to the window\n"
					+ "randomy select them");
		
		parser.addArgument("--BSseq", "-bs")
			.action(Arguments.storeTrue())
			.help("Bisulphite mode: Mark bases as methylated (M/m) or unmethylated (U/u). Requires -fa");
		
		parser.addArgument("--noFormat", "-nf")
			.action(Arguments.storeTrue())
			.help("Do not format output with non ascii chars (colour, bold, etc.)");

		parser.addArgument("--nonInteractive", "-ni")
			.action(Arguments.storeFalse())
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
