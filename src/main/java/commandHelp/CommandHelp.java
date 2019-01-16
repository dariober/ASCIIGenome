package commandHelp;

import java.util.ArrayList;
import java.util.List;

import com.google.common.base.CharMatcher;
import com.google.common.base.Splitter;

import exceptions.InvalidCommandLineException;

public class CommandHelp {

	private String name;
	private String printName; // Name to use to print on screen or on docs.
	private String args= "";
	private String briefDescription;
	private String additionalDescription= "";
	protected Section inSection= Section.GENERAL;
	private int LINE_LEN= 70;
	
	/* C o n s t r u c t o r */
	
	public CommandHelp(){}
	
	/* M e t h o d s */

	public String printBriefHelp(){
		
		String INDENT= "      ";
		
		String helpStr= this.name + " " + this.args + "\n";
		for(String line : this.wrapLines(this.briefDescription, LINE_LEN - INDENT.length())){
			helpStr += (INDENT + line + "\n");
		}
		return helpStr;
		
	}
	
	public String printCommandHelp(){
	
		String INDENT= "    ";
		
		String helpStr= this.name + " " + this.args + "\n";
		String fullDescr= this.briefDescription + " " + this.additionalDescription;
		for(String line : this.wrapLines(fullDescr, LINE_LEN - INDENT.length())){
			helpStr += (INDENT + line + "\n");
		}
		return helpStr.trim() + "\n";
		
	}
	
	/** Strip some of the reStructuredText formatting from string x.
	 * */
	private String stripReStTextFormat(String x){
		x= x.replaceAll(":: *\\n", ":\n");
		x= x.replaceAll(":code:", "");
		return x;
	}
	
	/** Wrap lines once maxLen is exceeded.
	 * Use ~ to separate words that should stay on the same line and to add spaces.*/
	private List<String> wrapLines(String text, int maxLen){
		
		text= this.stripReStTextFormat(text);
		
		// Words separated by \n only are split into different words since we put a space after \n.
		text= text.replaceAll("\n", "\n "); 
		
		Iterable<String> words = Splitter.on(" ").trimResults(CharMatcher.is(' ')).omitEmptyStrings().split(text);
		
		String line= "";
		List<String> lines= new ArrayList<String>();
		for(String w : words){
			line += (w + " ");
			if(line.trim().length() >= maxLen || w.endsWith("\n")){
				// squiggle is used to mark consecutive spaces that should not be removed. To actually add a squiggle
				// escape it with \~. Here we temprarrily replace \~ with something.
				line= line.trim().replaceAll("\\\\~", "=squiggle=").replaceAll("~", " ").replaceAll("=squiggle=", "~");
				lines.add(line);
				line= "";
			} 
		}
		if(!line.isEmpty()){
			line= line.trim().replaceAll("\\\\~", "=squiggle=").replaceAll("~", " ").replaceAll("=squiggle=", "~");
			lines.add(line);
		}
		return lines;
	}
	
	/* S e t t e r   and   G e t t e r s */
	protected void setName(String name) throws InvalidCommandLineException {
		if(!CommandList.cmds().contains(name)){
			throw new InvalidCommandLineException();
		}
		this.name = name;
	}
	public String getName(){
		return this.name;
	}
	
	protected String getArgs() {
		return this.args;
	}	
	protected void setArgs(String args) {
		this.args = args;
	}
	
	protected void setAdditionalDescription(String description) {
		this.additionalDescription = description;
	}
	protected String getAdditionalDescription() {
		return this.additionalDescription;
	}
	
	protected void setBriefDescription(String briefDescription) {
		this.briefDescription = briefDescription;
	}
	protected String getBriefDescription() {
		return this.briefDescription;
	}

	protected String getPrintName() {
		if(this.printName == null || this.printName.isEmpty()){
			return this.getName();
		}
		return printName;
	}

	protected void setPrintName(String printName) {
		this.printName = printName;
	}


}
