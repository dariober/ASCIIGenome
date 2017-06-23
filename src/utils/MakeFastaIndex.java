package utils;


import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.List;

import exceptions.InvalidCommandLineException;

public class MakeFastaIndex {

	final private String fastaFile;
			
	public MakeFastaIndex(String fastaFile) {
		this.fastaFile= fastaFile;
	}
    
	public List<fastaIndexRecord> makeIndex() throws InvalidCommandLineException, FileNotFoundException {
		
		List<fastaIndexRecord> faidx= new ArrayList<fastaIndexRecord>();
		
		return faidx;
	}

	private class fastaIndexRecord {
		String seqname;
		int seqLen= 0;
		int byteOffset= 0;
		int lineLen= 0;
	}
	
}
