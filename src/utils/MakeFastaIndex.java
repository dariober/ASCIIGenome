package utils;


import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
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
		
		BufferedReader br= new BufferedReader(new FileReader(new File(fastaFile)));
		
		int currLineLen= -1;
		int prevLineLen= -1;
		long byteOffset= 0;
		int seqLen= 0;
		String line;
	
//		while( (br.read()) ){
//			
//		}
		
		return faidx;
	}

	private class fastaIndexRecord {
		String seqname;
		int seqLen= 0;
		int byteOffset= 0;
		int lineLen= 0;
	}
	
}
