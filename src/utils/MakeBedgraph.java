package utils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.net.MalformedURLException;

import org.apache.commons.math.stat.descriptive.DescriptiveStatistics;

import samTextViewer.Utils;

public class MakeBedgraph {
	
	public MakeBedgraph(String input, int dataCol, File output) throws MalformedURLException, IOException{
		
		BufferedReader br= Utils.reader(input);
		File tmpBedgraph= Utils.createTempFile(".asciigenome.", ".makeBedgraph.bedGraph");
		tmpBedgraph.deleteOnExit();
		BufferedWriter wr= new BufferedWriter(new FileWriter(tmpBedgraph.getAbsoluteFile()));
		String line;
		while((line = br.readLine()) != null){
			if(line.trim().startsWith("##FASTA")){
				break;
			}
			if(line.trim().startsWith("#")){
				wr.write(line + "\n");
				continue;
			}
			String[] arr= line.split("\t");
			String[] bdg= new String[arr.length];
			
			// List<String> lst= new ArrayList<String>();
		}
		br.close();
		wr.close();
	}
}
