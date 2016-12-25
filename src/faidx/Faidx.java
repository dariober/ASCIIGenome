package faidx;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.RandomAccessFile;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.zip.GZIPInputStream;
import java.util.zip.ZipException;

public class Faidx {

	public Faidx(File fasta) throws IOException, UnindexableFastaFileException {

		if(this.isCompressed(fasta)){
			System.err.println(fasta.getAbsolutePath() + " is gzip compressed. Indexing of gzip file is not supported.");
			throw new UnindexableFastaFileException();			
		}
		
		RandomAccessFile fa = new RandomAccessFile(fasta, "r");

		Set<String> seqNames= new HashSet<String>();
		List<FaidxRecord> records= new ArrayList<FaidxRecord>();
		try{
			FaidxRecord faidxRecord = null;
			String line= "";
			boolean isFirstSeqLine= false;
			long prevOffset= 0;
			boolean isLast= false; // True when line is expected to be the last one of sequence
			while( (line= fa.readLine()) != null ){
				if(line.trim().isEmpty()){ // FIXME: Empty lines should be allowed if they occur at the end of the file 
					                       // or before a new sequence name.
					isLast= true;
					continue;
				}
				if(line.startsWith(">")){
					isLast= false;
					if(faidxRecord != null){
						records.add(faidxRecord);
					}
					faidxRecord= new FaidxRecord();
					faidxRecord.makeSeqNameFromRawLine(line);
					
					if(seqNames.contains(faidxRecord.getSeqName())){
						System.err.println(fasta.getAbsolutePath() + ": Duplicate sequence name found for " + faidxRecord.getSeqName());
						throw new UnindexableFastaFileException();
					} else {
						seqNames.add(faidxRecord.getSeqName());
					}
					faidxRecord.byteOffset= fa.getFilePointer();
					isFirstSeqLine= true;
				} else {
					if(isLast){
						System.err.println(fasta.getAbsolutePath() + ": Different line length in " + faidxRecord.getSeqName());
						throw new UnindexableFastaFileException();
					}
					int seqLen= line.replaceAll("\\s", "").length();
					faidxRecord.seqLength += seqLen;
					if(isFirstSeqLine){
						faidxRecord.lineLength= seqLen;
						faidxRecord.lineFullLength= (int)(fa.getFilePointer() - prevOffset);
						isFirstSeqLine= false;
					} else if(faidxRecord.lineLength != seqLen){
						isLast= true;
					}
				}
				prevOffset= fa.getFilePointer();
			}
			records.add(faidxRecord);
		} finally {
			fa.close();
		}

		// Write out index
		BufferedWriter wr= new BufferedWriter(new FileWriter(new File(fasta.getAbsolutePath() + ".fai")));
		for(FaidxRecord rec : records){
			wr.write(rec.toString() + "\n");	
		}
		wr.close();
	}
	
	private boolean isCompressed(File fasta) throws IOException{
		try{
		InputStream fileStream = new FileInputStream(fasta);
		InputStream gzipStream = new GZIPInputStream(fileStream);
		gzipStream.close();
		} catch(ZipException e){
			return false;
		}
		return true;
	}
	
}
