package faidx;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.nio.ByteBuffer;
import java.nio.channels.FileChannel;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.zip.GZIPInputStream;
import java.util.zip.ZipException;

import com.google.common.base.CharMatcher;

public class Faidx {

	/** Read input fasta and write out the corresponding .fai index. See tests/faidx for examples. 
	 * Basically just do:
	 * 
	 * new Faidx(new File("genome.fa"));
	 *  
	 * Index will be genome.fa.fai
	 * */
	public Faidx(File fasta) throws IOException, UnindexableFastaFileException {

		if(this.isCompressed(fasta)){
			// System.err.println(fasta.getAbsolutePath() + " is gzip compressed. Indexing of gzip file is not supported.");
			throw new UnindexableFastaFileException();			
		}
		
		FileChannel fileChannel = FileChannel.open(Paths.get(fasta.getAbsolutePath()));
		int noOfBytesRead = 0;
		StringBuilder sb= new StringBuilder();
		
		boolean isFirstSeqLine= false;
		long currOffset= 0;
		long prevOffset= 0;
		boolean isLast= false; // True when line is expected to be the last one of sequence

		Set<String> seqNames= new HashSet<String>();
		List<FaidxRecord> records= new ArrayList<FaidxRecord>();
		FaidxRecord faidxRecord = null;

		while(noOfBytesRead != -1){
			ByteBuffer buffer = ByteBuffer.allocate(100000);
			noOfBytesRead = fileChannel.read(buffer);
			buffer.flip();
			
			while ( buffer.hasRemaining() ) {
				char x= (char)buffer.get();
				
				if( ! CharMatcher.ASCII.matches(x) ){
					throw new UnindexableFastaFileException();			
				}
				
				currOffset++;
				sb.append(x);
				if(x == '\n'){ // One full line read.
					
					String line= sb.toString();
					sb.setLength(0);
					if(line.trim().isEmpty()){
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
						faidxRecord.byteOffset= currOffset;
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
							faidxRecord.lineFullLength= (int)(currOffset - prevOffset);
							isFirstSeqLine= false;
						} else if(faidxRecord.lineLength != seqLen){
							isLast= true;
						}
					}
					prevOffset= currOffset;

				} // End of processing one full line 
			
			} // End reading chunk of bytes 
			
		} // End of reading channel.

		records.add(faidxRecord); // Add last record
		
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
