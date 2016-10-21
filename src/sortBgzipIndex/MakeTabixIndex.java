package sortBgzipIndex;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.Reader;
import java.net.URL;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.zip.GZIPInputStream;

import org.apache.commons.validator.routines.UrlValidator;
import exceptions.InvalidRecordException;
import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.tribble.index.Index;
import htsjdk.tribble.index.tabix.TabixFormat;
import htsjdk.tribble.index.tabix.TabixIndexCreator;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderVersion;
import utils.BedLine;
import utils.BedLineCodec;
import utils.GtfLine;

public class MakeTabixIndex {

	private File sqliteFile;
	
	/** Sort, block compress and index the input with format fmt to the given output file.
	 * Input is either a local file, possibly compressed, or a URL.
	 * @throws InvalidRecordException 
	 * @throws IOException 
	 * @throws SQLException 
	 * @throws ClassNotFoundException 
	 * */
	public MakeTabixIndex(String intab, File bgzfOut, TabixFormat fmt) throws IOException, InvalidRecordException, ClassNotFoundException, SQLException{

		try{
			// Try to block compress and create index assuming the file is sorted
			blockCompressAndIndex(intab, bgzfOut, fmt);
		} catch(IllegalArgumentException e){
			// If intab is not sorted, sort it first. 
			File sorted= File.createTempFile("asciigenome.", ".sorted.tmp");
			sortByChromThenPos(intab, sorted, fmt);
			blockCompressAndIndex(sorted.getAbsolutePath(), bgzfOut, fmt);
			sorted.delete();
		}
	}

	/**
	 * Block compress input file and create associated tabix index. 
	 * @throws IOException 
	 * @throws InvalidRecordException 
	 * */
	private void blockCompressAndIndex(String intab, File bgzfOut, TabixFormat fmt) throws IOException, InvalidRecordException {
				
		LineIterator lin= utils.IOUtils.openURIForLineIterator(intab);

		BlockCompressedOutputStream writer = new BlockCompressedOutputStream(bgzfOut);
		long filePosition= writer.getFilePointer();
			
		TabixIndexCreator indexCreator=new TabixIndexCreator(fmt);
		
		while(lin.hasNext()){
			
			String line = lin.next().trim();
			if(line.isEmpty() || line.startsWith("#")){
				continue;
			}
			if(line.startsWith("##FASTA")){
				break;
			}			
			
			addLineToIndex(line, indexCreator, filePosition, fmt);
			
			writer.write(line.getBytes());
			writer.write('\n');
			filePosition = writer.getFilePointer();
		}
		writer.flush();
		
		Index index = indexCreator.finalizeIndex(writer.getFilePointer());
		index.writeBasedOnFeatureFile(bgzfOut);
		writer.close();
		
	}

	private void addLineToIndex(String line, TabixIndexCreator indexCreator, long filePosition, TabixFormat fmt) throws InvalidRecordException {

		if(fmt.equals(TabixFormat.BED)){
			
			BedLineCodec bedCodec= new BedLineCodec();
			BedLine bed = bedCodec.decode(line);
			if(bed == null) {
				return;
			}
			indexCreator.addFeature(bed, filePosition);
		
		} else if(fmt.equals(TabixFormat.GFF)){
			GtfLine gtf= new GtfLine(line.split("\t"));
			indexCreator.addFeature(gtf, filePosition);
		
		} else if(fmt.equals(TabixFormat.VCF)) {

			// NB: bgzf file will be headerless!
			VCFCodec vcfCodec= new VCFCodec();
		    VCFHeader header= new VCFHeader();
		    vcfCodec.setVCFHeader(header, VCFHeaderVersion.VCF4_0);
			VariantContext vcf = vcfCodec.decode(line);
			indexCreator.addFeature(vcf, filePosition);
			
		} else {
			System.err.println("Unexpected TabixFormat: " + fmt.sequenceColumn + " " + fmt.startPositionColumn);
			throw new InvalidRecordException();
		}
		
	}

	/** Sort file by columns chrom (text) and pos (int). chromIdx and posIdx are 1-based  
	 * indexes for the chrom and pos column. For bed use 1 and 2 respectively. For use GTF/GFF  1 and 4.
	 * Comment lines, starting with #, are returned as they are. Reading stops if the line ##FASTA is found.
	 * */
	private void sortByChromThenPos(String unsorted, File sorted, TabixFormat fmt) throws SQLException, InvalidRecordException, IOException, ClassNotFoundException{

		int chromIdx= 1;
		int posIdx= 2;
		if(fmt.equals(TabixFormat.BED)){
			//
		} else if(fmt.equals(TabixFormat.GFF)){
			posIdx= 4;
		} else if(fmt.equals(TabixFormat.VCF)){
			posIdx= 2;
		} else {
			System.err.println("Invalid format found");
			throw new InvalidRecordException();
		}
		
		Connection conn= this.createSQLiteDb("data");
		PreparedStatement stmtInsert= conn.prepareStatement("INSERT INTO data (contig, pos, line) VALUES (?, ?, ?)");
		
		BufferedReader br= null;
		InputStream gzipStream= null;
		UrlValidator urlValidator = new UrlValidator();
		if(unsorted.endsWith(".gz")){
			if(urlValidator.isValid(unsorted)) {
				gzipStream = new GZIPInputStream(new URL(unsorted).openStream());
			} else {
				gzipStream = new GZIPInputStream(new FileInputStream(unsorted));
			}
			Reader decoder = new InputStreamReader(gzipStream, "UTF-8");
			br = new BufferedReader(decoder);
		} else if(urlValidator.isValid(unsorted)) {
			InputStream instream= new URL(unsorted).openStream();
			Reader decoder = new InputStreamReader(instream, "UTF-8");
			br = new BufferedReader(decoder);
		} else {
			br = new BufferedReader(new FileReader(unsorted));
		}
		
		BufferedWriter wr= new BufferedWriter(new FileWriter(sorted));
		String line;
		while((line = br.readLine()) != null){
			if(line.trim().startsWith("##FASTA")){
				break;
			}
			if(line.trim().startsWith("#")){
				wr.write(line + "\n");
				continue;
			}
			String[] tabs= line.split("\t"); 
			stmtInsert.setString(1, tabs[chromIdx-1]);
			stmtInsert.setInt(2, Integer.parseInt(tabs[posIdx-1]));
			stmtInsert.setString(3, line.replaceAll("\n$", ""));
			stmtInsert.executeUpdate();
		}
		stmtInsert.close();
		br.close();

		PreparedStatement stmtSelect= conn.prepareStatement("SELECT * FROM data ORDER BY contig, pos");
		
		ResultSet rs= stmtSelect.executeQuery();
		
		while(rs.next()){
			wr.write(rs.getString("line") + "\n");
		}
		conn.commit();
		stmtSelect.close();
		wr.close();
		this.sqliteFile.delete();
	}

	/** Create a tmp sqlite db and return the connection to it. 
	 */
	private Connection createSQLiteDb(String tablename) throws IOException {
		
	    this.sqliteFile= File.createTempFile("asciigenome.", ".tmp.sqlite");
	    this.sqliteFile.deleteOnExit();
	    
	    try {
			Class.forName("org.sqlite.JDBC");
		} catch (ClassNotFoundException e) {
			e.printStackTrace();
		}
	    try {
			Connection conn = DriverManager.getConnection("jdbc:sqlite:" + this.sqliteFile);
	        Statement stmt = conn.createStatement();
	        String sql = "CREATE TABLE " + tablename +  
	        		   " (" +
	                   "contig text, " +
	        		   "pos int," +
	                   "line text" + // This is the row line as read from input file
	                   ")"; 
	        stmt.executeUpdate(sql);
	        stmt= conn.createStatement();

	        // http://stackoverflow.com/questions/1711631/improve-insert-per-second-performance-of-sqlite
	        stmt.execute("PRAGMA journal_mode = MEMORY"); // This is not to leave tmp journal file o disk
			conn.setAutoCommit(false); // This is important: By default each insert is committed 
						     		   // as it is executed, which is slow. Let's commit in bulk at the end instead.
			stmt.close();
	        return conn;
		} catch (SQLException e) {
			e.printStackTrace();
		}
	    return null;
	}

	
}
