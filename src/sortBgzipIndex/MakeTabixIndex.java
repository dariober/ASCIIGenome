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
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.zip.GZIPInputStream;

import org.apache.commons.validator.routines.UrlValidator;

import com.google.common.base.Joiner;

import exceptions.InvalidRecordException;
import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.tribble.TribbleException.MalformedFeatureFile;
import htsjdk.tribble.index.Index;
import htsjdk.tribble.index.tabix.TabixFormat;
import htsjdk.tribble.index.tabix.TabixIndexCreator;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.tribble.util.TabixUtils;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
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
		
		File tmp = File.createTempFile("asciigenome", "makeTabixIndex.tmp.gz");
		File tmpTbi= new File(tmp.getAbsolutePath() + TabixUtils.STANDARD_INDEX_EXTENSION);
		tmp.deleteOnExit();
		tmpTbi.deleteOnExit();
		
		try{
			// Try to block compress and create index assuming the file is sorted
			blockCompressAndIndex(intab, tmp, fmt);
		} catch(Exception e){
			// If intab is not sorted, sort it first. 
			File sorted= File.createTempFile("asciigenome.", ".sorted.tmp");
			sorted.deleteOnExit();
			sortByChromThenPos(intab, sorted, fmt);
			blockCompressAndIndex(sorted.getAbsolutePath(), tmp, fmt);
			sorted.delete();
		}
		
		// This renaming and the use of File tmp allows to block compress and index an inout file in place.
		// Original intab file is overwritten of course!
		tmp.renameTo(bgzfOut);
		File bgzfOutTbi= new File(bgzfOut.getAbsolutePath() + TabixUtils.STANDARD_INDEX_EXTENSION);
		tmpTbi.renameTo(bgzfOutTbi);
		
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
		
		boolean first= true;
		
		// This is relevant to vcf files only: Prepare header and codec
		// ------------------------------------------------------------
		VCFHeader vcfHeader= null;
		VCFCodec vcfCodec= null;
		if(fmt.equals(TabixFormat.VCF)){
			
			try{
				VCFFileReader vcfr= new VCFFileReader(new File(intab), false);
			    vcfHeader= vcfr.getFileHeader(); // new VCFHeader();
			    vcfr.close();
			} catch(MalformedFeatureFile e){
				vcfHeader= new VCFHeader();
			}
			vcfCodec= new VCFCodec();
		    vcfCodec.setVCFHeader(vcfHeader, this.getVCFHeaderVersion(vcfHeader));
//		    writeVCFHeader(vcfHeader, writer);
//		    filePosition= writer.getFilePointer();
		}
		// ------------------------------------------------------------

		while(lin.hasNext()){
			
			String line = lin.next().trim();
			
			if(line.isEmpty() || line.startsWith("track ")){
				continue;
			}
			if(line.startsWith("#")){
				writer.write((line + "\n").getBytes());
				filePosition = writer.getFilePointer();
				continue;
			}
			if(line.startsWith("##FASTA")){
				break;
			}			
			
			if(first && ! fmt.equals(TabixFormat.VCF)){
				String dummy= this.makeDummyLine(line, fmt);
				addLineToIndex(dummy, indexCreator, filePosition, fmt, null, null);
				
				writer.write(dummy.getBytes());
				writer.write('\n');
				filePosition = writer.getFilePointer();
				first= false;
			}
			addLineToIndex(line, indexCreator, filePosition, fmt, vcfHeader, vcfCodec);
			
			writer.write(line.getBytes());
			writer.write('\n');
			filePosition = writer.getFilePointer();
		}

		writer.flush();
		
		Index index = indexCreator.finalizeIndex(writer.getFilePointer());
		index.writeBasedOnFeatureFile(bgzfOut);
		writer.close();
		CloserUtil.close(lin);
	}

	/** Set vcfHeader and vcfCodec to null if reading non-vcf line.
	 * */
	private void addLineToIndex(String line, TabixIndexCreator indexCreator, long filePosition, TabixFormat fmt, VCFHeader vcfHeader, VCFCodec vcfCodec) throws InvalidRecordException {

		if(fmt.equals(TabixFormat.BED)){
			
			BedLineCodec bedCodec= new BedLineCodec();
			BedLine bed = bedCodec.decode(line);
			indexCreator.addFeature(bed, filePosition);
		
		} else if(fmt.equals(TabixFormat.GFF)){
			GtfLine gtf= new GtfLine(line.split("\t"));
			indexCreator.addFeature(gtf, filePosition);
		
		} else if(fmt.equals(TabixFormat.VCF)) {

			VariantContext vcf = vcfCodec.decode(line);
			indexCreator.addFeature(vcf, filePosition);
			
		} else {
			System.err.println("Unexpected TabixFormat: " + fmt.sequenceColumn + " " + fmt.startPositionColumn);
			throw new InvalidRecordException();
		}	
	}
	
	private VCFHeaderVersion getVCFHeaderVersion(VCFHeader vcfHeader){
		Iterator<VCFHeaderLine> iter = vcfHeader.getMetaDataInInputOrder().iterator();
		while(iter.hasNext()){
			VCFHeaderLine hl = iter.next();
			if(hl.getKey().equals("fileformat")){
				return VCFHeaderVersion.toHeaderVersion(hl.getValue());
			}
		}
		return null;
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
		PreparedStatement stmtInsert= conn.prepareStatement("INSERT INTO data (contig, pos, posEnd, line) VALUES (?, ?, ?, ?)");
		
		BufferedReader br= null;
		InputStream gzipStream= null;
		UrlValidator urlValidator = new UrlValidator();
		if(unsorted.endsWith(".gz") || unsorted.endsWith(".bgz")){
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
			if(fmt.equals(TabixFormat.VCF)){
				stmtInsert.setInt(3, 0);
			} else {
				stmtInsert.setInt(3, Integer.parseInt(tabs[posIdx]));
			}
			stmtInsert.setString(4, line.replaceAll("\n$", ""));
			stmtInsert.executeUpdate();
		}
		stmtInsert.close();
		br.close();

		PreparedStatement stmtSelect= conn.prepareStatement("SELECT * FROM data ORDER BY contig, pos, posEnd");
		
		ResultSet rs= stmtSelect.executeQuery();
		
		while(rs.next()){
			wr.write(rs.getString("line") + "\n");
		}
		conn.commit();
		stmtSelect.close();
		wr.close();
		this.sqliteFile.delete();
	}

	/** Create a dummy line overcome the problem of first line ignored by tabix idnex creator.
	 * This is a horrible hack and it should be fixed in Tabix.
	 * See issue #38
	 * */
	private String makeDummyLine(String line, TabixFormat fmt){
		
		String[] feature= line.split("\t");
		
		List<String> dummy= new ArrayList<String>();
		dummy.add(feature[0]); // chrom
		
		if(fmt.equals(TabixFormat.BED)){
			dummy.add("0");
			dummy.add("1");
			dummy.add("__ignore_me__");
			
		} else if(fmt.equals(TabixFormat.GFF)){
			dummy.add("__ignore_me__");
			dummy.add("__ignore_me__");
			dummy.add("1");
			dummy.add("2");
			dummy.add(".");
			dummy.add(".");
			dummy.add(".");
		
		} else {
			return "";
		}
		
		return Joiner.on("\t").join(dummy);
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
	        		   "posEnd int," +
	                   "line text" + // This is the row line as read from input file
	                   ")"; 
	        stmt.executeUpdate(sql);
	        stmt= conn.createStatement();

	        // http://stackoverflow.com/questions/1711631/improve-insert-per-second-performance-of-sqlite
	        stmt.execute("PRAGMA journal_mode = OFF"); // This is not to leave tmp journal file o disk
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
