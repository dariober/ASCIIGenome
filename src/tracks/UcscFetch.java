package tracks;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.nio.file.Files;
import java.nio.file.attribute.PosixFilePermission;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Calendar;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import org.apache.commons.io.IOUtils;

import com.google.common.base.Joiner;
import com.google.common.base.Splitter;

import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;


import exceptions.InvalidCommandLineException;
import exceptions.InvalidRecordException;
import exceptions.UnableToExecuteUtilException;
import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.tribble.index.Index;
import htsjdk.tribble.index.tabix.TabixFormat;
import htsjdk.tribble.index.tabix.TabixIndexCreator;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.tribble.util.TabixUtils;
import samTextViewer.ArgParse;
import samTextViewer.Main;
import utils.GtfLine;

/** Class to extract annotation from UCSC Genome Browser using genePredToGtf.
 * */
public class UcscFetch {

	private String db;
	private String table;
	private File gtf;
	// This are the dirs in resources where genePredToGtf is to be found. 
	// Executables from http://hgdownload.soe.ucsc.edu/admin/exe/
	static String UTILS_DIR = "/ucscUtils/";
	static final List<String> ARCHIT_DIR= Arrays.asList(new String[] {"linux.x86_64.v287", "macOSX.x86_64", "macOSX.i386", "macOSX.ppc"});
    static File tmpdir= null;
	
	public UcscFetch(String dbNameTableName) throws IOException, InterruptedException, 
		UnableToExecuteUtilException, InvalidCommandLineException, InvalidRecordException {

		tmpdir= Files.createTempDirectory("asciigenome").toAbsolutePath().toFile();
		tmpdir.deleteOnExit();
				
		File genePredToGtfExec= copyExecutableToLocal();
		if(genePredToGtfExec == null){
	        throw new UnableToExecuteUtilException();
		}
		
		checkOrCreateHgConfFile();
		
		File unsortedGtf= executeGenePredToGtf(genePredToGtfExec, dbNameTableName);
		unsortedGtf.deleteOnExit();
		File sortedGtf= new File(unsortedGtf.getAbsoluteFile() + ".sorted.gtf");
		
		try {
			sortGtf(unsortedGtf, sortedGtf);
		} catch (ClassNotFoundException e) {
			e.printStackTrace();
		} catch (SQLException e) {
			e.printStackTrace();
		}
		unsortedGtf.delete();
		
		this.gtf= new File(unsortedGtf.getAbsolutePath() + ".gz");
		this.gtf.deleteOnExit(); // This is kind of redundant as blockCompressAndIndex(true) will delete it anyway  
		blockCompressAndIndex(sortedGtf.getAbsolutePath(), 
				this.gtf.getAbsolutePath(), true);
		sortedGtf.delete();
	}	

	/*      M e t h o d s     */

	/**Copy the resource file src to local file dest.
	 * src is a String like: "/ucscUtils/linux.x86_64.v287/genePredToGtf" 
	 * */
	private void copyResourceFile(String src, String dest) throws IOException{
		
		// NB: You can't use BufferedREader/Writer because this is a binary file not text.
		InputStream res= Main.class.getResourceAsStream(src);
		byte[] bytes = IOUtils.toByteArray(res);
		FileOutputStream fout = new FileOutputStream(dest);
		fout.write(bytes);
		fout.close();
	}
	
	
	/** Copy genePredToGtf executables from jar to a local tmp dir. The version of genePredToGtf
	 * to copy is found by trial and error by executing them one by one until you find the one that returns
	 * the right exit code.
	 * @throws IOException 
	 * @throws InterruptedException 
	 * */
	private File copyExecutableToLocal() throws IOException, InterruptedException{

        Set<PosixFilePermission> perms = new HashSet<PosixFilePermission>();
		perms.add(PosixFilePermission.OWNER_EXECUTE);
		perms.add(PosixFilePermission.OWNER_WRITE);
		perms.add(PosixFilePermission.OWNER_READ);
		perms.add(PosixFilePermission.GROUP_EXECUTE);
		perms.add(PosixFilePermission.GROUP_READ);
		perms.add(PosixFilePermission.GROUP_WRITE);
		perms.add(PosixFilePermission.OTHERS_EXECUTE);
		perms.add(PosixFilePermission.OTHERS_READ);
		perms.add(PosixFilePermission.OTHERS_WRITE);

		for(String x : ARCHIT_DIR){

			File dest= new File(tmpdir, "genePredToGtf");
			dest.deleteOnExit();
			String src= UTILS_DIR + x + "/genePredToGtf";  
			
			this.copyResourceFile(src, dest.getAbsolutePath());
			Files.setPosixFilePermissions(dest.toPath(), perms);
			
			ProcessBuilder probuilder = new ProcessBuilder(dest.getAbsolutePath());
            Process process = probuilder.start();
            process.waitFor();
            if(process.exitValue() == 255 || process.exitValue() == 0){
            	return dest;
            }
		}
		return null;
        
	}
	
	/**  See if file $HOME/.hg.conf exists, if not create it. See also
	 * http://genomewiki.ucsc.edu/index.php/Genes_in_gtf_or_gff_format 
	 * @throws IOException 
	 */
	private static void checkOrCreateHgConfFile() throws IOException{
	
		File hgConf= new File(System.getProperty("user.home"), ".hg.conf");
		
		if(hgConf.isDirectory()) { 
		    throw new RuntimeException(hgConf.getAbsolutePath() + " is a directory!");
		}
		
		if(!hgConf.exists()) { 

			String timeStamp = new SimpleDateFormat("dd/MM/yyyy HH:mm:ss").format(Calendar.getInstance().getTime());
			
			BufferedWriter wr= new BufferedWriter(new FileWriter(hgConf));
			wr.write("# " + timeStamp + " file created by " + ArgParse.PROG_NAME + "\n"
					+ "db.host=genome-mysql.cse.ucsc.edu\n"
					+ "db.user=genomep\n"
					+ "db.password=password\n"
					+ "central.db=hgcentral\n");
			wr.close();
			
			Set<PosixFilePermission> perms = new HashSet<PosixFilePermission>();
			perms.add(PosixFilePermission.OWNER_READ);
			perms.add(PosixFilePermission.OWNER_WRITE);

			Files.setPosixFilePermissions(hgConf.toPath(), perms);
			
		}
	}

	/** Split dbNameTableName, e.g. "hg19:knownGene" into a list as [db, table] or throw exception. 
	 * @throws InvalidCommandLineException  
	 * */
	private List<String> parseCmdArgs(String dbNameTableName) throws InvalidCommandLineException{
		Iterator<String> x = Splitter.on(":").omitEmptyStrings().trimResults().split(dbNameTableName).iterator();
		List<String> args= new ArrayList<String>();
		while(x.hasNext()){
			args.add(x.next());
		}
		if(args.size() != 2){
			System.err.println("Expected a string with database and table name separated by colon, e.g. hg19:knownGene. Got " + dbNameTableName);	
			throw new InvalidCommandLineException();
		}
		return args;
	}
	
	private File executeGenePredToGtf(File genePredToGtfExec, String dbNameTableName) throws InvalidCommandLineException, IOException, InterruptedException {
		
		List<String> args= parseCmdArgs(dbNameTableName);
		this.db= args.get(0);
		this.table= args.get(1);
		
		File outgtf= new File(tmpdir.getAbsolutePath(), this.db + ":" + this.table + ".gtf");
		outgtf.deleteOnExit();
		
    	ProcessBuilder probuilder = new ProcessBuilder(genePredToGtfExec.getAbsolutePath(), "-utr", db, table, outgtf.getAbsolutePath());
    	System.err.println("Querying UCSC...");
    	// System.err.println(Joiner.on(" ").join(probuilder.command()));
    	Process process= probuilder.start();
    	process.waitFor();
    	if(process.exitValue() != 0){
    		//System.err.println("I got problems executing genePredToGtf. Exit code: " + process.exitValue());
    		//System.err.println(Joiner.on(" ").join(probuilder.command()));
    		throw new InvalidCommandLineException();
    	}
    	return outgtf;
    }
	
	/** Sort gtf file using Sqlite
	 * */
	public static void sortGtf(File unsortedGtf, File sortedGtf) throws SQLException, InvalidRecordException, IOException, ClassNotFoundException{

		// http://stackoverflow.com/questions/1711631/improve-insert-per-second-performance-of-sqlite
		GtfSqlite gtfsql= new GtfSqlite();
		gtfsql.getSqlitedb().deleteOnExit();

		Statement stmt= gtfsql.getConn().createStatement();
		stmt.execute("PRAGMA journal_mode = MEMORY"); // This is not to leave tmp journal file o disk
		gtfsql.getConn().setAutoCommit(false); // This is important: By default each insert is committed 
											   // as it is executed, which is slow. Let's commit in bulk at the end instead.

		PreparedStatement stmtInsert= gtfsql.getConn().prepareStatement(gtfsql.getSqlInsert());
		
		BufferedReader br= new BufferedReader(new FileReader(unsortedGtf));
		String gtfLine;
		while((gtfLine = br.readLine()) != null){
			if(gtfLine.trim().startsWith("#")){
				continue;
			}
			if(gtfLine.trim().startsWith("##FASTA")){
				break;
			}
			gtfsql.add(gtfLine, stmtInsert);
		}
		stmtInsert.close();
		br.close();

		BufferedWriter wr= new BufferedWriter(new FileWriter(sortedGtf.getAbsolutePath()));
		
		PreparedStatement stmtSelect= gtfsql.getConn().prepareStatement("SELECT * FROM gtf ORDER BY contig, start, end");
		
		ResultSet rs= stmtSelect.executeQuery();
		
		// ResultSet rs= stmt.executeQuery("SELECT * FROM gtf ORDER BY contig, start, end");
		while(rs.next()){

			StringBuilder sb= new StringBuilder();
			sb.append(rs.getString("contig"));   sb.append("\t");
			sb.append(rs.getString("source"));   sb.append("\t");
			sb.append(rs.getString("feature"));  sb.append("\t");
			sb.append(rs.getInt("start"));       sb.append("\t");
			sb.append(rs.getInt("end"));         sb.append("\t");
			sb.append(rs.getString("score"));    sb.append("\t");
			sb.append(rs.getString("strand"));   sb.append("\t");
			sb.append(rs.getString("frame"));    sb.append("\t");
			sb.append(rs.getString("attribute"));
			wr.write(sb.append("\n").toString());
		}
		gtfsql.getConn().commit();
		stmt.close();
		wr.close();
		gtfsql.getSqlitedb().delete();
	}
	
	/**
	 * Block compress input file and create associated tabix index. Newly created file and index are
	 * deleted on exit if deleteOnExit true.
	 * @throws IOException 
	 * @throws InvalidRecordException 
	 * */
	private void blockCompressAndIndex(String in, String bgzfOut, boolean deleteOnExit) throws IOException, InvalidRecordException {
				
		File inFile= new File(in);
		File outFile= new File(bgzfOut);
		
		LineIterator lin= utils.IOUtils.openURIForLineIterator(inFile.getAbsolutePath());

		BlockCompressedOutputStream writer = new BlockCompressedOutputStream(outFile);
		long filePosition= writer.getFilePointer();
				
		TabixIndexCreator indexCreator=new TabixIndexCreator(TabixFormat.GFF);		
		
		while(lin.hasNext()){
			String line = lin.next();			
			GtfLine gtf= new GtfLine(line.split("\t"));
			writer.write(line.getBytes());
			writer.write('\n');
			indexCreator.addFeature(gtf, filePosition);
			filePosition = writer.getFilePointer();
		}
		writer.flush();
				
		File tbi= new File(bgzfOut + TabixUtils.STANDARD_INDEX_EXTENSION);
		if(tbi.exists() && tbi.isFile()){
			writer.close();
			throw new RuntimeException("Index file exists: " + tbi);
		}
		Index index = indexCreator.finalizeIndex(writer.getFilePointer());
		index.writeBasedOnFeatureFile(outFile);
		writer.close();
		
		if(deleteOnExit){
			outFile.deleteOnExit();
			File idx= new File(outFile.getAbsolutePath() + TabixUtils.STANDARD_INDEX_EXTENSION);
			idx.deleteOnExit();
		}
	}
	
	public File genePredToGtf(){
		return this.gtf;
	}
	
}
