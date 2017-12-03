package tracks;

import java.io.File;
import java.io.IOException;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.PreparedStatement;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.List;

import com.google.common.base.Splitter;

import exceptions.InvalidRecordException;
import samTextViewer.Utils;

public class GtfSqlite {

	private File sqlitedb= null; 
	private Connection conn= null;
	private String sqlInsert= "INSERT INTO gtf (" +
            "contig, " +
            "source, " +
            "feature, " + 
            "start, " +
            "end, " +
            "score, " +
            "strand, " +
            "frame, " +
            "attribute" + 
			") VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)";

	/* G E T T E R S   AND   S E T T E R S */
	protected File getSqlitedb() {
		return sqlitedb;
	}

	protected String getSqlInsert() {
		return sqlInsert;
	}

	
	protected Connection getConn() {
		return conn;
	}
	
	/* C O N S T R U C T O R */
	public GtfSqlite() throws IOException, SQLException{
		this.createSQLiteGtfDb();
	}
	
	private void createSQLiteGtfDb() throws IOException {
		
	    // Get a tmp file name: Delete tmp file, keep the name.
	    File tmpfile= Utils.createTempFile(".gtf.", ".tmp.db");
	    tmpfile.deleteOnExit();
	    String sqlitedbName= tmpfile.getName();
	    tmpfile.delete();
	    
	    File sqlitedb= new File(sqlitedbName);

	    try {
			Class.forName("org.sqlite.JDBC");
		} catch (ClassNotFoundException e) {
			e.printStackTrace();
		}
	    try {
			this.conn = DriverManager.getConnection("jdbc:sqlite:" + sqlitedb.getName());
	        Statement stmt = conn.createStatement();
	        String sql = "CREATE TABLE gtf" + 
	        		   " (" +
	                   "contig text, " +
	                   "source text, " +
	                   "feature text, " + 
	                   "start int, " +
	                   "end int, " +
	                   "score text, " +
	                   "strand text, " +
	                   "frame text, " +
	                   "attribute text" + 
	                   ")"; 
	        stmt.executeUpdate(sql);
	        stmt.close();   
	        this.sqlitedb= sqlitedb;
		} catch (SQLException e) {
			e.printStackTrace();
		}

	}
	
	protected void add(String gtfLine, PreparedStatement stmt) throws SQLException, InvalidRecordException {
		// It's important to get the prepared statement from outside so it is not recompiled
		// for each execution of .add().
		
		List<String> gtfRecord= Splitter.on("\t").limit(9).splitToList(gtfLine);
		if(gtfRecord.size() < 8){
			System.err.println("Invalid gtf record. Less than 8 fields found:");
			System.err.println(gtfLine);
			throw new InvalidRecordException();
		}
		if(gtfRecord.size() == 0){
			// Gtf attributes are missing add a dummy one;
			gtfRecord.add("NA");
		}

		stmt.setString(1,  gtfRecord.get(0));
		stmt.setString(2,  gtfRecord.get(1));
		stmt.setString(3,  gtfRecord.get(2));
		stmt.setInt(4,     Integer.parseInt(gtfRecord.get(3)));
		stmt.setInt(5,     Integer.parseInt(gtfRecord.get(4)));
		stmt.setString(6,  gtfRecord.get(5));
		stmt.setString(7,  gtfRecord.get(6));
		stmt.setString(8,  gtfRecord.get(7));
		stmt.setString(9,  gtfRecord.get(8));
		stmt.executeUpdate();
	}
}
