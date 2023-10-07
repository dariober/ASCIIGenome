package sortBgzipIndex;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;

import org.sqlite.SQLiteException;


import exceptions.InvalidRecordException;
import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.FileExtensions;
import htsjdk.tribble.TribbleException.MalformedFeatureFile;
import htsjdk.tribble.index.Index;
import htsjdk.tribble.index.tabix.TabixFormat;
import htsjdk.tribble.index.tabix.TabixIndexCreator;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import samTextViewer.Utils;
import utils.BedLine;
import utils.BedLineCodec;
import utils.GtfLine;

public class MakeTabixIndex {

    private File sqliteFile;
    final private String columnSeparator;
    final private TabixFormat tabixFormat;
    
    /** Sort, block compress and index the input with format fmt to the given output file.
     * Input is either a local file, possibly compressed, or a URL.
     * @throws InvalidRecordException 
     * @throws IOException 
     * @throws SQLException 
     * @throws ClassNotFoundException 
     * */
    public MakeTabixIndex(String intab, File bgzfOut, TabixFormat tabixFormat) throws IOException, InvalidRecordException, ClassNotFoundException, SQLException{
        
        this.tabixFormat = tabixFormat;
        this.columnSeparator = Utils.getColumnSeparator(intab);
        
        File tmp = Utils.createTempFile(".asciigenome", "makeTabixIndex.tmp.gz", true);
        File tmpTbi= new File(tmp.getAbsolutePath() + FileExtensions.TABIX_INDEX);
        tmpTbi.deleteOnExit();

        try{
            if(!this.columnSeparator.equals("\t")) {
                throw new InvalidRecordException();
            }
            // Try to block compress and create index assuming the file is sorted
            blockCompressAndIndex(intab, tmp);
        } catch(Exception e){
            // If intab is not sorted, sort it first.
            File sorted= Utils.createTempFile(".asciigenome.", ".sorted.tmp", true);
            sortByChromThenPos(intab, sorted);
            blockCompressAndIndex(sorted.getAbsolutePath(), tmp);
            Files.delete(Paths.get(sorted.getAbsolutePath()));
        }
        
        // This renaming and the use of File tmp allows to block compress and index an input file in place.
        // Original intab file is overwritten of course!
        if(bgzfOut.exists()){
            Files.delete(Paths.get(bgzfOut.getAbsolutePath()));
        }
        Files.move(Paths.get(tmp.getAbsolutePath()), Paths.get(bgzfOut.getAbsolutePath()));

        File bgzfOutTbi= new File(bgzfOut.getAbsolutePath() + FileExtensions.TABIX_INDEX);
        if(bgzfOutTbi.exists()){
            Files.delete(Paths.get(bgzfOutTbi.getAbsolutePath()));
        }
        Files.move(Paths.get(tmpTbi.getAbsolutePath()), Paths.get(bgzfOutTbi.getAbsolutePath()));
        
    }

    /**
     * Block compress input file and create associated tabix index. 
     * @throws IOException 
     * @throws InvalidRecordException 
     * */
    private void blockCompressAndIndex(String intab, File bgzfOut) throws IOException, InvalidRecordException {
        BlockCompressedOutputStream writer = new BlockCompressedOutputStream(bgzfOut);
        long filePosition= writer.getFilePointer();
            
        TabixIndexCreator indexCreator=new TabixIndexCreator(this.tabixFormat);
        
        // This is relevant to vcf files only: Prepare header and codec
        // ------------------------------------------------------------
        VCFHeader vcfHeader= null;
        VCFCodec vcfCodec= null;
        if(this.tabixFormat.equals(TabixFormat.VCF)){
            try{
                VCFFileReader vcfr= new VCFFileReader(new File(intab), false);
                vcfHeader= vcfr.getFileHeader();
                vcfr.close();
            } catch(MalformedFeatureFile e){
                vcfHeader= new VCFHeader();
            }
            vcfCodec= new VCFCodec();
            vcfCodec.setVCFHeader(vcfHeader, Utils.getVCFHeaderVersion(vcfHeader));
        }
        // ------------------------------------------------------------

        int nWarnings= 10;
        LineIterator lin= utils.IOUtils.openURIForLineIterator(intab);
        boolean dataLinesFound = false;
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
            if(line.startsWith("##FASTA") && dataLinesFound){
                break;
            }			
            try {
                addLineToIndex(line, indexCreator, filePosition, vcfHeader, vcfCodec);
            } catch(Exception e){
                if(e.getMessage().contains("added out sequence of order") || 
                        e.getMessage().contains("Features added out of order")){
                    // Get a string marker for out-of-order from htsjdk/tribble/index/tabix/TabixIndexCreator.java 
                    writer.close();
                    throw new InvalidRecordException();
                }
                if(nWarnings >= 0){
                    System.err.println("Warning: " + e.getMessage() + ". Skipping:\n" + line);
                }
                if(nWarnings == 0){
                    System.err.println("Additional warnings will not be show.");
                }
                nWarnings--;
            }
            writer.write(line.getBytes());
            writer.write('\n');
            filePosition = writer.getFilePointer();
            dataLinesFound = true;
        }

        writer.flush();
        
        Index index = indexCreator.finalizeIndex(writer.getFilePointer());
        index.writeBasedOnFeatureFile(bgzfOut);
        writer.close();
        CloserUtil.close(lin);
    }

    /** Set vcfHeader and vcfCodec to null if reading non-vcf line.
     * */
    private void addLineToIndex(String line, TabixIndexCreator indexCreator, long filePosition, VCFHeader vcfHeader, VCFCodec vcfCodec) throws InvalidRecordException {
        if(this.tabixFormat.equals(TabixFormat.BED)){
            BedLineCodec bedCodec= new BedLineCodec();
            BedLine bed = bedCodec.decode(line);
            indexCreator.addFeature(bed, filePosition);
        } else if(this.tabixFormat.equals(TabixFormat.GFF)){
            GtfLine gtf= new GtfLine(line.split("\t"));
            indexCreator.addFeature(gtf, filePosition);
        } else if(this.tabixFormat.equals(TabixFormat.VCF)) {
            VariantContext vcf = vcfCodec.decode(line);
            indexCreator.addFeature(vcf, filePosition);
        } else {
            System.err.println("Unexpected TabixFormat: " + this.tabixFormat.sequenceColumn + " " + this.tabixFormat.startPositionColumn);
            throw new InvalidRecordException();
        }	
    }
    
    /** Sort file by columns chrom (text) and pos (int). chromIdx and posIdx are 1-based  
     * indexes for the chrom and pos column. For bed use 1 and 2 respectively. For use GTF/GFF  1 and 4.
     * Comment lines, starting with #, are returned as they are. Reading stops if the line ##FASTA is found.
     * */
    private void sortByChromThenPos(String unsorted, File sorted) throws SQLException, InvalidRecordException, IOException, ClassNotFoundException{

        int chromIdx= 1;
        int posIdx= 2;
        if(this.tabixFormat.equals(TabixFormat.BED)){
            //
        } else if(this.tabixFormat.equals(TabixFormat.GFF)){
            posIdx= 4;
        } else if(this.tabixFormat.equals(TabixFormat.VCF)){
            posIdx= 2;
        } else {
            System.err.println("Invalid format found");
            throw new InvalidRecordException();
        }
        
        Connection conn= null;
        try{
            this.sqliteFile= Utils.createTempFile(".asciigenome.", ".tmp.sqlite", true);
            conn= this.createSQLiteDb(this.sqliteFile, "data");
        } catch(SQLiteException e){
            this.sqliteFile= File.createTempFile(".asciigenome.", ".tmp.sqlite");
            this.sqliteFile.deleteOnExit();
            conn= this.createSQLiteDb(this.sqliteFile, "data");
        }
        PreparedStatement stmtInsert= conn.prepareStatement("INSERT INTO data (contig, pos, posEnd, line) VALUES (?, ?, ?, ?)");
    
        BufferedReader br= Utils.reader(unsorted);
        BufferedWriter wr= new BufferedWriter(new FileWriter(sorted));
        String line;
        int n = 0;
        boolean dataLinesfound = false;
        while((line = br.readLine()) != null){
            if(line.trim().startsWith("##FASTA") && dataLinesfound && this.tabixFormat.equals(TabixFormat.GFF)){
                break;
            }
            if(line.trim().startsWith("#")){
                wr.write(line + "\n");
                continue;
            }
            if(line.trim().isEmpty()){
                continue;
            }
            if(this.tabixFormat.equals(TabixFormat.BED) && !this.columnSeparator.equals("\t")) {
                line = line.replace(this.columnSeparator, "\t");
            }
            String[] tabs = line.split("\t");
            if(n == 0 && this.tabixFormat.equals(TabixFormat.BED)) {
                // Allow first uncommented line to fail
                n++;
                try {
                    Integer.parseInt(tabs[1]);
                    Integer.parseInt(tabs[2]);
                } catch(NumberFormatException e) {
                    continue;
                }
            }
            stmtInsert.setString(1, tabs[chromIdx-1]);
            stmtInsert.setInt(2, Integer.parseInt(tabs[posIdx-1]));
            if(this.tabixFormat.equals(TabixFormat.VCF)){
                stmtInsert.setInt(3, 0);
            } else {
                stmtInsert.setInt(3, Integer.parseInt(tabs[posIdx]));
            }
            stmtInsert.setString(4, line.replaceAll("\n$", ""));
            stmtInsert.executeUpdate();
            dataLinesfound = true;
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
        conn.close();
        Files.delete(Paths.get(this.sqliteFile.getAbsolutePath()));
        
    }

    /** Create a tmp sqlite db and return the connection to it. 
     * @throws SQLException 
     */
    private Connection createSQLiteDb(File sqliteFile, String tablename) throws SQLException {
        // this.sqliteFile= Utils.createTempFile(".asciigenome.", ".tmp.sqlite");
        
        try {
            Class.forName("org.sqlite.JDBC");
        } catch (ClassNotFoundException e) {
            e.printStackTrace();
        }
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
    }
}
