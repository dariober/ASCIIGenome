package sortBgzipIndex;

import com.google.common.base.Splitter;
import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.FileExtensions;
import htsjdk.tribble.index.Index;
import htsjdk.tribble.index.tabix.TabixFormat;
import htsjdk.tribble.index.tabix.TabixIndexCreator;
import htsjdk.tribble.readers.LineIterator;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.List;
import org.sqlite.SQLiteException;
import samTextViewer.Utils;

public class MakeTabixIndex {

  private File sqliteFile;
  private final char columnSeparator;
  private final TabixFormat tabixFormat;

  public MakeTabixIndex(String intab, File bgzfOut, TabixFormat tabixFormat)
      throws SQLException, IOException, ClassNotFoundException {
    this.columnSeparator = '\t';
    this.tabixFormat = tabixFormat;
    new MakeTabixIndex(intab, bgzfOut, tabixFormat, '\t');
  }

  public MakeTabixIndex(String intab, File bgzfOut, TabixFormat tabixFormat, char columnSeparator)
      throws IOException, ClassNotFoundException, SQLException {

    this.tabixFormat = tabixFormat;
    this.columnSeparator = columnSeparator;

    if (columnSeparator != '\t'
        && (tabixFormat == TabixFormat.BED
            || tabixFormat == TabixFormat.GFF
            || tabixFormat == TabixFormat.SAM
            || tabixFormat == TabixFormat.VCF)) {
      throw new IOException(
          "This tabix format must have tab delimter. Got:'" + this.columnSeparator + "'");
    }

    File tmp = Utils.createTempFile(".asciigenome", "makeTabixIndex.tmp.gz", true);
    File tmpTbi = new File(tmp.getAbsolutePath() + FileExtensions.TABIX_INDEX);
    tmpTbi.deleteOnExit();

    try {
      // Try to block compress and create index assuming the file is sorted
      blockCompressAndIndex(intab, tmp);
    } catch (Exception e) {
      // If intab is not sorted, sort it first.
      PrintWriter pw = new PrintWriter(tmp);
      pw.close();
      File sorted = Utils.createTempFile(".asciigenome.", ".sorted.tmp", true);
      sortByChromThenPos(intab, sorted);
      blockCompressAndIndex(sorted.getAbsolutePath(), tmp);
      Files.delete(Paths.get(sorted.getAbsolutePath()));
    }

    // This renaming and the use of File tmp allows to block compress and index an input file in
    // place.
    // Original intab file is overwritten of course!
    if (bgzfOut.exists()) {
      Files.delete(Paths.get(bgzfOut.getAbsolutePath()));
    }
    Files.move(Paths.get(tmp.getAbsolutePath()), Paths.get(bgzfOut.getAbsolutePath()));

    File bgzfOutTbi = new File(bgzfOut.getAbsolutePath() + FileExtensions.TABIX_INDEX);
    if (bgzfOutTbi.exists()) {
      Files.delete(Paths.get(bgzfOutTbi.getAbsolutePath()));
    }
    Files.move(Paths.get(tmpTbi.getAbsolutePath()), Paths.get(bgzfOutTbi.getAbsolutePath()));
  }

  private void blockCompressAndIndex(String intab, File bgzfOut) throws IOException {
    BlockCompressedOutputStream writer = new BlockCompressedOutputStream(bgzfOut);
    long filePosition = writer.getFilePointer();

    TabixIndexCreator indexCreator = new TabixIndexCreator(this.tabixFormat);

    LineIterator lin = utils.IOUtils.openURIForLineIterator(intab);
    boolean dataLinesFound = false;
    int numHeaderLinesToSkip = tabixFormat.numHeaderLinesToSkip;
    while (lin.hasNext()) {

      String line = lin.next().trim();
      if (line.isEmpty() || line.startsWith("track ")) {
        continue;
      }
      if (line.startsWith("#")) {
        writer.write((line + "\n").getBytes());
        filePosition = writer.getFilePointer();
        continue;
      }
      if (line.startsWith("##FASTA") && dataLinesFound) {
        break;
      }
      if (numHeaderLinesToSkip > 0) {
        writer.write((line + "\n").getBytes());
        filePosition = writer.getFilePointer();
        numHeaderLinesToSkip -= 1;
        continue;
      }
      addLineToIndex(line, indexCreator, filePosition);
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

  /** Set vcfHeader and vcfCodec to null if reading non-vcf line. */
  private void addLineToIndex(String line, TabixIndexCreator indexCreator, long filePosition) {
    List<String> parts = Splitter.on(this.columnSeparator).splitToList(line);
    int start = Integer.parseInt(parts.get(tabixFormat.startPositionColumn - 1));
    GenericFeature feature =
        new GenericFeature(parts.get(tabixFormat.sequenceColumn - 1), start, start + 1);
    indexCreator.addFeature(feature, filePosition);
  }

  private void sortByChromThenPos(String unsorted, File sorted)
      throws SQLException, IOException, ClassNotFoundException {

    Connection conn;
    try {
      this.sqliteFile = Utils.createTempFile(".asciigenome.", ".tmp.sqlite", true);
      conn = this.createSQLiteDb("data");
    } catch (SQLiteException e) {
      this.sqliteFile = File.createTempFile(".asciigenome.", ".tmp.sqlite");
      this.sqliteFile.deleteOnExit();
      conn = this.createSQLiteDb("data");
    }
    PreparedStatement stmtInsert =
        conn.prepareStatement("INSERT INTO data (contig, pos, line) VALUES (?, ?, ?)");

    BufferedReader br = Utils.reader(unsorted);
    BufferedWriter wr = new BufferedWriter(new FileWriter(sorted));
    String line;
    int headerLinesToSkip = this.tabixFormat.numHeaderLinesToSkip;
    boolean dataLinesfound = false;
    while ((line = br.readLine()) != null) {
      if (line.trim().startsWith("##FASTA")
          && dataLinesfound
          && this.tabixFormat.equals(TabixFormat.GFF)) {
        break;
      }
      if (line.trim().startsWith(String.valueOf(this.tabixFormat.metaCharacter))) {
        wr.write(line + "\n");
        continue;
      }
      if (headerLinesToSkip > 0) {
        wr.write(line + "\n");
        headerLinesToSkip -= 1;
        continue;
      }
      if (line.trim().isEmpty()) {
        continue;
      }

      List<String> tabs = Splitter.on(this.columnSeparator).splitToList(line);
      stmtInsert.setString(1, tabs.get(tabixFormat.sequenceColumn - 1));
      stmtInsert.setInt(2, Integer.parseInt(tabs.get(tabixFormat.startPositionColumn - 1)));
      stmtInsert.setString(3, line.replaceAll("\n$", ""));
      stmtInsert.executeUpdate();
      dataLinesfound = true;
    }
    stmtInsert.close();
    br.close();

    PreparedStatement stmtSelect = conn.prepareStatement("SELECT * FROM data ORDER BY contig, pos");

    ResultSet rs = stmtSelect.executeQuery();

    while (rs.next()) {
      wr.write(rs.getString("line") + "\n");
    }
    conn.commit();
    stmtSelect.close();
    wr.close();
    conn.close();
    Files.delete(Paths.get(this.sqliteFile.getAbsolutePath()));
  }

  private Connection createSQLiteDb(String tablename) throws SQLException, ClassNotFoundException {
    Class.forName("org.sqlite.JDBC");
    Connection conn = DriverManager.getConnection("jdbc:sqlite:" + this.sqliteFile);
    Statement stmt = conn.createStatement();
    String sql =
        "CREATE TABLE "
            + tablename
            + " ("
            + "contig text, "
            + "pos int,"
            + "line text"
            + // This is the row line as read from input file
            ")";
    stmt.executeUpdate(sql);
    stmt = conn.createStatement();

    // http://stackoverflow.com/questions/1711631/improve-insert-per-second-performance-of-sqlite
    stmt.execute("PRAGMA journal_mode = OFF"); // This is not to leave tmp journal file o disk
    conn.setAutoCommit(false); // This is important: By default each insert is committed
    // as it is executed, which is slow. Let's commit in bulk at the end instead.
    stmt.close();
    return conn;
  }
}
