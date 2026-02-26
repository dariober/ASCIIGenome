package sortBgzipIndex;

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
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.Arrays;
import org.sqlite.SQLiteException;
import samTextViewer.Utils;
import utils.BedLine;
import utils.BedLineCodec;

public class MakeTabixIndex {

  private File sqliteFile;
  private String columnSeparator;
  private TabixFormat tabixFormat;

  public MakeTabixIndex(String intab, File bgzfOut, TabixFormat tabixFormat)
      throws IOException, InvalidRecordException, ClassNotFoundException, SQLException {
      new MakeTabixIndex(intab, bgzfOut, tabixFormat, Utils.guessSeparator(Path.of(intab)).toString());
  }

  /**
   * Sort, block compress and index the input with format fmt to the given output file. Input is
   * either a local file, possibly compressed, or a URL.
   */
  public MakeTabixIndex(String intab, File bgzfOut, TabixFormat tabixFormat, String columnSeparator)
      throws IOException, InvalidRecordException, ClassNotFoundException, SQLException {

    this.tabixFormat = tabixFormat;
    this.columnSeparator = columnSeparator;

    File tmp = Utils.createTempFile(".asciigenome", "makeTabixIndex.tmp.gz", true);
    File tmpTbi = new File(tmp.getAbsolutePath() + FileExtensions.TABIX_INDEX);
    tmpTbi.deleteOnExit();

    try {
      if (!this.columnSeparator.equals("\t")) {
        throw new InvalidRecordException();
      }
      // Try to block compress and create index assuming the file is sorted
      blockCompressAndIndex(intab, tmp);
    } catch (Exception e) {
      // If intab is not sorted, sort it first.
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

  /**
   * Block compress input file and create associated tabix index.
   *
   * @throws IOException
   * @throws InvalidRecordException
   */
  private void blockCompressAndIndex(String intab, File bgzfOut)
      throws IOException, InvalidRecordException {
    BlockCompressedOutputStream writer = new BlockCompressedOutputStream(bgzfOut);
    long filePosition = writer.getFilePointer();

    TabixIndexCreator indexCreator = new TabixIndexCreator(this.tabixFormat);

    // This is relevant to vcf files only: Prepare header and codec
    // ------------------------------------------------------------
    VCFHeader vcfHeader;
    VCFCodec vcfCodec = null;
    if (this.tabixFormat.equals(TabixFormat.VCF)) {
      try {
        VCFFileReader vcfr = new VCFFileReader(new File(intab), false);
        vcfHeader = vcfr.getFileHeader();
        vcfr.close();
      } catch (MalformedFeatureFile e) {
        vcfHeader = new VCFHeader();
      }
      vcfCodec = new VCFCodec();
      vcfCodec.setVCFHeader(vcfHeader, Utils.getVCFHeaderVersion(vcfHeader));
    }
    // ------------------------------------------------------------

    LineIterator lin = utils.IOUtils.openURIForLineIterator(intab);
    boolean dataLinesFound = false;
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
      addLineToIndex(line, indexCreator, filePosition, vcfCodec);
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
  private void addLineToIndex(
      String line,
      TabixIndexCreator indexCreator,
      long filePosition,
      VCFCodec vcfCodec)
      throws InvalidRecordException {
    if (this.tabixFormat.equals(TabixFormat.BED)) {
      BedLineCodec bedCodec = new BedLineCodec();
      BedLine bed = bedCodec.decode(line);
      indexCreator.addFeature(bed, filePosition);
    } else if (this.tabixFormat.equals(TabixFormat.GFF)) {
      GtfLine gtf = new GtfLine(line.split("\t"));
      indexCreator.addFeature(gtf, filePosition);
    } else if (this.tabixFormat.equals(TabixFormat.VCF)) {
      VariantContext vcf = vcfCodec.decode(line);
      indexCreator.addFeature(vcf, filePosition);
    } else if (this.tabixFormat.flags == TabixFormat.GENERIC_FLAGS) {
      String[] parts = line.split(this.columnSeparator);
      GenericFeature feature = new GenericFeature(parts[tabixFormat.sequenceColumn - 1], tabixFormat.startPositionColumn, tabixFormat.endPositionColumn);
      indexCreator.addFeature(feature, filePosition);
    } else {
      System.err.println(
          "Unexpected TabixFormat: "
              + this.tabixFormat.sequenceColumn
              + " "
              + this.tabixFormat.startPositionColumn);
      throw new InvalidRecordException();
    }
  }

  /**
   * Sort file by columns chrom (text) and pos (int). chromIdx and posIdx are 1-based indexes for
   * the chrom and pos column. For bed use 1 and 2 respectively. For use GTF/GFF 1 and 4. Comment
   * lines, starting with #, are returned as they are. Reading stops if the line ##FASTA is found.
   */
  private void sortByChromThenPos(String unsorted, File sorted)
      throws SQLException, InvalidRecordException, IOException, ClassNotFoundException {

    //int chromIdx = this.tabixFormat.sequenceColumn;
    //int posIdx = this.tabixFormat.startPositionColumn;

    Connection conn = null;
    try {
      this.sqliteFile = Utils.createTempFile(".asciigenome.", ".tmp.sqlite", true);
      conn = this.createSQLiteDb("data");
    } catch (SQLiteException e) {
      this.sqliteFile = File.createTempFile(".asciigenome.", ".tmp.sqlite");
      this.sqliteFile.deleteOnExit();
      conn = this.createSQLiteDb("data");
    }
    PreparedStatement stmtInsert =
        conn.prepareStatement("INSERT INTO data (contig, pos, posEnd, line) VALUES (?, ?, ?, ?)");

    BufferedReader br = Utils.reader(unsorted);
    BufferedWriter wr = new BufferedWriter(new FileWriter(sorted));
    String line;
    int n = 0;
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
//      if ((this.tabixFormat.equals(TabixFormat.BED) || this.tabixFormat.flags == TabixFormat.GENERIC_FLAGS) && !this.columnSeparator.equals("\t")) {
//        line = line.replace(this.columnSeparator, "\t");
//      }
      String[] tabs = line.split(this.columnSeparator);
      if (n == 0 && this.tabixFormat.equals(TabixFormat.BED)) {
        // Allow first uncommented line to fail
        n++;
        try {
          Integer.parseInt(tabs[tabixFormat.startPositionColumn - 1]);
          Integer.parseInt(tabs[tabixFormat.endPositionColumn - 1]);
        } catch (NumberFormatException e) {
          continue;
        }
      }
      stmtInsert.setString(1, tabs[tabixFormat.sequenceColumn - 1]);
      stmtInsert.setInt(2, Integer.parseInt(tabs[tabixFormat.startPositionColumn - 1]));
      if (this.tabixFormat.equals(TabixFormat.VCF)) {
        stmtInsert.setInt(3, 0);
      } else {
        stmtInsert.setInt(3, Integer.parseInt(tabs[tabixFormat.startPositionColumn - 1]));
      }
      stmtInsert.setString(4, line.replaceAll("\n$", ""));
      stmtInsert.executeUpdate();
      dataLinesfound = true;
    }
    stmtInsert.close();
    br.close();

    PreparedStatement stmtSelect =
        conn.prepareStatement("SELECT * FROM data ORDER BY contig, pos, posEnd");

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

  private Connection createSQLiteDb(String tablename) throws SQLException {
    try {
      Class.forName("org.sqlite.JDBC");
    } catch (ClassNotFoundException e) {
      e.printStackTrace();
    }
    Connection conn = DriverManager.getConnection("jdbc:sqlite:" + this.sqliteFile);
    Statement stmt = conn.createStatement();
    String sql =
        "CREATE TABLE "
            + tablename
            + " ("
            + "contig text, "
            + "pos int,"
            + "posEnd int,"
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
