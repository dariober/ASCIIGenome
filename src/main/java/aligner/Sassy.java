package aligner;

import static java.nio.file.StandardCopyOption.REPLACE_EXISTING;

import com.google.common.base.Joiner;
import com.google.common.base.Splitter;
import htsjdk.samtools.*;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.SequenceUtil;
import java.io.*;
import java.nio.charset.Charset;
import java.nio.file.*;
import java.util.*;
import java.util.stream.Stream;
import java.util.stream.StreamSupport;
import org.apache.commons.io.FileUtils;
import samTextViewer.GenomicCoords;
import samTextViewer.Utils;

// import static java.nio.file.StandardCopyOption.REPLACE_EXISTING;

public class Sassy {
  private final Path referenceFasta;
  private final Path execPath;
  private final SAMFileHeader samFileHeader;
  private final Path workDir;
  private final int chromOffset;
  private Stream<SAMRecord> samRecords;

  public Sassy(GenomicCoords gc, Path workDir, Path execPath) throws IOException {
    this.workDir = workDir;
    if (execPath != null) {
      this.execPath = execPath;
    } else {
      this.execPath = this.findExecPath();
    }
    this.referenceFasta =
        Paths.get(
            workDir.toString(), gc.getChrom() + "_" + gc.getFrom() + "-" + gc.getTo() + ".fa");
    this.chromOffset = gc.getFrom() - 1;
    this.writeFasta(
        Paths.get(gc.getFastaFile()), gc.getChrom(), gc.getFrom(), gc.getTo(), this.referenceFasta);
    this.samFileHeader = new SAMFileHeader();
    this.samFileHeader.addSequence(gc.getSamSeqDict().getSequence(gc.getChrom()));
  }

  private Path findExecPath() throws IOException {

    try {
      Utils.execSystemCommand(new String[] {}, List.of("sassy", "--version"));
      return Paths.get("sassy");
    } catch (Exception e) {
      //
    }

    String os = System.getProperty("os.name").toLowerCase();
    String sassy;
    if (os.contains("mac")) {
      sassy = "sassy-aarch64-apple-darwin";
    } else if (os.contains("linux")) {
      sassy = "sassy-x86_64-unknown-linux-gnu";
    } else {
      throw new RuntimeException("Unable to find sassy executable");
    }
    try (InputStream instream =
        Sassy.class.getResourceAsStream(Paths.get("/sassy/v0.2.0", sassy).toString())) {
      if (instream == null) {
        throw new RuntimeException("Null resource");
      }
      Path out = Paths.get(this.workDir.toString(), "sassy");
      Files.copy(instream, out, REPLACE_EXISTING);
      out.toFile().setExecutable(true);
      return out;
    }
  }

  public void search(List<String> opts) throws IOException, InterruptedException {
    String sassy = this.getExecPath().toString();
    List<String> cmd = new ArrayList<>();
    cmd.add(sassy);
    cmd.add("search");
    cmd.addAll(opts);
    cmd.add(this.referenceFasta.toString());
    this.replacePatternOptWithPatternFileOpt(cmd);

    SAMProgramRecord pr = new SAMProgramRecord("sassy");
    pr.setCommandLine(Joiner.on(" ").join(cmd));
    this.samFileHeader.addProgramRecord(pr);

    Map<String, byte[]> querySequence;
    if (cmd.contains("-f") || cmd.contains("--pattern-fasta")) {
      int i = cmd.contains("-f") ? cmd.indexOf("-f") + 1 : cmd.indexOf("--pattern-fasta") + 1;
      querySequence = this.fastaFileToMap(Path.of(cmd.get(i)));
    } else if (cmd.contains("-l") || cmd.contains("--pattern-file")) {
      int i = cmd.contains("-l") ? cmd.indexOf("-l") + 1 : cmd.indexOf("--pattern-file") + 1;
      querySequence = this.patternFileToMap(Path.of(cmd.get(i)));
    } else {
      throw new RuntimeException(
          "Please provide one of:\n"
              + "  -p [comma separated patterns]\n"
              + "  -l [file with one pattern per line]\n"
              + "  -f [fasta file of patterns]");
    }

    Stream<String> lines = Utils.execSystemCommandStream(Stream.empty(), cmd);
    Iterator<String> iter = lines.iterator();

    if (!iter.hasNext()) {
      lines.close();
      throw new RuntimeException("No output from command");
    }

    String colnames = iter.next();
    if (!colnames.contains("cigar")) {
      lines.close();
      throw new RuntimeException("Unexpected column names: " + colnames);
    }

    this.samRecords =
        StreamSupport.stream(Spliterators.spliteratorUnknownSize(iter, Spliterator.ORDERED), false)
            .map(l -> searchToSAMRecord(l, querySequence))
            .onClose(lines::close);
  }

  /** Replace -p/--pattern option with --pattern-file */
  private void replacePatternOptWithPatternFileOpt(List<String> opts) throws IOException {
    if (opts.contains("-p") || opts.contains("--pattern")) {
      // Replace the -p/--pattern option with --pattern-file
      int pidx = opts.contains("-p") ? opts.indexOf("-p") : opts.indexOf("--pattern");
      opts.set(pidx, "--pattern-file");

      // Get patterns
      String patterns = opts.get(pidx + 1).replaceAll(",", "\n");
      // Write patterns to temp file:
      Path patternFile = Paths.get(this.workDir.toString(), "query.txt");
      FileUtils.write(patternFile.toFile(), patterns + "\n", Charset.defaultCharset());
      // Replace pattern arg with pattern file
      opts.set(pidx + 1, patternFile.toString());
    }
  }

  public long writeSAMFile(Path samfile) {
    try (SAMFileWriter writer =
        new SAMFileWriterFactory().makeSAMWriter(this.samFileHeader, false, samfile)) {
      return this.samRecords.peek(writer::addAlignment).count();
    }
  }

  private void writeFasta(Path inFasta, String chrom, int from, int to, Path outFasta)
      throws IOException {
    int lineWidth = 60;
    int chunkSize = lineWidth * 50;
    try (IndexedFastaSequenceFile fastaReader = new IndexedFastaSequenceFile(inFasta);
        BufferedWriter writer = Files.newBufferedWriter(outFasta)) {

      writer.write(">" + chrom);
      writer.newLine();

      int pos = from;

      while (pos <= to) {
        int end = Math.min(pos + chunkSize - 1, to);

        ReferenceSequence chunkSeq = fastaReader.getSubsequenceAt(chrom, pos, end);
        byte[] bases = chunkSeq.getBases();

        for (int i = 0; i < bases.length; i += lineWidth) {
          int lineEnd = Math.min(i + lineWidth, bases.length);
          writer.write(new String(bases, i, lineEnd - i));
          writer.newLine();
        }

        pos += chunkSize;
      }
    }
  }

  private Map<String, byte[]> patternFileToMap(Path patternFile) {
    Map<String, byte[]> map = new HashMap<>();
    try (BufferedReader br = new BufferedReader(new FileReader(String.valueOf(patternFile)))) {
      String line;
      int n = 1;
      while ((line = br.readLine()) != null) {
        map.put(String.valueOf(n), line.trim().getBytes());
        n++;
      }
    } catch (IOException e) {
      throw new RuntimeException(e);
    }
    return map;
  }

  private Map<String, byte[]> fastaFileToMap(Path fasta) throws IOException {

    Map<String, byte[]> map = new HashMap<>();

    try (ReferenceSequenceFile fastaReader =
        ReferenceSequenceFileFactory.getReferenceSequenceFile(fasta)) {

      ReferenceSequence seq;

      while ((seq = fastaReader.nextSequence()) != null) {
        String header = seq.getName();
        String name = header.split("\\s+", 2)[0];

        if (map.containsKey(name)) {
          throw new IllegalArgumentException("Duplicate FASTA sequence name: " + name);
        }
        map.put(name, seq.getBases());
      }
    }

    return map;
  }

  private SAMRecord searchToSAMRecord(String line, Map<String, byte[]> querySequence) {
    List<String> lst = Splitter.on("\t").splitToList(line.trim());

    String patternName = lst.get(0).split("\\s+", 2)[0];

    SAMRecord rec = new SAMRecord(this.samFileHeader);

    rec.setReadName(patternName);
    rec.setReferenceName(lst.get(1));
    rec.setAttribute("NM", Integer.parseInt(lst.get(2)));
    rec.setReadNegativeStrandFlag(lst.get(3).equals("-"));
    rec.setAlignmentStart(Integer.parseInt(lst.get(4)) + 1 + this.chromOffset);
    if (rec.getReadNegativeStrandFlag()) {
      Cigar cigar = TextCigarCodec.decode(lst.get(7));
      List<CigarElement> elements = new ArrayList<>(cigar.getCigarElements());
      Collections.reverse(elements);
      rec.setCigar(new Cigar(elements));
    } else {
      rec.setCigarString(lst.get(7));
    }
    rec.setMappingQuality(255);

    byte[] readBytes = querySequence.get(patternName).clone();
    if (rec.getReadNegativeStrandFlag()) {
      SequenceUtil.reverseComplement(readBytes);
    }
    rec.setReadBases(readBytes);

    return rec;
  }

  public Path getWorkDir() {
    return this.workDir;
  }

  public Path getExecPath() {
    return this.execPath;
  }

  public SAMFileHeader getSamFileHeader() {
    return this.samFileHeader;
  }

  public Stream<SAMRecord> getSamRecords() {
    return this.samRecords;
  }
}
