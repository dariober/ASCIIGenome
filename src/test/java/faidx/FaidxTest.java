package faidx;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import org.apache.commons.io.FileUtils;
import org.junit.Test;

public class FaidxTest {

  @Test
  public void canIndexValidSeqs() throws IOException, UnindexableFastaFileException {

    File fasta = new File("test_data/faidx/indexable.fa");
    File fai = new File(fasta.getAbsoluteFile() + ".fai");
    fai.deleteOnExit();

    /* Expected index file created with samtools 1.3.1:
       samtools faidx indexable.fa
       mv indexable.fa.fai indexable.fa.fai.expected
    */
    List<String> expected =
        Arrays.asList(
            FileUtils.readFileToString(new File("test_data/faidx/indexable.fa.fai.expected"))
                .split("\n"));

    new Faidx(fasta);
    assertTrue(fai.isFile());

    List<String> observed = Arrays.asList(FileUtils.readFileToString(fai).split("\n"));

    assertEquals(expected.size(), observed.size());

    for (int i = 0; i < expected.size(); i++) {
      assertEquals(expected.get(i), observed.get(i));
    }
  }

  @Test
  public void canIndexGfx() throws IOException, UnindexableFastaFileException {
    File fasta = new File("test_data/faidx/fasta.gff");
    new Faidx(fasta);
    File fai = new File(fasta.getAbsoluteFile() + ".fai");
    fai.deleteOnExit();
    List<String> observed = Arrays.asList(FileUtils.readFileToString(fai).split("\n"));
    assertEquals(2, observed.size());
    assertEquals("1\t1195\t194\t50\t51", observed.get(0));
    assertEquals("2\t1262\t1416\t50\t51", observed.get(1));
  }

  @Test
  public void canHandleWindowsLineEndings() throws IOException, UnindexableFastaFileException {

    File fasta = new File("test_data/faidx/indexable.crlf.fa");
    File fai = new File(fasta.getAbsoluteFile() + ".fai");
    fai.deleteOnExit();

    /* Expected index file created with samtools 1.3.1:
       samtools faidx indexable.crlf.fa
       mv indexable.crlf.fa.fai indexable.crlf.fa.fai.expected
    */
    List<String> expected =
        Arrays.asList(
            FileUtils.readFileToString(new File("test_data/faidx/indexable.crlf.fa.fai.expected"))
                .split("\n"));

    new Faidx(fasta);
    assertTrue(fai.isFile());

    List<String> observed = Arrays.asList(FileUtils.readFileToString(fai).split("\n"));

    assertEquals(expected.size(), observed.size());

    for (int i = 0; i < expected.size(); i++) {
      assertEquals(expected.get(i), observed.get(i));
    }
  }

  @Test
  public void canHandleEmptySequence() throws IOException, UnindexableFastaFileException {
    File fasta = new File("test_data/faidx/empty.fa");
    File fai = new File(fasta.getAbsoluteFile() + ".fai");
    fai.deleteOnExit();

    new Faidx(fasta);

    String[] observed = FileUtils.readFileToString(fai).split("\n");
    for (String x : observed) {
      if (x.startsWith("empty")) { // Empty sequences just check sequence length is 0
        assertEquals("0", x.split("\t")[1]);
      }
    }

    // Can retrieve seqs
    ReferenceSequenceFile faSeqFile =
        ReferenceSequenceFileFactory.getReferenceSequenceFile(new File("test_data/faidx/empty.fa"));
    assertEquals(
        "ACTGNNNNNNNNNNNNN", new String(faSeqFile.getSubsequenceAt("seq", 1, 17).getBases()));
    assertEquals("AACCGGTTNN", new String(faSeqFile.getSubsequenceAt("seq2", 1, 10).getBases()));
    assertEquals(
        "GGGAAATTTNNNCCC", new String(faSeqFile.getSubsequenceAt("seq3", 1, 15).getBases()));

    faSeqFile.close();
  }

  @Test(expected = UnindexableFastaFileException.class)
  public void exceptionOnDuplicateName() throws IOException, UnindexableFastaFileException {
    File fasta = new File("test_data/faidx/dups.fa");
    new Faidx(fasta);
    assertFalse(new File(fasta.getAbsolutePath() + ".fai").exists());
  }

  @Test(expected = UnindexableFastaFileException.class)
  public void exceptionOnDifferentLineLength() throws IOException, UnindexableFastaFileException {
    File fasta = new File("test_data/faidx/lineLen.fa");
    new Faidx(fasta);
    assertFalse(new File(fasta.getAbsolutePath() + ".fai").exists());
  }

  @Test(expected = UnindexableFastaFileException.class)
  public void exceptionOnGzipInput() throws IOException, UnindexableFastaFileException {
    File fasta = new File("test_data/faidx/indexable.fa.gz");
    new Faidx(fasta);
    assertFalse(new File(fasta.getAbsolutePath() + ".fai").exists());
  }

  @Test(expected = UnindexableFastaFileException.class)
  public void exceptionOnNonASCIIchars() throws IOException, UnindexableFastaFileException {
    File fasta = new File("test_data/faidx/nonascii.fa");
    new Faidx(fasta);
    assertFalse(new File(fasta.getAbsolutePath() + ".fai").exists());
  }

  //	@Test
  //	public void testFileChannel() throws IOException{
  //
  //		FileChannel fileChannel = FileChannel.open(Paths.get("chr1.fa"));
  //		long n= 0;
  //		int noOfBytesRead = 0;
  //
  //		long t0= System.nanoTime();
  //
  //		while(noOfBytesRead != -1){
  //			ByteBuffer buffer = ByteBuffer.allocate(10000);
  //			noOfBytesRead = fileChannel.read(buffer);
  //			buffer.flip();
  //			while ( buffer.hasRemaining() ) {
  //				char x= (char)buffer.get();
  //				n++;
  //			}
  //		}
  //		long t1= System.nanoTime();
  //		System.err.println((float)(t1-t0) / 1e6);
  //		System.err.println("nchars: " + n);
  //	}

  //	@Test
  //	public void readArrayTest() throws IOException, InterruptedException {
  //
  //		final int BUFFER_SIZE = 100;
  //		char[] arr = new char[BUFFER_SIZE];
  //
  //	    int result = 0;
  //	    try (BufferedReader reader = new BufferedReader(new FileReader("chr1.fa"))) {
  //	        int charsRead;
  //			long t0= System.nanoTime();
  //	        while ((charsRead = reader.read(arr)) != -1) {
  //	            for (int i = 0; i < charsRead; i++) {
  //	            	result += arr[i];
  //	            }
  //	        }
  //			long t1= System.nanoTime();
  //			System.err.println((t1-t0) / 1000000.0);
  //	    }
  //	    System.err.println(result);
  //	}

  //	@Test
  //	public void testRead() throws IOException, UnindexableFastaFileException{
  //
  //		BufferedReader fa= new BufferedReader(new FileReader(new File("chr1.fa")));
  //		long t0= System.nanoTime();
  //		int c;
  //		while( (c = fa.read()) != -1 ){
  //			if(c == '>'){ System.out.println(c); };
  //		}
  //		long t1= System.nanoTime();
  //		System.err.println(t1-t0); // ~ 7000 ms
  //
  //	}
  //
  //	@Test
  //	public void testReadLine() throws IOException{
  //
  //		BufferedReader fa= new BufferedReader(new FileReader(new File("chr1.fa")));
  //
  //		String line;
  //		long t0= System.nanoTime();
  //		while( (line = fa.readLine()) != null ){
  //			if(line.contains(">")){ System.out.println(line); }
  //		}
  //		long t1= System.nanoTime();
  //		System.err.println(t1-t0); // ~ 700 ms
  //	}

  //	@Test
  //	public void simpleFasta() throws IOException, UnindexableFastaFileException{
  //		File fasta= new File("/Users/berald01/Desktop/asciigenome_demo/genome.fa");
  //		long t0= System.currentTimeMillis();
  //		new Faidx(fasta);
  //		long t1= System.currentTimeMillis();
  //		System.err.println(t1-t0);
  //	}

}
