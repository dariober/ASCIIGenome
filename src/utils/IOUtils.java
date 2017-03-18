/*
The MIT License (MIT)

Copyright (c) 2014 Pierre Lindenbaum

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.


History:
* 2014 creation

*/
package utils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.ByteArrayInputStream;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FilterOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.io.PushbackInputStream;
import java.io.Reader;
import java.io.Writer;
import java.net.URL;
import java.nio.charset.Charset;
import java.util.LinkedHashSet;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import htsjdk.samtools.Defaults;
import htsjdk.samtools.util.AbstractIterator;
import htsjdk.samtools.util.BlockCompressedInputStream;
import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.tribble.readers.LineIteratorImpl;
import htsjdk.tribble.readers.LineReader;
import htsjdk.tribble.readers.LineReaderUtil;

public class IOUtils {
	
	public static void copyTo(final File f,final OutputStream fous) throws IOException
		{
		FileInputStream fin=null;
		try {
			fin =new FileInputStream(f);
			copyTo(fin,fous);
			fous.flush();
		} finally {
			CloserUtil.close(fin);
			}
		}
	
	public static void copyTo(final File f,final Writer fous) throws IOException
		{
		FileReader fin=null;
		try {
			fin =new FileReader(f);
			copyTo(fin,fous);
			fous.flush();
		} finally {
			CloserUtil.close(fin);
			}
		}
	
	public static void copyTo(InputStream in,File f) throws IOException
		{
		FileOutputStream fous=new FileOutputStream(f);
		copyTo(in,fous);
		fous.flush();
		fous.close();
		}
	public static void copyTo(InputStream in,OutputStream out) throws IOException
		{
		byte buffer[]=new byte[2048];
		int nRead;
		while((nRead=in.read(buffer))!=-1)
			{
			out.write(buffer,0,nRead);
			}
		out.flush();
		}

	public static void copyTo(Reader in,Writer out) throws IOException
		{
		char buffer[]=new char[2048];
		int nRead;
		while((nRead=in.read(buffer))!=-1)
			{
			out.write(buffer,0,nRead);
			}
		out.flush();
		}
	
	public static boolean isRemoteURI(String uri)
		{	
	    return  uri.startsWith("http://") ||
				uri.startsWith("https://") ||
				uri.startsWith("ftp://")
				;
		}
	private static InputStream tryBGZIP(InputStream in) throws IOException
		{
		byte buffer[]=new byte[2048];
		
		PushbackInputStream push_back=new PushbackInputStream(in,buffer.length+10);
		int nReads=push_back.read(buffer);
		push_back.unread(buffer, 0, nReads);
		
		try
			{
			BlockCompressedInputStream bgz=new BlockCompressedInputStream(new ByteArrayInputStream(buffer, 0, nReads));
			bgz.read();
			bgz.close();
			return new BlockCompressedInputStream(push_back);
			}
		catch(Exception err)
			{
			return new GZIPInputStream(push_back);
			}
		}
	
	public static InputStream openURIForReading(String uri) throws IOException
		{
		if(isRemoteURI(uri))
			{
			URL url=new URL(uri);
			InputStream in=url.openStream();
			//do we have .... azdazpdoazkd.vcf.gz?param=1&param=2
			int question=uri.indexOf('?');
			if(question!=-1) uri=uri.substring(0, question);
			if(uri.endsWith(".gz") || uri.endsWith(".bgz"))
				{
				return tryBGZIP(in);
				}
			return in;
			}
		if(uri.startsWith("file://"))
			{
			uri=uri.substring(7);
			}
		return openFileForReading(new File(uri));
		}
	
	public static BufferedReader openURIForBufferedReading(String uri) throws IOException
		{
		return  new BufferedReader(new InputStreamReader(openURIForReading(uri), Charset.forName("UTF-8")));
		}


	
	@SuppressWarnings("resource")
	public static InputStream openFileForReading(File file) throws IOException
		{
		IOUtil.assertFileIsReadable(file);
		InputStream in= new FileInputStream(file);
		if(file.getName().endsWith(".gz") || file.getName().endsWith(".bgz"))
			{
			in=tryBGZIP(in);
			}
		return in;
		}
	
	public static BufferedReader openFileForBufferedReading(File file) throws IOException
		{
		return  new BufferedReader(new InputStreamReader(openFileForReading(file), Charset.forName("UTF-8")));
		}

    public static BufferedWriter openFileForBufferedWriting(final File file)  throws IOException
    	{
        return new BufferedWriter(new OutputStreamWriter(openFileForWriting(file)), Defaults.BUFFER_SIZE);
    	}
    
    public static OutputStream openFileForWriting(final File file) throws IOException
    	{
        if (file.getName().endsWith(".vcf.gz"))
        	{
            return new BlockCompressedOutputStream(file);
        	}
        else if (file.getName().endsWith(".gz"))
	    	{
	        return new GZIPOutputStream(new FileOutputStream(file));
	    	}
        else
        	{
            return new FileOutputStream(file);
        	}         
    	}
    
    /** open a printwriter, compress if it ends with *.gz  */
    public static PrintWriter openFileForPrintWriter(final File file) throws IOException
		{
	    if (file.getName().endsWith(".gz"))
	    	{
	        return new PrintWriter(openFileForWriting(file));
	    	}
	    else
	    	{
	        return new PrintWriter(file);
	    	}         
		}

    
    public static LineReader openFileForLineReader(File file) throws IOException
		{
		return  LineReaderUtil.fromBufferedStream(openFileForReading(file));
		}
    /** @return a LineIterator that should be closed with CloserUtils */
    public static LineIterator openFileForLineIterator(File file) throws IOException
  		{
  		return  new LineIteratorImpl(openFileForLineReader(file));
  		}
    
    public static LineReader openStdinForLineReader() throws IOException
		{
		return  LineReaderUtil.fromBufferedStream(System.in);
		}
    
    public static BufferedReader openStdinForBufferedReader() throws IOException
		{
		return openStreamForBufferedReader(System.in);
		}
    public static BufferedReader openStreamForBufferedReader(InputStream in) throws IOException
		{
		return  new BufferedReader(new InputStreamReader(in));
		}

    /** @return a LineIterator that should be closed with CloserUtils */
    public static LineIterator openStreamForLineIterator(final InputStream in) throws IOException
  		{
  		return  new LineIteratorImpl(openStreamForLineReader(in));
  		}
    
    /** @return a LineIterator that should be closed with CloserUtils */
    public static LineIterator openStdinForLineIterator() throws IOException
  		{
  		return  openStreamForLineIterator(System.in);
  		}

    public static LineReader openStreamForLineReader(final InputStream in) throws IOException
		{
		return  LineReaderUtil.fromBufferedStream(in);
		}

    
    public static LineReader openURIForLineReader(String uri) throws IOException
		{
		return  LineReaderUtil.fromBufferedStream(openURIForReading(uri));
		}
    /** @return a LineIterator that should be closed with CloserUtils */
    public static LineIterator openURIForLineIterator(String uri) throws IOException
  		{
  		return  new LineIteratorImpl(openURIForLineReader(uri));
  		}
    
    /** read String from DataInputStream
     * motivation: readUTF can't print lines larger than USHORTMAX
     *  */
    public static String readString(DataInputStream in) throws IOException
    	{
    	int llength=in.readInt();
    	if(llength==-1) return null;
    	byte a[]=new byte[llength];
    	in.readFully(a, 0, llength);
    	return new String(a);
    	}
    
    /** write String to DataOutputStream
     * motivation: DataInputStream.readUTF can't print lines larger than USHORTMAX
     *  */
    public static void writeString(DataOutputStream os,String s) throws IOException
		{
		if(s==null)
			{
			os.writeInt(-1);
			}
		else
			{
			byte array[]=s.getBytes();
			os.writeInt(array.length);
			os.write(array);
			}
		}
    
    /** return 'inputs' as set, if a filename ends with '*.list'
     * is is considered as a file:one file per line
     * @param inputs files
     * @return set of files
     */
	public static LinkedHashSet<String> unrollFiles(java.util.Collection<String> inputs)
		{
		LinkedHashSet<String> vcfFiles= new LinkedHashSet<>(inputs.size()+1);
		for(String file : inputs)
			{
			if(file.isEmpty())
				{
				//ignore
				}
			else if(IOUtil.isUrl(file))
				{
				vcfFiles.add(file);
				}
			else if(file.endsWith(".list"))
				{
				File f = new File(file);
				IOUtil.assertFileIsReadable(f);
				 for (final String s : IOUtil.readLines(f))
				 	{
					if (s.endsWith("#")) continue;
					if (s.trim().isEmpty()) continue;
					vcfFiles.add(s);
					}
				}
			else
				{
				vcfFiles.add(file);
				}
			
			}
		return vcfFiles;
		}
	
    /** return 'inputs' as set, if a filename ends with '*.list'
     * is is considered as a file:one file per line
     * @param inputs files
     * @return set of files
     */
	public static LinkedHashSet<File> unrollFileCollection(java.util.Collection<File> inputs)
		{
		LinkedHashSet<File> vcfFiles= new LinkedHashSet<>(inputs.size()+1);
		for(File file : inputs)
			{
			if(file.getName().endsWith(".list"))
				{
				IOUtil.assertFileIsReadable(file);
				 for (final String s : IOUtil.readLines(file))
				 	{
					if (s.endsWith("#")) continue;
					if (s.trim().isEmpty()) continue;
					vcfFiles.add(new File(s));
					}
				}
			else
				{
				vcfFiles.add(file);
				}
			
			}
		return vcfFiles;
		}

	
	/** converts a BufferedReader to a line Iterator */
	public static  LineIterator toLineIterator(BufferedReader r)
		{
		return new BuffReadIter(r);
		}
	private static class BuffReadIter 
		extends AbstractIterator<String>
		implements LineIterator
		{
		private BufferedReader in;
		BuffReadIter(BufferedReader in)
			{
			this.in=in;
			}
		@Override
		protected String advance()
			{
			if(in==null) return null;
			String s;
			try {
				s = in.readLine();
				if(s==null) {CloserUtil.close(this.in);this.in=null;}
				return s;
			} catch (Exception e) {
				throw new RuntimeException(e);
				}
			}
		}
	
	/** Prevent an output stream to be closed. 
	 *  Wrap into a stream that will be flushed instead of being closed 
	 *  For example calling vcfwriter.close in a zip entry. */
	public static OutputStream uncloseableOutputStream(final OutputStream os) {
		return new FilterOutputStream(os) {
			@Override
			public void close() throws IOException {
				this.flush();
				os.flush();
				}
		};
	}
	
	}
