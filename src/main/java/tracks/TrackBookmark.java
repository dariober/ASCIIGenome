package tracks;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.Reader;
import java.io.UnsupportedEncodingException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.List;
import java.util.zip.GZIPInputStream;

import com.google.common.base.Splitter;
import com.google.common.collect.Lists;

import exceptions.InvalidGenomicCoordsException;
import exceptions.InvalidRecordException;
import htsjdk.tribble.index.tabix.TabixFormat;
import htsjdk.tribble.readers.TabixReader;
import samTextViewer.GenomicCoords;
import samTextViewer.Utils;
import sortBgzipIndex.MakeTabixIndex;

public class TrackBookmark extends TrackIntervalFeature {
	
	private final String trackName= "Bookmarks"; 
	
	public TrackBookmark(GenomicCoords gc, String nameForBookmark) throws IOException, ClassNotFoundException, InvalidRecordException, SQLException, InvalidGenomicCoordsException{
		super(gc);
				
		this.setTrackTag(trackName);
				
		// Prepare bookmark file
		// =====================
	    // First write out the current position as plain text. Then gzip and index.
		File bookmarkPlain= Utils.createTempFile(".asciigenome.bookmarks.", ".gff", true);
		BufferedWriter wr = new BufferedWriter(new FileWriter(bookmarkPlain));
		wr.write(this.positionToGffLine(gc, nameForBookmark) + "\n"); // gc.getChrom() + "\t" + (gc.getFrom() - 1) + "\t" + gc.getTo() + "\t" + nameForBookmark + "\n");
		wr.close();

		File bookmark= new File(bookmarkPlain + ".gz");
		bookmark.deleteOnExit();
		(new File(bookmark.getAbsolutePath() + ".tbi")).deleteOnExit();
		this.setFilename(bookmark.getAbsolutePath());
		this.setWorkFilename(bookmark.getAbsolutePath());
		
		new MakeTabixIndex(bookmarkPlain.getAbsolutePath(), bookmark, TabixFormat.GFF);
		bookmarkPlain.delete();
		
		this.setTabixReader(new TabixReader(bookmark.getAbsolutePath()));
		this.setTrackFormat(TrackFormat.GTF);
		this.setGc(gc);
	}
		
	/** Add current genomic position to track.  
	 * @throws IOException 
	 * @throws SQLException 
	 * @throws InvalidRecordException 
	 * @throws ClassNotFoundException 
	 * @throws InvalidGenomicCoordsException 
	 * */
	@Override
	public void addBookmark(GenomicCoords gc, String nameForBookmark) throws IOException, ClassNotFoundException, InvalidRecordException, SQLException, InvalidGenomicCoordsException{
		
		// Adding a bookmark means re-preparing the bgzip file again. 
		
		// First write the current position to a new tmp file, then append to this file
		// all the bookmarks previously added. Then bgzip and compress replacing the old bookmark file.
		File plainNew= new File(this.getWorkFilename() + ".add");
		plainNew.deleteOnExit();
		BufferedWriter wr = new BufferedWriter(new FileWriter(plainNew));
		
		InputStream fileStream = new FileInputStream(this.getWorkFilename());
		Reader decoder = new InputStreamReader(new GZIPInputStream(fileStream), "UTF-8");
		BufferedReader br= new BufferedReader(decoder);

		String line;
		while( (line = br.readLine()) != null) {
			if(line.contains("\t__ignore_me__")){ // Hack to circumvent issue #38
				continue;
			}			
			wr.write(line + "\n");
		}
		// Add new bookamrk
		wr.write(this.positionToGffLine(gc, nameForBookmark) + "\n");
		br.close();
		wr.close();

		// Recompress and index replacing the original bgzip file
		new MakeTabixIndex(plainNew.getAbsolutePath(), new File(this.getWorkFilename()), TabixFormat.GFF);
		plainNew.delete();
		this.tabixReader= new TabixReader(this.getWorkFilename());
		// Update track.
		this.update();
	}
	
	/** Convert current position to printable GFF line with given feature name.
	 * */
	private String positionToGffLine(GenomicCoords gc, String nameForBookmark){
		
		nameForBookmark= nameForBookmark.replaceAll("\"", "_"); // This is a hack to prevent troubles getting back name from GFF raw line.
		
		return gc.getChrom() 
				+ "\tASCIIGenome\tbookmark\t" + gc.getFrom() + "\t" 
				+ gc.getTo() 
				+ "\t.\t.\t.\tID=\"" + nameForBookmark + "\"";
	}

	
	/** Remove the bookmark matching the exact coordinates of the current position. 
	 * Bookmarks partially overlapping are not removed. 
	 * @param bookmarkRegion 
	 * @throws IOException 
	 * @throws SQLException 
	 * @throws InvalidRecordException 
	 * @throws ClassNotFoundException 
	 * @throws InvalidGenomicCoordsException 
	 * */
	public void removeBookmark(GenomicCoords bookmarkRegion) throws IOException, ClassNotFoundException, InvalidRecordException, SQLException, InvalidGenomicCoordsException {
		// To remove a bookmark, iterate through the bgzip file writing records to a tmp file.
		// The record(s) matching this position is not written. 
		// 
		File plainNew= new File(this.getWorkFilename() + ".remove");
		plainNew.deleteOnExit();
		BufferedWriter wr = new BufferedWriter(new FileWriter(plainNew));
		
		InputStream fileStream = new FileInputStream(this.getWorkFilename());
		Reader decoder = new InputStreamReader(new GZIPInputStream(fileStream), "UTF-8");
		BufferedReader br= new BufferedReader(decoder);

		String line;
		while( (line = br.readLine()) != null) {
			
			if( ! line.startsWith("#")){ // In case of comment lines
				List<String>pos= Lists.newArrayList(Splitter.on("\t").split(line));
				
				if(bookmarkRegion.getChrom().equals(pos.get(0)) && 
				   bookmarkRegion.getFrom() == Integer.parseInt(pos.get(3)) && 
				   bookmarkRegion.getTo() == Integer.parseInt(pos.get(4))){
					continue;
				}
			}
			wr.write(line + "\n");
		}
		wr.close();
		br.close();
		
		// Recompress and index replacing the original bgzip file
		new MakeTabixIndex(plainNew.getAbsolutePath(), new File(this.getWorkFilename()), TabixFormat.GFF);
		plainNew.delete();
		this.tabixReader= new TabixReader(this.getWorkFilename());
		// Update track.
		this.update();
	}

	/** Save track to file.
	 * @param append 
	 * @throws IOException 
	 * */
	public void save(String filename, boolean append) throws IOException {
		
		filename= Utils.tildeToHomeDir(filename);
		
		BufferedWriter wr = new BufferedWriter(new FileWriter(filename, append));
		
		InputStream fileStream = new FileInputStream(this.getWorkFilename());
		Reader decoder = new InputStreamReader(new GZIPInputStream(fileStream), "UTF-8");
		BufferedReader br= new BufferedReader(decoder);

		String line;
		while( (line = br.readLine()) != null) {
			if(line.contains("\t__ignore_me__")){ // Hack to circumvent issue #38
				continue;
			}
			wr.write(line + "\n");
		}
		wr.close();
		br.close();
		
	}
		
	public List<String> asList() throws UnsupportedEncodingException, IOException {
		List<String> marks= new ArrayList<String>();
		
		InputStream fileStream = new FileInputStream(this.getWorkFilename());
		Reader decoder = new InputStreamReader(new GZIPInputStream(fileStream), "UTF-8");
		BufferedReader br= new BufferedReader(decoder);

		String line;
		int i= 1;
		while( (line = br.readLine()) != null) {
			if(line.contains("\t__ignore_me__")){ // Hack to circumvent issue #38
				continue;
			}			
			List<String> lst= Lists.newArrayList(Splitter.on("\t").omitEmptyStrings().split(line));
			String reg= lst.get(0) + ":" + Integer.parseInt(lst.get(3)) + "-" + lst.get(4); 
			line= i + ":\t" + reg + "\t" + line;
			marks.add(line);
			i++;
		}
		br.close();

		return marks;
		
	}
}
