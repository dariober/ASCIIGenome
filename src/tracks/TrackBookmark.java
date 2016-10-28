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
import java.sql.SQLException;
import java.util.zip.GZIPInputStream;

import exceptions.InvalidGenomicCoordsException;
import exceptions.InvalidRecordException;
import htsjdk.tribble.index.tabix.TabixFormat;
import htsjdk.tribble.readers.TabixReader;
import samTextViewer.GenomicCoords;
import sortBgzipIndex.MakeTabixIndex;

public class TrackBookmark extends TrackIntervalFeature {
	
	
	public TrackBookmark(GenomicCoords gc, String nameForBookmark) throws IOException, ClassNotFoundException, InvalidRecordException, SQLException, InvalidGenomicCoordsException{
		super(gc);
				
		this.setGc(gc);
		this.setTrackTag("Bookmarks");
		
		
		// Prepare bookmark file
		// =====================
	    // First write out the current position as plain text. Then gzip and index.
		File bookmarkPlain= File.createTempFile("asciigenome.bookmarks.", ".bed"); // new File("bookmark");
		bookmarkPlain.deleteOnExit();
		BufferedWriter wr = new BufferedWriter(new FileWriter(bookmarkPlain));
		wr.write(this.positionToBedLine(nameForBookmark) + "\n");
		wr.close();

		File bookmark= new File(bookmarkPlain + ".gz");
		bookmark.deleteOnExit();
		(new File(bookmark.getAbsolutePath() + ".tbi")).deleteOnExit();
		this.setFilename(bookmark.getAbsolutePath());
		
		new MakeTabixIndex(bookmarkPlain.getAbsolutePath(), bookmark, TabixFormat.BED);
		bookmarkPlain.delete();
		
		this.tabixReader= new TabixReader(bookmark.getAbsolutePath());
		this.setType(TrackFormat.BED);
		this.update();
	}
		
	/** Add current genomic position to track.  
	 * @throws IOException 
	 * @throws SQLException 
	 * @throws InvalidRecordException 
	 * @throws ClassNotFoundException 
	 * @throws InvalidGenomicCoordsException 
	 * */
	public void add(String nameForBookmark) throws IOException, ClassNotFoundException, InvalidRecordException, SQLException, InvalidGenomicCoordsException{
		
		// Stub
		File plainNew= new File(this.getFilename() + "2");
		plainNew.deleteOnExit();
		BufferedWriter wr = new BufferedWriter(new FileWriter(plainNew));
		
		GZIPInputStream gzipStream;
		InputStream fileStream = new FileInputStream(this.getFilename());
		gzipStream = new GZIPInputStream(fileStream);
		Reader decoder = new InputStreamReader(gzipStream, "UTF-8");
		BufferedReader br= new BufferedReader(decoder);

		String line;
		while( (line = br.readLine()) != null) {
			wr.write(line + "\n");
		}
		wr.write(this.positionToBedLine(nameForBookmark) + "\n");
		br.close();
		wr.close();
		new MakeTabixIndex(plainNew.getAbsolutePath(), new File(this.getFilename()), TabixFormat.BED);
		plainNew.delete();
		this.tabixReader= new TabixReader(this.getFilename());
		this.update();
	}
	
	/** Convert current position to printable bed line with given feature name.
	 * */
	private String positionToBedLine(String bedFeatureName){
		return this.getGc().getChrom() + "\t" + (this.getGc().getFrom() - 1) + "\t" + this.getGc().getTo() + "\t" +
				bedFeatureName;
	}
	
}
