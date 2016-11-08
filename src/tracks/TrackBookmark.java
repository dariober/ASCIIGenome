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

import com.google.common.base.Joiner;
import com.google.common.base.Splitter;
import com.google.common.collect.Lists;

import exceptions.InvalidGenomicCoordsException;
import exceptions.InvalidRecordException;
import htsjdk.tribble.index.tabix.TabixFormat;
import htsjdk.tribble.readers.TabixReader;
import samTextViewer.GenomicCoords;
import sortBgzipIndex.MakeTabixIndex;

public class TrackBookmark extends TrackIntervalFeature {
	
	private final String trackName= "Bookmarks"; 
	
	public TrackBookmark(GenomicCoords gc, String nameForBookmark) throws IOException, ClassNotFoundException, InvalidRecordException, SQLException, InvalidGenomicCoordsException{
		super(gc);
				
		this.setTrackTag(trackName);
		
		
		// Prepare bookmark file
		// =====================
	    // First write out the current position as plain text. Then gzip and index.
		File bookmarkPlain= File.createTempFile("asciigenome.bookmarks.", ".bed"); // new File("bookmark");
		bookmarkPlain.deleteOnExit();
		BufferedWriter wr = new BufferedWriter(new FileWriter(bookmarkPlain));
		wr.write(gc.getChrom() + "\t" + (gc.getFrom() - 1) + "\t" + 
				gc.getTo() + "\t" + nameForBookmark + "\n");
		wr.close();

		File bookmark= new File(bookmarkPlain + ".gz");
		bookmark.deleteOnExit();
		(new File(bookmark.getAbsolutePath() + ".tbi")).deleteOnExit();
		this.setFilename(bookmark.getAbsolutePath());
		this.setWorkFilename(bookmark.getAbsolutePath());
		
		new MakeTabixIndex(bookmarkPlain.getAbsolutePath(), bookmark, TabixFormat.BED);
		bookmarkPlain.delete();
		
		this.tabixReader= new TabixReader(bookmark.getAbsolutePath());
		this.setType(TrackFormat.BED);
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
	public void addBookmark(String nameForBookmark) throws IOException, ClassNotFoundException, InvalidRecordException, SQLException, InvalidGenomicCoordsException{
		
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
			wr.write(line + "\n");
		}
		// Add new bookamrk
		wr.write(this.positionToBedLine(nameForBookmark) + "\n");
		br.close();
		wr.close();

		// Recompress and index replacing the original bgzip file
		new MakeTabixIndex(plainNew.getAbsolutePath(), new File(this.getWorkFilename()), TabixFormat.BED);
		plainNew.delete();
		this.tabixReader= new TabixReader(this.getWorkFilename());
		// Update track.
		this.update();
	}
	
	/** Convert current position to printable bed line with given feature name.
	 * */
	private String positionToBedLine(String bedFeatureName){
		return this.getGc().getChrom() + "\t" + (this.getGc().getFrom() - 1) + "\t" + this.getGc().getTo() + "\t" +
				bedFeatureName;
	}

	/** Remove the bookmark matching the exact coordinates of the current position. 
	 * Bookmarks partially overlapping are not removed. 
	 * @throws IOException 
	 * @throws SQLException 
	 * @throws InvalidRecordException 
	 * @throws ClassNotFoundException 
	 * @throws InvalidGenomicCoordsException 
	 * */
	public void removeBookmark() throws IOException, ClassNotFoundException, InvalidRecordException, SQLException, InvalidGenomicCoordsException {
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
				
				if(this.getGc().getChrom().equals(pos.get(0)) && 
				   (this.getGc().getFrom() - 1) == Integer.parseInt(pos.get(1)) && 
				   this.getGc().getTo() == Integer.parseInt(pos.get(2))){
					continue;
				}
			}
			wr.write(line + "\n");
		}
		wr.close();
		br.close();
		
		// Recompress and index replacing the original bgzip file
		new MakeTabixIndex(plainNew.getAbsolutePath(), new File(this.getWorkFilename()), TabixFormat.BED);
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
		
		BufferedWriter wr = new BufferedWriter(new FileWriter(filename, append));
		
		InputStream fileStream = new FileInputStream(this.getWorkFilename());
		Reader decoder = new InputStreamReader(new GZIPInputStream(fileStream), "UTF-8");
		BufferedReader br= new BufferedReader(decoder);

		String line;
		while( (line = br.readLine()) != null) {
			wr.write(line + "\n");
		}
		wr.close();
		br.close();
		
	}
	
	/** Generate a string that parsed by -x/--exec re-creates the current bookmarks.
	 * @throws IOException 
	 * @throws UnsupportedEncodingException 
	 * */
	@Override
	public String settingsToString() throws UnsupportedEncodingException, IOException{
		List<String> set= new ArrayList<String>();

		InputStream fileStream = new FileInputStream(this.getWorkFilename());
		Reader decoder = new InputStreamReader(new GZIPInputStream(fileStream), "UTF-8");
		BufferedReader br= new BufferedReader(decoder);

		String rawLine;
		while( (rawLine = br.readLine()) != null) {
			set.add(this.rawLineToBookmarkCmd(rawLine));
		}
		br.close();
		
		set.add("colorTrack " + this.getTitleColour() + " ^" + trackName);
		set.add("trackHeight " + this.getyMaxLines() + " ^" + trackName);
		set.add("grep -i " + this.getShowRegex() + " -e " + this.getHideRegex() + " ^" + trackName);
		if(this.isHideTitle()){
			set.add("hideTitle ^" + trackName);
		}
		
		return Joiner.on(" && ").join(set);
	}
	
	/** Parse a raw line read from bed file to return a command string that reproduces this 
	 * bookmark entry. 
	 * */
	private String rawLineToBookmarkCmd(String rawLine){
		List<String>line= Lists.newArrayList(Splitter.on("\t").split(rawLine));
		String cmd= "goto " + line.get(0) + ":" + (Integer.parseInt(line.get(1))+1) + "-" + line.get(2) + " && ";
		cmd += "bookmark " + line.get(3);
		return cmd;
	}

	public List<String> asList() throws UnsupportedEncodingException, IOException {
		List<String> marks= new ArrayList<String>();
		
		InputStream fileStream = new FileInputStream(this.getWorkFilename());
		Reader decoder = new InputStreamReader(new GZIPInputStream(fileStream), "UTF-8");
		BufferedReader br= new BufferedReader(decoder);

		String line;
		int i= 1;
		while( (line = br.readLine()) != null) {
			List<String> lst= Lists.newArrayList(Splitter.on("\t").omitEmptyStrings().split(line));
			String reg= lst.get(0) + ":" + lst.get(1) + "-" + lst.get(2); 
			line= i + ":\t" + reg + "\t" + line;
			marks.add(line + "\t");
			i++;
		}
		br.close();

		return marks;
		
	}
}
