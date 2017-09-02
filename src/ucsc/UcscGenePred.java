package ucsc;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.Reader;
import java.io.UnsupportedEncodingException;
import java.net.URL;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.zip.GZIPInputStream;

import com.google.common.base.Joiner;
import com.google.common.base.Splitter;
import com.google.common.collect.Lists;

import exceptions.InvalidCommandLineException;
import exceptions.InvalidGenomicCoordsException;
import exceptions.InvalidRecordException;
import htsjdk.tribble.index.tabix.TabixFormat;
import samTextViewer.Utils;
import sortBgzipIndex.MakeTabixIndex;
import tracks.IntervalFeature;
import tracks.TrackFormat;

public class UcscGenePred {

	private final String ucscGoldenPath= "http://hgdownload.soe.ucsc.edu/goldenPath";
	private String tabixFile;
	
	/* C O N S T R U C T O R S */
	
	public UcscGenePred(){
		
	}

	/** Fetch genePred data and convert it to tabix indexed GTF. 
	 * @param maxRecords stop reading after this many records. If maxRecrods < 0 read entire data. 
	 * Use this option to test whether connection is suitable genePred input.   
	 * */
	public UcscGenePred(String databaseTableOrFile, int maxRecords) throws IOException, InvalidCommandLineException, InvalidGenomicCoordsException, ClassNotFoundException, InvalidRecordException, SQLException{

		if(maxRecords < 0){
			maxRecords= Integer.MAX_VALUE;
		}
		
		// Prepare reader 
		File genePred= new File(databaseTableOrFile); 

		BufferedReader br; 

		if(genePred.isFile()){
			try{
				InputStream fileStream = new FileInputStream(genePred);
				Reader decoder = new InputStreamReader(new GZIPInputStream(fileStream), "UTF-8");
				br= new BufferedReader(decoder);
			} catch (Exception e){
				br= new BufferedReader(new FileReader(genePred));
			}

		} else if(Utils.urlFileExists(databaseTableOrFile)) {

				br= this.readFromUrl(databaseTableOrFile);
			
		} else {
			// Read from UCSC if not local file or remote file.
			List<String> conn= Lists.newArrayList((Splitter.on(":").omitEmptyStrings().split(databaseTableOrFile)));
			if(conn.size() != 2){
				System.err.println(databaseTableOrFile + " is not a local file and does not look like a string in the form database:table");
				throw new InvalidCommandLineException();
			}
			br= this.readGenePredFromUcsc(conn.get(0), conn.get(1));
		} 
		
		// Prepare output gtf file
		String source= (new File(databaseTableOrFile).getName()).replaceAll(".txt.gz$|.txt$", "");
		File tmpGtf= File.createTempFile("asciigenome.", ".ucsc." + source + ".gtf");
		tmpGtf.deleteOnExit();

		BufferedWriter wr= new BufferedWriter(new FileWriter(tmpGtf));
		
		String line= null;
		int i= 0;
		while((line = br.readLine()) != null && i < maxRecords){
			if(line.startsWith("#")){ // Are comment lines allowed at all?
				continue;
			}
			List<IntervalFeature> gene= this.genePredRecordToGene(line, source);
			i++;
			for(IntervalFeature x : gene){
				wr.write(x.getRaw() + "\n");
			}
		}
		if(maxRecords == Integer.MAX_VALUE){
			// System.err.println(i + " genePred records converted to " + j + " GTF features (" + Math.rint((t1-t0) / 1000) + " s).");
		}
		br.close();
		wr.close();
		
		// Tabix conversion
		String infile= tmpGtf.getAbsolutePath();
		File outfile= new File(infile + ".gz");
		outfile.deleteOnExit();
		File expectedTbi= new File(outfile.getAbsolutePath() + ".tbi"); 
		expectedTbi.deleteOnExit();
		
		new MakeTabixIndex(infile, outfile, TabixFormat.GFF);
		tmpGtf.delete();
		this.tabixFile= outfile.getAbsolutePath();
	}
	
	private BufferedReader readFromUrl(String databaseTableOrFile) throws UnsupportedEncodingException, IOException {

		URL url= new URL(databaseTableOrFile);

		InputStream fileStream = url.openConnection().getInputStream();
		try{
			Reader decoder = new InputStreamReader(new GZIPInputStream(fileStream), "UTF-8");
			BufferedReader br= new BufferedReader(decoder);
			return br;
		} catch (Exception e) {
			BufferedReader br= new BufferedReader(new InputStreamReader(fileStream));
			return br;
		}

	}

	/* M E T H O D S */
	
	protected List<IntervalFeature> genePredRecordToGene(String genePredLine, String source) throws InvalidGenomicCoordsException {
		List<String> genePredList= Lists.newArrayList(Splitter.on("\t").split(genePredLine));
		
		if(genePredList.size() == 16){ // Remove bin field.
			try{
				Integer.parseInt(genePredList.get(0));
				genePredList.remove(0);
			} catch (NumberFormatException e){
				//
			}
		}


		List<IntervalFeature> gene= new ArrayList<IntervalFeature>();

		gene.add(this.getTranscript(genePredList, source));
		
		List<IntervalFeature> exons = this.getExons(genePredList, source);
		gene.addAll(exons);

		// CDS // MEMO: frame is "-1" if no frame for exon. 
		// final List<String> exonFrames= Lists.newArrayList(Splitter.on(",").omitEmptyStrings().split(genePredList.get(14)));

		final int cdsStart= Integer.parseInt(genePredList.get(5)) + 1; // +1 because genePred is 0-based
		final int cdsEnd= Integer.parseInt(genePredList.get(6));
		final String cdsStartStat= genePredList.get(12);
		final String cdsEndStat= genePredList.get(13);
		List<IntervalFeature> cds = this.getCDS(exons, cdsStart, cdsEnd, cdsStartStat, cdsEndStat);
		gene.addAll(cds);
		List<IntervalFeature> startCodonsFwd = this.getStartCodonFwd(cds, cdsStartStat);
		gene.addAll(startCodonsFwd);
		List<IntervalFeature> startCodonsRev = this.getStartCodonRev(cds, cdsEndStat);
		gene.addAll(startCodonsRev);
		//gene.addAll(this.getLeftUTR(exons, cdsStart, cdsEnd));
		//gene.addAll(this.getRightUTR(exons, cdsStart, cdsEnd));
		return gene;
	}
	
	/** @param genePredList List where bin column has been already removed.  
	 * @throws InvalidGenomicCoordsException 
	 * */
	private IntervalFeature getTranscript(List<String>genePredList, String source) throws InvalidGenomicCoordsException{
		
		String[] gff= new String[9];

		gff[0]= genePredList.get(1);
		gff[1]= source;
		gff[2]= "transcript";
		gff[3]= Integer.toString(Integer.parseInt(genePredList.get(3)) + 1);
		gff[4]= genePredList.get(4);
		gff[5]= ".";
		gff[6]= genePredList.get(2); // Strand
		gff[7]= ".";
		gff[8]= "gene_id \"" + genePredList.get(11) + '"' +
				"; transcript_id \"" + genePredList.get(0) + '"' + 
				"; gene_name \"" + genePredList.get(11) + "\";";

		IntervalFeature tx = new IntervalFeature(Joiner.on("\t").join(gff), TrackFormat.GTF, null);
		return tx;
	}
	
	private List<IntervalFeature> getExons(List<String>genePredList, String source) throws InvalidGenomicCoordsException{
		
		int exonCount= Integer.parseInt(genePredList.get(7));
		
		Iterator<String> iter= Splitter.on(",").omitEmptyStrings().split(genePredList.get(8)).iterator();
		int[] exonStarts= new int[exonCount];
		int i= 0;
		while(iter.hasNext()){
			exonStarts[i]= Integer.parseInt(iter.next());
			i++;
		}
		iter= Splitter.on(",").omitEmptyStrings().split(genePredList.get(9)).iterator();
		int[] exonEnds= new int[exonCount];
		i= 0;
		while(iter.hasNext()){
			exonEnds[i]= Integer.parseInt(iter.next());
			i++;
		}
		
		List<IntervalFeature> exons= new ArrayList<IntervalFeature>();
		for(int j = 0; j < exonStarts.length; j++){

			String[] gff= new String[9];

			gff[0]= genePredList.get(1);
			gff[1]= source;
			gff[2]= "exon";
			gff[3]= Integer.toString(exonStarts[j] + 1);
			gff[4]= Integer.toString(exonEnds[j]);
			gff[5]= ".";
			gff[6]= genePredList.get(2); // Strand
			gff[7]= "."; // Frame: Exon do not have frame set. CDS do.
			gff[8]= "gene_id \"" + genePredList.get(11) + '"' +
					"; transcript_id \"" + genePredList.get(0) + '"' + 
					"; exon_number \"" + (j + 1) + '"' +
					"; exon_id \"" + genePredList.get(0) + "." + (j + 1) + '"' + 
					"; gene_name \"" + genePredList.get(11) + "\";";
			IntervalFeature x = new IntervalFeature(Joiner.on("\t").join(gff), TrackFormat.GTF, null);
			exons.add(x);
		}
		return exons;
	}
	
	private List<IntervalFeature> getCDS(List<IntervalFeature> exons, final int cdsStart, final int cdsEnd, String cdsStartStat, String cdsEndStat) throws InvalidGenomicCoordsException{

		if(cdsStart > cdsEnd){ // There are no CDS in this transcript
			return new ArrayList<IntervalFeature>();
		}
		
		// Iterate through exons checking whether at least part of it is containing in the interval cdsStart:cdsEnd.
		// If so, take the exon slice inside the interval cdsStart:cdsEnd and add it to the list of CDSs
		List<IntervalFeature> cds= new ArrayList<IntervalFeature>();
		
		for(int i= 0; i < exons.size(); i++){
			IntervalFeature exon= exons.get(i);
			if(exon.getTo() < cdsStart || exon.getFrom() > cdsEnd){ // Exon is not in interval cdsStart:cdsEnd
				continue;
			}
			int cdsFrom= exon.getFrom();
			if(cdsFrom < cdsStart){ // If only part of exon is CDS 
				cdsFrom= cdsStart;
			}
			int cdsTo= exon.getTo();
			if(cdsTo > cdsEnd){ // If only part of exon is CDS
				cdsTo= cdsEnd;
			}
			
			String attr= Lists.newArrayList(Splitter.on("\t").omitEmptyStrings().split(exon.getRaw())).get(8);

			// Build the interval feature object
			String[] gff= new String[9];

			gff[0]= exons.get(0).getChrom();
			gff[1]= exons.get(0).getSource();
			gff[2]= "CDS";
			gff[3]= Integer.toString(cdsFrom);
			gff[4]= Integer.toString(cdsTo);
			gff[5]= ".";
			gff[6]= String.valueOf(exons.get(0).getStrand());
			gff[7]= "."; // It's unclear to me how frames are assigned so leave it N/A.
			gff[8]= attr;
			IntervalFeature x = new IntervalFeature(Joiner.on("\t").join(gff), TrackFormat.GTF, null);
			cds.add(x);
		}
		// Adjust end of CDS. 
		// The last three bases of the CDS are stop_codon and should not be part of the CDS.
		// For transcripts on +: Subtract 3 from the end of the last CDS (rightmost)
		// For transcripts on -: Add 3 to the start of the first CDS (leftmost). 
		// If such adjustment result in 0 or negative length CDS, remove that CDS altogether 
		// and remove the remainder from the following CDS.
		if(exons.get(0).getStrand() == '+' && cdsEndStat.equals("cmpl")){
			
			IntervalFeature stopCds= cds.get(cds.size() - 1);
			int newStop= stopCds.getTo() - 3;
			int remainder= -(newStop - stopCds.getFrom()); 
			
			if(remainder > 0){
				// If remainder is > 0, this CDS doesn't exist at all and must be removed. This happens if the 
				// stop codon is split across two exons (rare but it happens).
                // We also need to chip off "remainder" from the previous CDS.
				cds.remove(stopCds);
				stopCds= cds.get(cds.size() - 1);
				newStop= stopCds.getTo() - remainder + 1; // I'm not sure why you need +1 to make it work!
			} // In theory you should keep going because also the next CDS might need to be removed!!
			
			if(newStop <= 0 || newStop < stopCds.getFrom()){ // Sanity check
				throw new InvalidGenomicCoordsException();
			} 
			
			// We create an intervalFeature from scratch that will replace the old one.
			List<String> raw= Lists.newArrayList(Splitter.on("\t").omitEmptyStrings().splitToList(stopCds.getRaw()));
			raw.set(4, Integer.toString(newStop)); // Replace end coord
			cds.set(cds.size()-1, new IntervalFeature(Joiner.on("\t").join(raw), TrackFormat.GTF, null)); // Replace last element.
		
		} else if(exons.get(0).getStrand() == '-' && cdsStartStat.equals("cmpl")){
			// same as above. This time apply to first CDS whose start has to be increased by 3
			IntervalFeature stopCds= cds.get(0);
			int newStop= stopCds.getFrom() + 3;
			int remainder= newStop - stopCds.getTo(); 

			if(remainder > 0){
				// If remainder is >= 0, this CDS doesn't exist at all and must be removed. This happens if the 
				// stop codon is split across two exons (rare but it happens).
                // We also need to chip off "remainder" from the next CDS.
				cds.remove(stopCds);
				stopCds= cds.get(0);
				newStop= stopCds.getFrom() + remainder - 1; // Not sure why -1 works!  
			} // In theory you should keep going because also the next CDS might need to be removed!!

			if(newStop <= 0 || newStop > stopCds.getTo()){ // Sanity check
				throw new InvalidGenomicCoordsException();
			} 
			// We create an intervalFeature from scratch that will replace the old one.
			List<String> raw= Lists.newArrayList(Splitter.on("\t").omitEmptyStrings().splitToList(stopCds.getRaw()));
			raw.set(3, Integer.toString(newStop)); // Replace start coord
			cds.set(0, new IntervalFeature(Joiner.on("\t").join(raw), TrackFormat.GTF, null)); // Replace last element.
			
		}
		return cds;
	}
	
	private List<IntervalFeature> getStartCodonFwd(List<IntervalFeature> cds, String cdsStartStat) throws InvalidGenomicCoordsException{

		List<IntervalFeature> codons= new ArrayList<IntervalFeature>();
		if(cds.size() == 0 || cds.get(0).getStrand() == '-' || ! cdsStartStat.equals("cmpl") ){ 
			return codons; 
		}
		
//		CC   CCCCC
//		AA   A

		// For tx on +: Get the first three bases of the first CDS,. If CDS length is < 3, get the remainder from
		// next CDS
		IntervalFeature c= cds.get(0);
		int cdnEnd= c.getFrom() + 2;
		int remainder= cdnEnd - c.getTo();

		if(remainder <= 0){ // Codon is fully contained in CDS. Easy.
			String[] gff= new String[9];
			gff[0]= c.getChrom();
			gff[1]= c.getSource();
			gff[2]= "start_codon";
			gff[3]= Integer.toString(c.getFrom());
			gff[4]= Integer.toString(cdnEnd);
			gff[5]= ".";
			gff[6]= String.valueOf(c.getStrand());
			gff[7]= ".";
			gff[8]= Lists.newArrayList(Splitter.on("\t").omitEmptyStrings().split(c.getRaw())).get(8);
			IntervalFeature x = new IntervalFeature(Joiner.on("\t").join(gff), TrackFormat.GTF, null);
			codons.add(x);
			return codons;
		}
		// Add to list this partial codon and get "remainder" from next cds
		String[] gff= new String[9];
		gff[0]= c.getChrom();
		gff[1]= c.getSource();
		gff[2]= "start_codon";
		gff[3]= Integer.toString(c.getFrom());
		gff[4]= Integer.toString(cdnEnd - remainder);
		gff[5]= ".";
		gff[6]= String.valueOf(c.getStrand());
		gff[7]= ".";
		gff[8]= Lists.newArrayList(Splitter.on("\t").omitEmptyStrings().split(c.getRaw())).get(8);
		IntervalFeature x1 = new IntervalFeature(Joiner.on("\t").join(gff), TrackFormat.GTF, null);
		codons.add(x1);
		
		c= cds.get(1);
		cdnEnd= c.getFrom() + remainder - 1;

		gff= new String[9];
		gff[0]= c.getChrom();
		gff[1]= c.getSource();
		gff[2]= "start_codon";
		gff[3]= Integer.toString(c.getFrom());
		gff[4]= Integer.toString(cdnEnd);
		gff[5]= ".";
		gff[6]= String.valueOf(c.getStrand());
		gff[7]= ".";
		gff[8]= Lists.newArrayList(Splitter.on("\t").omitEmptyStrings().split(c.getRaw())).get(8);
		IntervalFeature x2 = new IntervalFeature(Joiner.on("\t").join(gff), TrackFormat.GTF, null);
		codons.add(x2);
		
		return codons;
	}

	private List<IntervalFeature> getStartCodonRev(List<IntervalFeature> cds, String cdsEndStat) throws InvalidGenomicCoordsException{

		List<IntervalFeature> codons= new ArrayList<IntervalFeature>();
		if(cds.size() == 0 || cds.get(0).getStrand() == '+' || ! cdsEndStat.equals("cmpl")){ 
			return codons; 
		}
		
//		cccccccc      cc
//		       a      aa
			
		// Get the last three bases of the last CDS. If CDS length is < 3, get the remainder from
		// previous CDS
		IntervalFeature c= cds.get(cds.size()-1);
		int cdnStart= c.getTo() - 2;
		int remainder= c.getFrom() - cdnStart;

		if(remainder <= 0){ // Codon is fully contained in CDS. Easy.
			String[] gff= new String[9];
			gff[0]= c.getChrom();
			gff[1]= c.getSource();
			gff[2]= "start_codon";
			gff[3]= Integer.toString(cdnStart);
			gff[4]= Integer.toString(c.getTo());
			gff[5]= ".";
			gff[6]= String.valueOf(c.getStrand());
			gff[7]= ".";
			gff[8]= Lists.newArrayList(Splitter.on("\t").omitEmptyStrings().split(c.getRaw())).get(8);
			IntervalFeature x = new IntervalFeature(Joiner.on("\t").join(gff), TrackFormat.GTF, null);
			codons.add(x);
			return codons;
		}
		// Add to list this partial codon and get "remainder" from next cds
		String[] gff= new String[9];
		gff[0]= c.getChrom();
		gff[1]= c.getSource();
		gff[2]= "start_codon";
		gff[3]= Integer.toString(c.getFrom());
		gff[4]= Integer.toString(c.getTo());
		gff[5]= ".";
		gff[6]= String.valueOf(c.getStrand());
		gff[7]= ".";
		gff[8]= Lists.newArrayList(Splitter.on("\t").omitEmptyStrings().split(c.getRaw())).get(8);
		IntervalFeature x1 = new IntervalFeature(Joiner.on("\t").join(gff), TrackFormat.GTF, null);
		codons.add(x1);
		
		c= cds.get(cds.size() - 2);
		cdnStart= c.getTo() - remainder + 1;

		gff= new String[9];
		gff[0]= c.getChrom();
		gff[1]= c.getSource();
		gff[2]= "start_codon";
		gff[3]= Integer.toString(cdnStart);
		gff[4]= Integer.toString(c.getTo());
		gff[5]= ".";
		gff[6]= String.valueOf(c.getStrand());
		gff[7]= ".";
		gff[8]= Lists.newArrayList(Splitter.on("\t").omitEmptyStrings().split(c.getRaw())).get(8);
		IntervalFeature x2 = new IntervalFeature(Joiner.on("\t").join(gff), TrackFormat.GTF, null);
		codons.add(x2);
		
		return codons;
	}

//	private List<IntervalFeature> getStartCodonRev(List<IntervalFeature> cds){
//		
//	}
	
	private List<IntervalFeature> getLeftUTR(List<IntervalFeature> exons, final int cdsStart, final int cdsEnd) throws InvalidGenomicCoordsException{
		List<IntervalFeature> utr= new ArrayList<IntervalFeature>();
		if(cdsStart > cdsEnd){ // There is no UTR in this transcirpt 
			return utr;
		}
// Left UTR goes from txStart to cdsStart
//		EEEEEEEEEEE          EEEEEEEEEEE
//		                         CCCCCCC
//		UUUUUUUUUUU          UUUU
//		X                        X
		
		for(IntervalFeature exon : exons){
			
			if(exon.getFrom() >= cdsStart){ // There is no UTR 
				break;
			}
			
			int utrExonEnd= exon.getTo(); // Use this if the exon is completely to the left of cdsStart 
			if(exon.getFrom() < cdsStart && exon.getTo() >= cdsStart){ // Is the exon containing the cdsStart?
				utrExonEnd=  cdsStart - 1;
			} 
			String[] gff= new String[9];

			gff[0]= exon.getChrom();
			gff[1]= exon.getSource();
			gff[2]= exon.getStrand() == '+' ? "5UTR" : "3UTR";
			gff[3]= Integer.toString(exon.getFrom());
			gff[4]= Integer.toString(utrExonEnd);
			gff[5]= ".";
			gff[6]= String.valueOf(exons.get(0).getStrand());
			gff[7]= ".";
			gff[8]= Lists.newArrayList(Splitter.on("\t").omitEmptyStrings().split(exon.getRaw())).get(8);
			IntervalFeature x = new IntervalFeature(Joiner.on("\t").join(gff), TrackFormat.GTF, null);
			utr.add(x);
		}
		return utr; 		
	}
		

	private List<IntervalFeature> getRightUTR(List<IntervalFeature> exons, final int cdsStart, final int cdsEnd) throws InvalidGenomicCoordsException{
		List<IntervalFeature> utr= new ArrayList<IntervalFeature>();
		if(cdsStart > cdsEnd){ // There is no UTR in this tranx
			return utr;
		}
		
// Right UTR goes from cdsEnd to txEnd
//EEEE  EEEEEEEEEEE          EEEEEEEEEEE		                         
//		      UUUUU          UUUUUUUUUUU
//            X                        X
		
		for(IntervalFeature exon : exons){
			
			if(exon.getTo() <= cdsEnd){ // Exon is fully the left of cdsEnd, ie. is CDS
				continue;
			}
			
			int utrExonStart= exon.getFrom(); // Use this if the exon is completely to the right of cdsEnd
			if(exon.getFrom() < cdsEnd && exon.getTo() > cdsEnd){ // Is the exon containing the cdsStart?
				utrExonStart=  cdsEnd - 1;
			} 
			String[] gff= new String[9];

			gff[0]= exon.getChrom();
			gff[1]= exon.getSource();
			gff[2]= exon.getStrand() == '+' ? "3UTR" : "5UTR";
			gff[3]= Integer.toString(utrExonStart);
			gff[4]= Integer.toString(exon.getTo());
			gff[5]= ".";
			gff[6]= String.valueOf(exons.get(0).getStrand());
			gff[7]= ".";
			gff[8]= Lists.newArrayList(Splitter.on("\t").omitEmptyStrings().split(exon.getRaw())).get(8);;
			IntervalFeature x = new IntervalFeature(Joiner.on("\t").join(gff), TrackFormat.GTF, null);
			utr.add(x);
		}
		return utr; 		
	}
	
	/** Download table from database in UCSC. Return the tmp file to which the file is downloaded.
	 * @throws IOException 
	 * */
	private BufferedReader readGenePredFromUcsc(String database, String table) throws IOException {
		return this.readFromUrl(this.ucscGoldenPath + "/" + database + "/database/" + table + ".txt.gz");
	}

	public String getTabixFile() {
		return this.tabixFile;
	}
	
}
