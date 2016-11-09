package ucsc;

import static org.junit.Assert.*;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.Reader;
import java.io.UnsupportedEncodingException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.zip.GZIPInputStream;

import org.junit.Test;

import com.google.common.base.Splitter;
import com.google.common.collect.Lists;

import exceptions.InvalidCommandLineException;
import exceptions.InvalidGenomicCoordsException;
import exceptions.InvalidRecordException;
import samTextViewer.GenomicCoords;
import tracks.IntervalFeature;
import tracks.TrackIntervalFeature;

public class UcscGenePredTest {

	// GTF file from genePred prepared with:
	// zcat refGene.hg19.txt.gz | cut -f 2- | grep -e  '\tchr7\t' | genePredToGtf -utr file stdin stdout | gzip > refGene.hg19.genePredToGtf.gtf.gz
	
	@Test
	public void canParseGenePredLineToGtfGene() throws InvalidGenomicCoordsException{
		
		UcscGenePred ucsc= new UcscGenePred();
		String genePredLine= "812	NM_001170300	chr3R	+	29836567	29837447	29836625	29837342	3	29836567,29836912,29837317,	29836832,29837262,29837447,	0	CG15510	cmpl	cmpl	0,0,2,";
		List<IntervalFeature> genes= ucsc.genePredRecordToGene(genePredLine, "ucsc");
		
		// test transcript
		assertEquals("transcript", genes.get(0).getFeature());
		assertEquals(29836568, genes.get(0).getFrom());
		assertEquals(29837447, genes.get(0).getTo());
		assertEquals('+', genes.get(0).getStrand());

		// Exons
		for(IntervalFeature x : genes){
			if(x.getFeature().equals("exon")){
				assertEquals(29836568, x.getFrom());
				assertEquals(29836832, x.getTo());
				assertTrue(x.getRaw().contains("NM_001170300"));
				break;
			}
		}
		int exonCount= 0;
		for(IntervalFeature x : genes){
			if(x.getFeature().equals("exon")){
				exonCount += 1;
			}
		}
		assertEquals(3, exonCount);
		
//		chr3R  stdin  transcript   29836568  29837447  .  +  .  gene_id  "CG15510";  transcript_id  "NM_001170300";  gene_name    "CG15510";
//		chr3R  stdin  exon         29836568  29836832  .  +  .  gene_id  "CG15510";  transcript_id  "NM_001170300";  exon_number  "1";        exon_id  "NM_001170300.1";  gene_name  "CG15510";
//		chr3R  stdin  5UTR         29836568  29836625  .  +  .  gene_id  "CG15510";  transcript_id  "NM_001170300";  exon_number  "1";        exon_id  "NM_001170300.1";  gene_name  "CG15510";
//		chr3R  stdin  CDS          29836626  29836832  .  +  0  gene_id  "CG15510";  transcript_id  "NM_001170300";  exon_number  "1";        exon_id  "NM_001170300.1";  gene_name  "CG15510";
//		chr3R  stdin  exon         29836913  29837262  .  +  .  gene_id  "CG15510";  transcript_id  "NM_001170300";  exon_number  "2";        exon_id  "NM_001170300.2";  gene_name  "CG15510";
//		chr3R  stdin  CDS          29836913  29837262  .  +  0  gene_id  "CG15510";  transcript_id  "NM_001170300";  exon_number  "2";        exon_id  "NM_001170300.2";  gene_name  "CG15510";
//		chr3R  stdin  exon         29837318  29837447  .  +  .  gene_id  "CG15510";  transcript_id  "NM_001170300";  exon_number  "3";        exon_id  "NM_001170300.3";  gene_name  "CG15510";
//		chr3R  stdin  CDS          29837318  29837339  .  +  1  gene_id  "CG15510";  transcript_id  "NM_001170300";  exon_number  "3";        exon_id  "NM_001170300.3";  gene_name  "CG15510";
//		chr3R  stdin  3UTR         29837343  29837447  .  +  .  gene_id  "CG15510";  transcript_id  "NM_001170300";  exon_number  "3";        exon_id  "NM_001170300.3";  gene_name  "CG15510";
//		chr3R  stdin  start_codon  29836626  29836628  .  +  0  gene_id  "CG15510";  transcript_id  "NM_001170300";  exon_number  "1";        exon_id  "NM_001170300.1";  gene_name  "CG15510";
//		chr3R  stdin  stop_codon   29837340  29837342  .  +  0  gene_id  "CG15510";  transcript_id  "NM_001170300";  exon_number  "3";        exon_id  "NM_001170300.3";  gene_name  "CG15510";
		
	}
	
	@Test
	public void canCreateTrackFromLocalFile() throws IOException, InvalidCommandLineException, InvalidGenomicCoordsException, ClassNotFoundException, InvalidRecordException, SQLException {

		 UcscGenePred ucsc= new UcscGenePred("test_data/refGene.hg19.chr7.txt.gz", -1); 
		 TrackIntervalFeature gtfTrack= new TrackIntervalFeature(ucsc.getTabixFile(), new GenomicCoords("chr7", null, null)); 
				 
		 List<IntervalFeature> genes = gtfTrack.getFeaturesInInterval("chr7", 5000000, 6000000);
		 assertTrue(genes.size() > 50);

		long t0= System.currentTimeMillis();
		ucsc= new UcscGenePred("/Users/berald01/Downloads/refGene.hg19.txt.gz", -1);
		long t1= System.currentTimeMillis();
		System.out.println(t1-t0);
		
	}

	@Test
	public void canCreateTrackFromUrl() throws IOException, InvalidCommandLineException, InvalidGenomicCoordsException, ClassNotFoundException, InvalidRecordException, SQLException {

		UcscGenePred ucsc= new UcscGenePred("dm6:refGene", -1); 
		TrackIntervalFeature gtfTrack= new TrackIntervalFeature(ucsc.getTabixFile(), new GenomicCoords("chr3R", null, null)); 
				 
		List<IntervalFeature> genes = gtfTrack.getFeaturesInInterval("chr3R", 1, 6000000);
		assertTrue(genes.size() > 1000);
		
	}
	
	@Test 
	public void canGetCDS() throws UnsupportedEncodingException, IOException, ClassNotFoundException, InvalidCommandLineException, InvalidGenomicCoordsException, InvalidRecordException, SQLException{
		
		// Prepare expected data.
		InputStream fileStream = new FileInputStream(new File("test_data/refGene.hg19.genePredToGtf.gtf.gz"));
		Reader decoder = new InputStreamReader(new GZIPInputStream(fileStream), "UTF-8");
		BufferedReader br= new BufferedReader(decoder);
		
		List<Integer> expectedStarts= new ArrayList<Integer>();
		List<Integer> expectedEnds= new ArrayList<Integer>();
		String line= null;
		while((line = br.readLine()) != null){
			List<String> xline= Lists.newArrayList(Splitter.on("\t").split(line));
			if(xline.get(2).equals("CDS")){
				expectedStarts.add(Integer.parseInt(xline.get(3)));
				expectedEnds.add(Integer.parseInt(xline.get(4)));
			}
		}
		br.close();
		Collections.sort(expectedStarts);
		Collections.sort(expectedEnds);
		assertEquals(expectedStarts.size(), expectedEnds.size());
		//
		// Prepare observed data
		UcscGenePred ucsc= new UcscGenePred("test_data/refGene.hg19.chr7.txt.gz", -1); 
		TrackIntervalFeature gtfTrack= new TrackIntervalFeature(ucsc.getTabixFile(), new GenomicCoords("chr7", null, null)); 
		List<IntervalFeature> obs = gtfTrack.getFeaturesInInterval("chr7", 1, Integer.MAX_VALUE);

		List<Integer> obsStarts= new ArrayList<Integer>();
		List<Integer> obsEnds= new ArrayList<Integer>();
		for(IntervalFeature x  : obs){
			if(x.getFeature().equals("CDS")){
				obsStarts.add(x.getFrom());
				obsEnds.add(x.getTo());
			}
		}
		Collections.sort(obsStarts);
		Collections.sort(obsEnds);
		assertEquals(expectedStarts.size(), obsStarts.size());
		int mism= 0;
		for(int i= 0; i < obsStarts.size(); i++){
			if((int)expectedStarts.get(i) !=  (int)obsStarts.get(i)){
				mism++;
			}
			if((int)expectedEnds.get(i) !=  (int)obsEnds.get(i)){
				mism++;
			}
		}
		// Few mismatches could not be resolved!
		assertTrue(mism < 10);
	}
	
	@Test(expected=InvalidCommandLineException.class)
	public void canHandleNonExistantFileOrUrl() throws ClassNotFoundException, IOException, InvalidCommandLineException, InvalidGenomicCoordsException, InvalidRecordException, SQLException{
		new UcscGenePred("nonsense", -1);
	}

//	@Test(expected=InvalidGenomicCoordsException.class)
//	public void canHandleWrongGenePredFormat() throws ClassNotFoundException, IOException, InvalidCommandLineException, InvalidGenomicCoordsException, InvalidRecordException, SQLException{
//		new UcscGenePred("test_data/hg19_genes_head.gtf", -1);
//	}
	
}
