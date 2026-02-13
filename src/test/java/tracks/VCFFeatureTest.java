package tracks;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import colouring.Config;
import com.github.lindenb.jvarkit.variant.bcf.BCFFileReader;
import exceptions.InvalidColourException;
import exceptions.InvalidConfigException;
import exceptions.InvalidGenomicCoordsException;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import org.junit.Before;
import org.junit.Test;
import samTextViewer.Utils;

public class VCFFeatureTest {

  private String ideogramToString(List<FeatureChar> fchars, boolean noFormat)
      throws InvalidColourException {
    StringBuilder sb = new StringBuilder();
    for (FeatureChar x : fchars) {
      sb.append(x.format(noFormat));
    }
    return sb.toString();
  }

  @Before
  public void setConfig() throws IOException, InvalidConfigException {
    new Config(null);
  }

  @Test
  public void handleMissingInfoLinesInHeader() {
    VCFFileReader reader = new VCFFileReader(new File("test_data/malformed_header3.vcf.gz"));
    reader.close();
  }

  @Test
  public void canFormatVCFLine() throws InvalidColourException {

    List<Double> rulerMap = new ArrayList<>();
    for (int i = 1; i < 100; i++) {
      rulerMap.add((double) i);
    }

    // Prepare header
    VCFFileReader reader = new VCFFileReader(new File("test_data/CHD.exon.2010_03.sites.vcf.gz"));
    VCFHeader vcfHeader = reader.getFileHeader();
    reader.close();
    VCFCodec vcfCodec = new VCFCodec();
    vcfCodec.setVCFHeader(vcfHeader, Utils.getVCFHeaderVersion(vcfHeader));

    String vcfLine = "1 10 . C G 23 PASS AA=.,foo;AC=.;AN=.DP=.".replaceAll(" ", "\t");
    VCFFeature vcff = new VCFFeature(vcfLine, vcfCodec);

    vcff.mapToScreen(rulerMap);
    assertEquals(1, vcff.getIdeogram(true, true).size());
    assertEquals('G', vcff.getIdeogram(true, true).get(0).getText());
    assertTrue(
        vcff.getIdeogram(true, true)
            .get(0)
            .format(false)
            .contains("[")); // Just check there is a formatting char

    // Deletion
    vcfLine = "1 10 . CTTG C 23 PASS AA=.;AC=.;AN=.DP=.".replaceAll(" ", "\t");
    vcff = new VCFFeature(vcfLine, vcfCodec);
    vcff.mapToScreen(rulerMap);
    assertEquals('D', vcff.getIdeogram(true, true).get(0).getText());
    assertEquals(3, vcff.getIdeogram(true, true).size());

    // Insertion
    vcfLine = "1 10 . C CTTG 23 PASS AA=.;AC=.;AN=.DP=.".replaceAll(" ", "\t");
    vcff = new VCFFeature(vcfLine, vcfCodec);
    vcff.mapToScreen(rulerMap);
    assertEquals("I", vcff.getIdeogram(true, true).get(0).format(true));
    assertEquals(3, vcff.getIdeogram(true, true).size());

    // Multiple alleles
    vcfLine = "1 10 . C CTTG,A 23 PASS AA=.;AC=.;AN=.DP=.".replaceAll(" ", "\t");
    vcff = new VCFFeature(vcfLine, vcfCodec);
    vcff.mapToScreen(rulerMap);
    assertEquals("|", vcff.getIdeogram(true, true).get(0).format(true));
  }

  @Test
  public void canFormatVCFLineStructVar()
      throws InvalidColourException {

    List<Double> rulerMap = new ArrayList<Double>();
    for (int i = 1; i < 100; i++) {
      rulerMap.add((double) i);
    }

    // Prepare header
    VCFFileReader reader =
        new VCFFileReader(new File("test_data/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz"));
    VCFHeader vcfHeader = reader.getFileHeader();
    reader.close();
    VCFCodec vcfCodec = new VCFCodec();
    vcfCodec.setVCFHeader(vcfHeader, Utils.getVCFHeaderVersion(vcfHeader));

    String vcfLine =
        "1 668630 DUP_delly_DUP20532 G <CN2> . PASS AC=64;AF=0.0127795;AFR_AF=0.0015;AMR_AF=0;AN=5008;CIEND=-150,150;CIPOS=-150,150;CS=DUP_delly;EAS_AF=0.0595;END=850204;EUR_AF=0.001;IMPRECISE;NS=2504;SAS_AF=0.001;SITEPOST=1;SVTYPE=DUP GT 0|0 0|0 0|0"
            .replaceAll(" ", "\t");
    VCFFeature ift = new VCFFeature(vcfLine, vcfCodec);
    ift.mapToScreen(rulerMap);
    assertEquals(850204, ift.getTo());
    assertEquals("|", ift.getIdeogram(true, true).get(0).format(true));
  }
}
