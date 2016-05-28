package tracks;

import static org.junit.Assert.*;

import java.io.File;

import org.junit.Test;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

public class IntervalFeatureVCFTest {

	@Test
	public void canConstructFromVariantContext() {
	
		VCFFileReader reader= new VCFFileReader(new File("test_data/CHD.exon.2010_03.sites.vcf.gz"));
		VariantContext varCtx = reader.iterator().next();
		IntervalFeatureVCF vcf= new IntervalFeatureVCF(varCtx);
		
		assertEquals("1", vcf.getChrom());
		assertEquals(1105468, vcf.getFrom());
		assertEquals(1105468, vcf.getTo());
		
		reader.close();

		
	}

}
