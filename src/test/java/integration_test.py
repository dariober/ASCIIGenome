#!/usr/bin/env python3

import shutil
import os
import sys
import subprocess
import unittest

ASCIIGenome = "../../../build/libs/ASCIIGenome --debug 2 -ni"


class shell:
    def __init__(self, cmd):
        p = subprocess.Popen(
            cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        stdout, stderr = p.communicate()
        self.returncode = p.returncode
        self.stdout = stdout.decode()
        self.stderr = stderr.decode()
        self.cmd = cmd


class TestCLI(unittest.TestCase):
    def setUp(self):
        sys.stderr.write("\n" + self.id().split(".")[-1] + "\n")  # Print test name
        if os.path.exists("test_out"):
            shutil.rmtree("test_out")
        os.mkdir("test_out")

    def tearDown(self):
        if os.path.exists("test_out"):
            shutil.rmtree("test_out")

    def testCanPrintVersion(self):
        cmd = f"{ASCIIGenome} -v"
        print(cmd)
        p = shell(cmd)
        self.assertEqual(p.returncode, 0)
        self.assertTrue(len(p.stdout) > 4)

    def testSetColorForFeatures(self):
        cmd = f"{ASCIIGenome} ../../../test_data/hg19_genes_head.gtf -x 'goto chr1:6267-17659 && featureColor -r DDX11L1 red -r WASH7P blue'"
        print(cmd)
        p = shell(cmd)
        self.assertEqual(p.returncode, 0)
        self.assertTrue(";9;" in p.stdout)  # 9 is int for red

    def testSetColorForFeaturesWithAwkScript(self):
        cmd = f"""{ASCIIGenome} ../../../test_data/hg19_genes_head.gtf -x "goto chr1:6267-17659 && featureColor -r '\$4 > 13000' red"
                  """
        print(cmd)
        p = shell(cmd)
        self.assertEqual(p.returncode, 0)
        self.assertTrue(";9;" in p.stdout)  # 9 is int for red

    def testCanSetConfigFromFile(self):
        cmd = f"""{ASCIIGenome} -c ../../main/resources/config/white_on_black.conf | grep -F '[48;5;0m'
                  """
        print(cmd)
        p = shell(cmd)
        self.assertEqual(p.returncode, 0)

    def testUsefulMessageIfCramWithoutGenome(self):
        cmd = f"""{ASCIIGenome} ../../../test_data/ds051.actb.cram
                  """
        print(cmd)
        p = shell(cmd)
        self.assertTrue(p.returncode != 0)
        self.assertTrue("CRAM input requires a genome file" in p.stderr)

    def testCanExplainSamFlags(self):
        cmd = f"""{ASCIIGenome} -nf -x 'explainSamFlag 2690'
                  """
        print(cmd)
        p = shell(cmd)
        self.assertEqual(p.returncode, 0)
        self.assertTrue("X    supplementary alignment [2048]" in p.stdout)
        self.assertTrue(".    read is PCR or optical duplicate [1024]" in p.stdout)

    def testCanAvoidHistoryPositionsNotCompatibleWithCurrentGenome(self):
        cmd = f"""{ASCIIGenome} -x 'goto FOOBAR && open ../../../test_data/pairs.sam && p'
                  """
        print(cmd)
        p = shell(cmd)
        self.assertEqual(p.returncode, 0)

    def testFindCaseInsensitive(self):
        cmd = f"""{ASCIIGenome} -nf -x 'print && find .actb' ../../../test_data/hg19_genes.gtf.gz
                  """
        print(cmd)
        p = shell(cmd)
        self.assertEqual(p.returncode, 0)
        self.assertTrue("LACTB" in p.stdout)

    def testFindLiteralMatch(self):
        cmd = f"""{ASCIIGenome} -nf -x 'print && find -F .ACTB' ../../../test_data/hg19_genes.gtf.gz
                  """
        print(cmd)
        p = shell(cmd)
        self.assertEqual(p.returncode, 0)
        self.assertTrue("LACTB" not in p.stdout)

    def testFindLiteralPlus(self):
        cmd = f"""{ASCIIGenome} -nf -x 'goto chr1:1 && print && find -F +' ../../../test_data/hg19_genes.gtf.gz
                  """
        print(cmd)
        p = shell(cmd)
        self.assertEqual(p.returncode, 0)
        self.assertTrue("DDX11L1" in p.stdout)

    def testFindCaseSensitive(self):
        cmd = f"""{ASCIIGenome} -nf -x 'print && find -c actb' ../../../test_data/hg19_genes.gtf.gz
                  """
        print(cmd)
        p = shell(cmd)
        self.assertEqual(p.returncode, 0)
        self.assertTrue("ACTB" not in p.stdout)

    def testGrepCaseInsensitive(self):
        cmd = f"""{ASCIIGenome} -nf -x 'goto chr1:11874 && print && grep -i ddx\d+' ../../../test_data/hg19_genes.gtf.gz
                  """
        print(cmd)
        p = shell(cmd)
        self.assertEqual(p.returncode, 0)
        self.assertTrue("DDX11L1" in p.stdout)

    def testCanReadVcfWithSequenceDictionaryButWithNoRecords(self):
        cmd = f"""{ASCIIGenome} ../../../test_data/norecords.vcf
                  """
        print(cmd)
        p = shell(cmd)
        self.assertEqual(p.returncode, 0)
        self.assertTrue("test_data/norecords.vcf#" in p.stdout)

    def testCanShowReadPairs(self):
        cmd = f"""{ASCIIGenome} -nf -x 'readsAsPairs' ../../../test_data/pairs.sam
                  """
        print(cmd)
        p = shell(cmd)
        self.assertEqual(p.returncode, 0)
        self.assertTrue("G~" in p.stdout and "~g" in p.stdout)

    def testCanInitializeFromVcf(self):
        cmd = f"""{ASCIIGenome} -nf ../../../test_data/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz -x 'goto MT:1234'
                  """
        print(cmd)
        p = shell(cmd)
        self.assertEqual(p.returncode, 0)
        self.assertTrue("MT:1234-" in p.stdout)

        # Because we asked for an invalid chromosome, ASCIIGenome exits with non-zero code.
        cmd = f"""{ASCIIGenome} -nf ../../../test_data/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz -x 'goto FOO'
                  """
        print(cmd)
        p = shell(cmd)
        self.assertEqual(p.returncode, 1)
        self.assertTrue(
            "Cannot find chromosome 'FOO' in sequence dictionary" in p.stderr
        )

        # `show genome` prints to stderr (is this desirable?)
        cmd = f"""{ASCIIGenome} -nf https://raw.githubusercontent.com/dariober/ASCIIGenome/master/test_data/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz -x 'show genome'
                  """
        print(cmd)
        p = shell(cmd)
        self.assertEqual(p.returncode, 0)
        self.assertTrue("MT 16569" in p.stderr)

    def testGenotypeMatrix(self):
        cmd = f"""{ASCIIGenome} -nf ../../../test_data/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz -x 'goto 1:1117997-1204429 && genotype -s HG00096'
                  """
        print(cmd)
        p = shell(cmd)
        self.assertEqual(p.returncode, 0)
        self.assertTrue("HG00096" in p.stdout)

    def testAwkWithGetFunc(self):
        # Note the escape mess here
        cmd = f"""{ASCIIGenome} -nf ../../../test_data/ds051.actb.bam -x "goto chr7:5570087-5570291 && awk 'get(\\\"MD\\\") == \\\"73G0G0\\\"'"
        """
        print(cmd)
        p = shell(cmd)
        self.assertEqual(p.returncode, 0)
        self.assertTrue("Reads: 1;" in p.stdout)

    def testAwkWithVcfFunc(self):
        cmd = f"""{ASCIIGenome} -nf ../../../test_data/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz -x "goto 1:200000-1000000 && awk 'get(\\\"AC\\\") == 35'"
                  """
        print(cmd)
        p = shell(cmd)
        self.assertEqual(p.returncode, 0)
        self.assertTrue("N: 3;" in p.stdout)

        # Keep where 3rd sample has GT == 1|1
        cmd = f"""{ASCIIGenome} -nf ../../../test_data/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz -x "goto 1:200000-1000000 && awk 'get(\\\"GT\\\", 3) == \\\"1|1\\\"'"
                  """
        print(cmd)
        p = shell(cmd)
        self.assertEqual(p.returncode, 0)
        self.assertTrue("N: 3;" in p.stdout)

    def testAwkWithGfxFunc(self):
        cmd = f"""{ASCIIGenome} ../../../test_data/hg19_genes.gtf.gz -x "awk 'get(\\\"gene_name\\\") ~ \\\"WASH\\\"' && goto chr1:10900-30900 && print -full"
                  """
        print(cmd)
        p = shell(cmd)
        self.assertEqual(p.returncode, 0)
        self.assertTrue("WASH7P" in p.stdout)
        self.assertTrue("DDX11L1" not in p.stdout)

        cmd = f"""{ASCIIGenome} ../../../test_data/Homo_sapiens.GRCh38.86.ENST00000331789.gff3 -x "awk 'get(\\\"ID\\\") ~ \\\"ENST\\\"' && goto 7:5525240-5531960 && print -full"
                  """
        print(cmd)
        p = shell(cmd)
        self.assertEqual(p.returncode, 0)
        self.assertTrue("ENST00000331789" in p.stdout)
        self.assertTrue("ENSP" not in p.stdout)

    def testAwkHeaderNames(self):
        # We enclose -x in single quotes so we don't need to escape $ in $POS
        cmd = f"""{ASCIIGenome} -nf ../../../test_data/ds051.actb.bam -x 'goto chr7:5570087-5570291 && awk "$POS > 5570000" actb.bam'
                  """
        print(cmd)
        p = shell(cmd)
        self.assertEqual(p.returncode, 0)
        self.assertTrue("Reads: 8;" in p.stdout)

    def testCanShowHideTrackSettings(self):
        cmd = f"""{ASCIIGenome} -nf ../../../test_data/ds051.actb.bam -x 'goto chr7:5568803-5568975 && show genome'
                  """
        print(cmd)
        p = shell(cmd)
        self.assertEqual(p.returncode, 0)
        self.assertTrue("chrM  16571" in p.stderr)

    def testCanResetGlobalOptions(self):
        cmd = f"""{ASCIIGenome} -nf ../../../test_data/ds051.actb.bam -x 'goto chr7:5568803-5568975 && dropTracks # && setConfig max_reads_in_stack 10000 && zo'
                  """
        print(cmd)
        p = shell(cmd)
        self.assertEqual(p.returncode, 0)
        before = p.stdout

        cmd = f"""{ASCIIGenome} -nf ../../../test_data/ds051.actb.bam -x 'goto chr7:5568803-5568975 && dropTracks # && setConfig max_reads_in_stack 3 && zo'
                  """
        print(cmd)
        p = shell(cmd)
        self.assertEqual(p.returncode, 0)
        after = p.stdout
        self.assertTrue(len(after) < len(before))

    def testCanHandleCoordsOutsideChromLimits(self):
        cmd = f"""{ASCIIGenome} -nf ../../../test_data/ds051.actb.bam -x 'goto chrM && zo 25 && :chrM:1-1000000'
             """
        print(cmd)
        p = shell(cmd)
        self.assertEqual(p.returncode, 0)
        self.assertTrue("chrM:1-16571" in p.stdout)

    def testCanUsePercentScreenCoords(self):
        cmd = f"""{ASCIIGenome} -nf ../../../test_data/ds051.actb.bam  -x 'goto chrM:1 && 0 .2 && 16555 && .1'
            """
        print(cmd)
        p = shell(cmd)
        self.assertEqual(p.returncode, 0)

    def testCanSetGtfAttributeForFeatureName(self):
        cmd = f"""{ASCIIGenome} -nf ../../../test_data/hg19_genes.gtf.gz -x 'nameForFeatures gene_name'
            """
        print(cmd)
        p = shell(cmd)
        self.assertEqual(p.returncode, 0)
        self.assertTrue("DDX11L1" in p.stdout)

    def testCanReadBgzipFasta(self):
        cmd = f"""{ASCIIGenome} -nf -fa ../../../test_data/chr7.fa.gz -r 'chr7:10000'
            """
        print(cmd)
        p = shell(cmd)
        self.assertEqual(p.returncode, 0)
        self.assertTrue("Nctaaccctaaccctaaccc" in p.stdout)


if __name__ == "__main__":
    unittest.main()
