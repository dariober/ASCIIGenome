#!/bin/bash

# These tests are far from comprehensive since they don't check the interactive session
## MEMO: use `less -R -S` to view colours on terminal with less

# Path to jar and data
# ====================
stvExe=~/Tritume/SamTextViewer.jar
cd ~/svn_git/Java-cafe/trunk/SamTextViewer/test_data/

## Get and prepare chr7.fa file, if not already available
if [ ! -e chr7.fa ]
    then
    wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr7.fa.gz &&
    gunzip chr7.fa.gz
    samtools faidx chr7.fa
fi

# Start tests
# ===========

echo "CAN LOAD BAM FILES"
java -Xmx500m -jar $stvExe -r chr7:5598650-5601530 ds051.actb.bam ear045.oxBS.actb.bam -ni
java -Xmx500m -jar $stvExe -r chr7:5598650-5601530 -fa chr7.fa ds051.actb.bam ear045.oxBS.actb.bam -ni 
java -Xmx500m -jar $stvExe -rpm -r chr7:5598650-5601530 ds051.actb.bam ear045.oxBS.actb.bam -ni 
java -Xmx500m -jar $stvExe ds051.actb.bam -r chr7:5566860 -m 10 -f 16 -ni 

echo "CAN SHOW BS DATA"
java -Xmx500m -jar $stvExe -bs -r chr7:5600000-5600179 -fa chr7.fa ds051.actb.bam ear045.oxBS.actb.bam -ni 

echo "HANDLE NO READS IN INPUT"
java -Xmx500m -jar $stvExe ds051.actb.bam -r chr7:5566860 -m 10 -f 16 -F 16 -ni 

echo "BED FILES"
java -Xmx500m -jar $stvExe refSeq.hg19.short.bed -ni
java -Xmx500m -jar $stvExe refSeq.hg19.short.sort.bed.gz -ni

echo "FROM URL"
java -Xmx500m -jar $stvExe http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeSydhTfbs/wgEncodeSydhTfbsGm12878P300bStdPk.narrowPeak.gz -ni

echo "TABIX FILES"
java -Xmx500m -jar $stvExe test.bedGraph.gz -ni

echo "BIGWIG FROM URL"
java -Xmx500m -jar $stvExe -r chr7:5494331-5505851 http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeSydhTfbs/wgEncodeSydhTfbsGm12878Nrf1IggmusSig.bigWig -ni

echo "TDF"
java -Xmx500m -jar $stvExe -r chr7:1-2149128 hg18_var_sample.wig.v2.1.30.tdf -ni

echo "GTF TABIX"
java -Xmx500m -jar $stvExe -r chr7:1-2149128 hg19.gencode_genes_v19.gtf.gz -ni

#
if [ ! Leishmania_major.ASM272v2.31.dna.genome.fa.fai ]
    then
    wget ftp://ftp.ensemblgenomes.org/pub/release-31/protists/fasta/leishmania_major/dna/Leishmania_major.ASM272v2.31.dna.genome.fa.gz &&
    gunzip Leishmania_major.ASM272v2.31.dna.genome.fa.gz &&
    samtools faidx Leishmania_major.ASM272v2.31.dna.genome.fa
fi
java -Xmx500m -jar $stvExe -r 36:1-2682151 -g Leishmania_major.ASM272v2.31.dna.genome.fa.fai ftp://ftp.ensemblgenomes.org/pub/release-31/protists/gtf/leishmania_major/Leishmania_major.ASM272v2.31.gtf.gz

echo -e "\n\nDONE\n\n"

echo -e "PRINT HELP"
java -Xmx500m -jar $stvExe -h &&


exit
