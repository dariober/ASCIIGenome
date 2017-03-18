#!/bin/bash

ASCIIGenome_jar=$1
test_data=$2

function checkZeroExit {
    if [[ $? == 0 ]]
        then
        echo -e "\033[38;5;2m PASSED \033[0m"
    else
        echo -e "\033[38;5;9m FAILED \033[0m"
        exit 1
    fi
}

if [[ -z $1 || -z $2 ]]
    then
    echo -e "DESCRIPTION
Test runner for ASCIIGenome. Beware these tests are not comprehensive and they 
just check that ASCIIGenome exits clean after executing them. 
Output files are written to current working directory. Make sure this is writable
and it can be potentially become littered.

USAGE
run_test.sh /path/to/ASCIIGenome.jar /path/to/test_data/"
    exit 1
fi

echo -e "\033[38;5;2m CAN SHOW HELP \033[0m"
java -Xmx2g -jar $ASCIIGenome_jar -h

echo -e "\033[38;5;2m CAN SHOW VERSION \033[0m"
java -Xmx2g -jar $ASCIIGenome_jar -v 

clear
echo -e "\033[38;5;2m CAN EXECUTE sys COMMAND \033[0m"
java -Xmx2g -jar $ASCIIGenome_jar -ni -x 'sys ls' &&
checkZeroExit

clear
echo -e "\033[38;5;2m CAN SET CONFIG \033[0m"
java -Xmx2g -jar $ASCIIGenome_jar -ni -x 'setConfig metal' &&
checkZeroExit


exit 0

echo "CAN SET CONFIGURATION"
java -jar 

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
java -Xmx500m -jar $stvExe ds051.actb.bam -r chr7:5566860 -x 'samtools -q 10 -F 16' -ni 

echo "CAN SHOW BS DATA"
java -Xmx500m -jar $stvExe -r chr7:5600000-5600179 -fa chr7.fa ds051.actb.bam ear045.oxBS.actb.bam -x 'samtools -q 10 && BSseq' -ni 

echo "HANDLE NO READS IN INPUT"
java -Xmx500m -jar $stvExe ds051.actb.bam -r chr7:5566860 -x ' -f 16 && -F 16' -ni # Note space between quote and -f

echo "BED FILES"
java -Xmx500m -jar $stvExe refSeq.hg19.short.bed -ni
java -Xmx500m -jar $stvExe refSeq.hg19.short.sort.bed.gz -ni

echo "FROM URL"
java -Xmx500m -jar $stvExe http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeSydhTfbs/wgEncodeSydhTfbsGm12878P300bStdPk.narrowPeak.gz -ni

#echo "FROM UCSC"
#java -Xmx500m -jar $stvExe dm6:refGene -ni

echo "BIGWIG FROM URL"
java -Xmx500m -jar $stvExe -r chr7:5494331-5505851 http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeSydhTfbs/wgEncodeSydhTfbsGm12878Nrf1IggmusSig.bigWig -ni

echo "TDF"
java -Xmx500m -jar $stvExe -r chr7:1-2149128 hg18_var_sample.wig.v2.1.30.tdf -ni

echo "GTF TABIX"
java -Xmx500m -jar $stvExe -r chr7:1-2149128 hg19.gencode_genes_v19.gtf.gz -ni

echo "FIND REGEX"
java -Xmx500m -jar $stvExe -r chr7:5772765-5772967 -fa chr7.fa -x 'seqRegex (?i)ac..tg' -ni 

echo "GRACEFULLY HANDLE INVALID INPUT"
java -Xmx500m -jar $stvExe refSeq.hg19.short.bed -x 'foo' -ni
java -Xmx500m -jar $stvExe refSeq.hg19.short.bed -x 'ylim 0 10 *' -ni
java -Xmx500m -jar $stvExe foo.bed -ni
java -Xmx500m -jar $stvExe invalid-1.bedgraph -ni

echo "BATCH FILE"
java -Xmx500m -jar $stvExe -b batch_actb.bed -x 'zo 3 && save >> deleteme.pdf' -g hg19 ear045.oxBS.actb.tdf hg19.gencode_genes_v19.gtf.gz batch_actb.bed > /dev/null
rm deleteme*

echo "SAVING AND LOADING SESSION"
java -Xmx500m -jar $stvExe -x 'sessionSave deleteme.txt' -ni -g hg19 ear045.oxBS.actb.tdf hg19.gencode_genes_v19.gtf.gz batch_actb.bed > /dev/null
java -Xmx500m -jar $stvExe -x deleteme.txt
rm deleteme.txt

echo "TESTING PNG COLOURS"
## Prepare a command string that will change the name and color of the 
## loaded track and save the output as png
cmd=""
for x in red light_red green light_green blue light_blue magenta light_magenta cyan light_cyan yellow light_yellow grey light_grey white black
do
    cmd="${cmd} editNames ^.*$ ${x} && colorTrack $x && save tmp.${x}.png &&"
done
echo "$cmd" > /tmp/cmd.txt

ASCIIGenome -ni -r chr7:5520653-5599501 -x /tmp/cmd.txt ear045.oxBS.actb.tdf
convert -append tmp.*.png /tmp/cols.png
rm tmp.*.png
open /tmp/cols.png # This is a gallery of colours.
rm /tmp/cols.png



echo -e "\n\nDONE\n\n"

echo -e "PRINT HELP"
java -Xmx500m -jar $stvExe -h

exit
