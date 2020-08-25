#!/bin/bash

set -e
set -uf -o pipefail

docstring="DESCRIPTION
Test commands work smoothly and do not throw unexpected exceptions.
USAGE
Move to the directory containing integration_test.sh and execute it by giving the path to the 
ASCIIGenome bash script:

    cd /path/to/test/
    ./integration_test.sh /path/to/ASCIIGenome
"

#if [ -z ${1+x} ]; 
#    then echo "var is unset"; 
#else 
#    echo "var is set to '$1'"; 
#fi
#exit

if [[ -z "${1+x}" || $1 == "-h" || $1 == "--help" ]]
then
    echo "$docstring"
    exit 1
fi

# ------------------------------------------------------------------------------

source bashTestFunctions.sh

gzip -c -d ../../../test_data/chr7.fa.gz > ../../../test_data/chr7.fa

ASCIIGenome="$1 --debug 2 -ni"

#pprint 'Can highlight pattern'
#$ASCIIGenome ../../../test_data/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz -r 1:1-400000 -x 'print -hl 200000'

pprint 'Can show version'
$ASCIIGenome -v
assertEquals 0 $?

pprint 'Set color for features'
$ASCIIGenome ../../../test_data/hg19_genes_head.gtf -x "goto chr1:6267-17659 && featureColor -r DDX11L1 red -r WASH7P blue" > tmp.txt
grep ';9;' tmp.txt # 9 is int for red
rm tmp.txt

pprint 'Set color for features with awk script'
$ASCIIGenome ../../../test_data/hg19_genes_head.gtf -x "goto chr1:6267-17659 && featureColor -r '\$4 > 13000' red" > tmp.txt
grep ';9;' tmp.txt # 9 is int for red
rm tmp.txt

pprint 'Can set config from file'
$ASCIIGenome -c ../../main/resources/config/white_on_black.conf | grep -F '[48;5;0m'
printf "\033c"
assertEquals 0 $?

pprint 'Can explain sam flags'
$ASCIIGenome -nf -x 'explainSamFlag 2690' > flags.tmp
grep '^ X ' flags.tmp | grep 'supplementary alignment' > /dev/null
assertEquals 0 $?
grep -v ' X ' flags.tmp | grep 'read is PCR or optical duplicate' > /dev/null
assertEquals 0 $?
rm flags.tmp

pprint 'Can avoid history positions not compatible with current genome'
$ASCIIGenome -x 'goto FOOBAR && open ../../../test_data/pairs.sam && p' > /dev/null
assertEquals 0 $?

pprint 'find - CASE INSENSITIVE'
$ASCIIGenome -nf -x 'print && find .actb' ../../../test_data/hg19_genes.gtf.gz | grep LACTB > /dev/null
assertEquals 0 $?

pprint 'find - LITERAL ".ACTB" does not match "LACTB"'
set +e
$ASCIIGenome -nf -x 'print && find -F .ACTB' ../../../test_data/hg19_genes.gtf.gz | grep LACTB 
assertEquals 1 $?
set -e

pprint 'find - LITERAL "+" - alone it would be an invalid regex'
$ASCIIGenome -nf -x 'goto chr1:1 && print && find -F +' ../../../test_data/hg19_genes.gtf.gz | grep DDX11L1 > /dev/null
assertEquals 0 $?

pprint 'find - CASE SENSITIVE'
set +e
$ASCIIGenome -nf -x 'print && find -c actb' ../../../test_data/hg19_genes.gtf.gz | grep ACTB > /dev/null
assertEquals 1 $?
set -e

pprint 'grep - CASE INSENSITIVE'
$ASCIIGenome -nf -x 'goto chr1:11874 && print && grep -i ddx\d+' ../../../test_data/hg19_genes.gtf.gz | grep DDX11L1 > /dev/null
assertEquals 0 $?

pprint 'Can read VCF with sequence dictionary but with no records'
$ASCIIGenome ../../../test_data/norecords.vcf | grep 'test_data/norecords.vcf#' > /dev/null
assertEquals 0 $?

pprint 'Can show read pairs'
$ASCIIGenome -x 'readsAsPairs' ../../../test_data/pairs.sam > /dev/null

pprint 'Init from VCF'
$ASCIIGenome ../../../test_data/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz -x 'show genome' > /dev/null
$ASCIIGenome https://raw.githubusercontent.com/dariober/ASCIIGenome/master/test_data/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz -x 'show genome' > /dev/null

pprint 'Genotype matrix'
$ASCIIGenome ../../../test_data/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz -x 'goto 1:1117997-1204429 && genotype -f {HOM} && genotype -n 1' > /dev/null
$ASCIIGenome ../../../test_data/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz -x 'goto 1:1117997-1204429 && genotype -s HG00096' > /dev/null

pprint 'Test awk with getSamTag()'
$ASCIIGenome ../../../test_data/ds051.actb.bam -x "goto chr7:5570087-5570291 && awk 'getSamTag(\"NM\") > 0'" > /dev/null

pprint 'Test awk with VCF functions'
$ASCIIGenome ../../../test_data/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz -x "goto 1:200000-1000000 && awk 'getInfoTag(\"AC\") > 0'" > /dev/null
$ASCIIGenome ../../../test_data/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz -x "goto 1:200000-1000000 && awk 'getFmtTag(\"GT\") == \"0|1\"'" > /dev/null

pprint 'Test awk with GTF/GFF'
$ASCIIGenome ../../../test_data/hg19_genes_head.gtf -x "awk 'getGtfTag(\"gene_name\") ~ \"DD\"'" > /dev/null
$ASCIIGenome ../../../test_data/hg19_genes_head.gtf -x "awk 'get(\"gene_name\") ~ \"DD\"'" > /dev/null
$ASCIIGenome ../../../test_data/Homo_sapiens.GRCh38.86.ENST00000331789.gff3 -x "awk 'getGffTag(\"ID\") ~ \"ENST\"'" > /dev/null
$ASCIIGenome ../../../test_data/Homo_sapiens.GRCh38.86.ENST00000331789.gff3 -x "awk 'get(\"ID\") ~ \"ENST\"'" > /dev/null

pprint 'Test header names. Note escape on $'
$ASCIIGenome ../../../test_data/ds051.actb.bam -x "goto chr7:5570087-5570291 && bookmark && awk '\$POS > 5500000' actb.bam" > /dev/null
$ASCIIGenome ../../../test_data/ds051.actb.bam -x "goto chr7:5570087-5570291 && bookmark && awk '\$START > 5500000' Book" > /dev/null
$ASCIIGenome ../../../test_data/CEU.exon.2010_06.genotypes.vcf.gz -x "awk 'length(\$REF) == 1 && length(\$ALT) > 1'" > /dev/null
$ASCIIGenome ../../../test_data/Homo_sapiens.GRCh38.86.ENST00000331789.gff3 -x "awk '\$SOURCE ~ \"havana\"'" > /dev/null

pprint 'Can show/hide track settings'
$ASCIIGenome ../../../test_data/ds051.actb.bam -x 'goto chr7:5568803-5568975 && show genome && show genome' > /dev/null

pprint 'Can reset global options'
$ASCIIGenome ../../../test_data/ds051.actb.bam -x 'goto chr7:5568803-5568975 && setConfig max_reads_in_stack 1000 && zo' > /dev/null

pprint 'Can handle coords outside chrom limits'
$ASCIIGenome ../../../test_data/ds051.actb.bam -x 'goto chrM && zo 25 && :chrM:1-1000000' > /dev/null

$ASCIIGenome ../../../test_data/ds051.actb.bam -x 'goto chr7:5568803-5568975 && zo && zi' > /dev/null

$ASCIIGenome ../../../test_data/ds051.actb.bam -fa ../../../test_data/chr7.fa -x 'goto chr7:5568803-5568975 && zo && zi' > /dev/null

pprint 'Use of PCT screen coords'
$ASCIIGenome ../../../test_data/ds051.actb.bam  -x 'goto chrM:1 && 0 .2 && 16555 && .1' > /dev/null

pprint 'Set gtf attribute for feature name'
$ASCIIGenome ../../../test_data/hg19_genes.gtf.gz -x 'nameForFeatures gene_name' > /dev/null

set +x
echo -e "\033[32mDONE\033[0m"
