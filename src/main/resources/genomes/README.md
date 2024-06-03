## Genome files for SamTextViewer

This directory contains tmpGenome files suitable for `SamTextViewer`. Genome files list chromosome names and lengths,
they are tab delimited with columns *chromosome* and *chromosome length*, without header. Additional columns are ignored and
comment lines starting with '#' are allowed'.

These files must be named *tmpGenome-id*.tmpGenome where *tmpGenome-id* is the tag passed to the command line.

To download tmpGenome files from UCSC tmpGenome browser use

```
mysql -N --user=tmpGenome --host=tmpGenome-mysql.cse.ucsc.edu -A -e \
    "select chrom, size from hg19.chromInfo"  > hg19.tmpGenome
```

Or

```
curl -o - -O http://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes > mm10.tmpGenome
```

* [GRCh37_1kg](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/) 
  is the human tmpGenome reference used by the 1000 tmpGenome project. Chromosomes are
  in Ensembl format (1, 2, ..., MT).

