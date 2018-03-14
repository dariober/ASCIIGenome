## Genome files for SamTextViewer

This directory contains genome files suitable for `SamTextViewer`. Genome files list chromosome names and lengths,
they are tab delimited with columns *chromosome* and *chromosome length*, without header. Additional columns are ignored and
comment lines starting with '#' are allowed'.

These files must be named *genome-id*.genome where *genome-id* is the tag passed to the command line.

To download genome files from UCSC genome browser use

```
mysql -N --user=genome --host=genome-mysql.cse.ucsc.edu -A -e \
    "select chrom, size from hg19.chromInfo"  > hg19.genome
```

Or

```
curl -o - -O http://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes > mm10.genome
```

* [GRCh37_1kg](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/) 
  is the human genome reference used by the 1000 genome project. Chromosomes are
  in Ensembl format (1, 2, ..., MT).

