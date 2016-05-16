<!-- MarkdownTOC -->

- [Text Only Genome Viewer!](#text-only-genome-viewer)
  - [Key Features](#key-features)
- [Usage](#usage)
  - [Quick start](#quick-start)
  - [Navigation](#navigation)
  - [Search](#search)
  - [Display](#display)
  - [Alignments](#alignments)
- [Supported input](#supported-input)
- [Tips gotchas and miscellanea](#tips-gotchas-and-miscellanea)
  - [Regular expressions](#regular-expressions)
- [Requirements and Installation](#requirements-and-installation)
  - [Installation quick start.](#installation-quick-start)
  - [A little more detail](#a-little-more-detail)
- [Credits](#credits)
- [DEPRECATED](#deprecated)
- [ASCIIGenome: A text based no GUI genome viewer](#asciigenome-a-text-based-no-gui-genome-viewer)
- [Usage](#usage-1)
  - [Quick start](#quick-start-1)
  - [Moving around the genome](#moving-around-the-genome)
  - [Display options](#display-options)
  - [Filtering reads](#filtering-reads)
  - [Searching features in annotation files](#searching-features-in-annotation-files)
  - [Genome option](#genome-option)
  - [Formatting of reads and features](#formatting-of-reads-and-features)

<!-- /MarkdownTOC -->

<!-- MEMO: Compile, package and uplaod to github releases
- Write-out jar from Eclipse
/Users/berald01/Tritume/ASCIIGenome.jar
 -->


Text Only Genome Viewer!
========================

```ASCIIGenome``` is a command-line genome browser running from terminal window and solely based on ASCII characters.
Since ```ASCIIGenome``` does not require a graphical interface it is particularly useful for 
quickly visualizing genomic data on remote servers. With some imagination ```ASCIIGenome``` is the Vim 
of genome viewers.

The closest program to ```ASCIIGenome``` is [samtools tview](http://samtools.sourceforge.net/tview.shtml) but 
```ASCIIGenome``` offers much more flexibility, similar to popular GUI viewers like [IGV](https://www.broadinstitute.org/igv/).


<img src="screenshots/ex3.png" width="600">

Key Features
------------

* Command line input and interaction, no graphical interface, minimal [installation and requirements](#requirements-and-installation)
* Can load multiple files in various [formats](#supported-input)
* Can access remote files view URL or ftp address
* Easy [navigation](#moving-around-the-genome), [search](#searching-features-in-annotation-files), and filtering options
* Support for BS-Seq alignment

Usage
=====

Quick start
-----------

* Minimal example, open and browse a bam file:

```
ASCIIGenome aln.bam
```

* Open several files, including the fasta reference genome and a remote file (from [ENCODE](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeSydhTfbs/)):

```
ASCIIGenome -fa genome.fa aln.bam scores.bigWig genes.gtf \
  http://.../wgEncodeSydhTfbsGm12878P300sc584IggmusPk.narrowPeak.gz
```

* Navigating, finding & filtering stuff

Once started, ```ASCIIGenome``` makes it easy to browse the genome. The picture below shows the distribution of transcripts on chromosome 36 of *Leishmania major*. It is clearly visible how transcripts in Leishmania tend to be grouped in blocks transcribed from the same direction (blue: forward strand, pink: reverse strand). Note how overlapping features are stacked on top of each other.

<img src="screenshots/leishmania_transcripts.png" width="700">

This screenshot has been produced by first loading the *L. major* GTF file:

```
ASCIIGenome ftp://ftp.ensemblgenomes.org/pub/release-31/protists/gtf/leishmania_major/Leishmania_major.ASM272v2.31.gtf.gz
```

Then move to visualize chromosome 36: `goto 36:1-2682151`, then make visible only the 'transcript' features: `visible '.*\ttranscript\t.*'`

Navigation
----------

```
f / b 
      Small step forward/backward 1/10 window
ff / bb
      Large step forward/backward 1/2 window
zi / zo [x]
      Zoom in / zoom out x times (default x= 1). Each zoom halves or doubles the window size
goto chrom:from-to
      Go to given region. E.g. "goto chr1:1-1000" or chr1:10 or chr1. goto keyword can be replaced with ':' (like goto in vim)
<from> [to]
      Go to position <from> or to region "from to" on current chromosome. E.g. 10 or "10 1000" or "10-1000"
+/-<int>[k,m]
      Move forward/backward by <int> bases. Suffixes k (kilo) and M (mega) allowed. E.g. -2m or +10k
p / n
      Go to previous/next visited position
```

Search
------

```
next <trackId>
      Move to the next feature in <trackId> on *current* chromosome
find <regex> [trackId]
      Find the next record in trackId matching regex. Use single quotes for strings containing spaces.
      For case insensitive matching prepend (?i) to regex. E.g. "next '(?i).*actb.*' myTrack#1"
```

Display
-------

```
visible [show regex] [hide regex] [track regex]
      In annotation tracks, only include rows captured by [show regex] and exclude [hide regex].
      Apply to tracks captured by [track regex]. With no optional arguments reset to default: "'.*' '^$' '.*'"
      Use '.*' to match everything and '^$' to hide nothing. E.g. "visible .*exon.* .*CDS.* .*gtf#.*"
trackHeight <int> [track regex]
      Set track height to int lines for all tracks captured by regex. Default regex: '.*'
ylim <min> <max> [track regex]
      Set limits of y axis for all track IDs captured by regex. Use na to autoscale to min and/or max.
      E.g. ylim 0 na. If regex is omitted all tracks will be captured. Default: "ylim na na .*"
dataCol <idx> [regex]
      Select data column for all bedgraph tracks captured by regex. <idx>: 1-based column index.
print / printFull
      Turn on/off the printing of bed/gtf features.
      print clip lines to fit the screen, printFull will wrap the long lines
showGenome
      Print the genome file
addTracks [file/url]...
      Add tracks
history
      Show visited positions
```


Alignments
----------

```
-f
    Required sam flags. Use 4096 for reads on top strand
-F
    Filtering sam flags. Use 4096 for reads on top strand
-q
    Minumum mapping quality for a read to be considered
-m
    Maximum number of lines to print for read tracks.
-rpm
    Toggle on/off the normalization of Reads Per Million for bam input. Default off
-ml
    Maximum number of lines to print for each methylation track
```

Supported input
===============

File name extensions matter as file types are usually recognized by their extension in case insensitive mode.

* **bam** files should be sorted and indexed, e.g. with `samtools sort` and `samtools index`. 
  Paths to remote URLs are supported but painfully slow.
* **bigWig** recognized by extension `.bw` or `.bigWig`. Remote URLs supported.
* **bedGraph** recognized by extension `.bedGraph` or `.bedgraph`
* **bed**, **gtf**, **gff** recognized by respective extensions. Remote URLs supported. 
* **tdf** This is very useful for quickly displaying very large intervals like tens of megabases or entire chromosomes see [tdf](https://www.broadinstitute.org/igv/TDF)
* All other extensions (e.g. txt, narrowPeak) will be treated as bed files, provided the format is actually bed!

Notable formats currently **not** supported:  cram (because it doesn't seem to be very popular), bigBed (for local files I think tabix indexing is preferable), vcf (too bad it's not recognized, it's because I don't use vcf much but supporting them should be reasonably easy). 

All plain text formats (bed, bedgraph, etc) can be read as gzipped and there is no need to decompress them.

Bedgraph files should be sorted by position, a `sort -k1,1 -k2,2n` will do. Unindexed bedGraph files are first bgzipped and indexed to temporary files which are deleted on exit. This can take time for large files so consider creating the index once for all with [tabix](http://www.htslib.org/doc/tabix.html), *e.g.*

```
bgzip my.bedgraph &&
tabix -p bed my.bedgraph.gz
```

Bed & gtf file are not required to be sorted or index but in this case they are loaded in memory. To save memory and time for large files you can again index them as above. Loaded in memory is typically fast for files of up to ~1/2 million records.

For input format specs see also [UCSC format](https://genome.ucsc.edu/FAQ/FAQformat.html) and for guidelines on the choice of format see [IGV recommendations](https://www.broadinstitute.org/igv/RecommendedFileFormats).


Tips gotchas and miscellanea
============================

* Use UP and DOWN arrow keys to display the previous/next executed command just like on Unix terminal. 
 Commands can be auto-completed with TAB, like Unix terminal but more rudimentary really. 

* On interactive prompt, no input, *i.e.* just pressing ENTER, repeats the previous command (useful for quick navigation).

* The `next` command does exactly that, it moves to the next feature. If there are no more features after the current position it
doesn't rewind to the beginning (use `:1` for that) and it doesn't move to another chromosome (use `-r chrom`). 

* **Performance** Alignment files are typically accessed very quickly but `ASCIIGenome` becomes slow when the window size grows
above a few hundreds of kilobases. Annotation files (bed, gff, gtf) are loaded in memory unless they are indexed with `tabix`. 


#### Regular expressions

*A programmer has a problem and thinks* "I know, I'll use regular expressions". *Now he has two problems.*

The options that take regular expressions assume some familiarity with regexes and they can be tricky at first. In particular remember that lines are scanned as a single string for matches. These are some pointers:

* Almost always you want to surround your pattern with `.*`, e.g. to find the ACTB gene use `.*ACTB.*`, but note that this will hit LACTB as well unless you are more specific (e.g. `.*"ACTB".*`) or you use fancier regex. E.g. `.*[^A-Z0-9]ACTB[^A-Z0-9].*` will match ACTB but not LACTB or ACTB9. 

* Just using 'ACTB' as pattern, as you would do with `grep`, will return nothing as no line should contain _only_ 'ACTB'.

* Use the `(?i)` modifier to match in case insensitve mode, e.g. '(?i).*actb.*'

* For reference, regexes are parsed by Java `String.match()` method.

Requirements and Installation
=============================

## Installation quick start. 

In the commands below replace version number with the latest from [releases](https://github.com/dariober/Java-cafe/releases):

```
wget https://github.com/dariober/Java-cafe/releases/download/v0.1.0/ASCIIGenome-0.1.0.zip
unzip ASCIIGenome-0.1.0.zip
cd ASCIIGenome-0.1.0/
chmod a+x ASCIIGenome
cp ASCIIGenome.jar /usr/local/bin/ # Or ~/bin/
cp ASCIIGenome /usr/local/bin/     # Or ~/bin/ 
```

A little more detail
--------------------

`ASCIIGenome.jar` requires **Java 1.7+** and this is (should be) the only requirement. There is virtually no installation needed as `ASCIIGenome` is pure Java and should work on most (all?) platforms. Download the zip file `ASCIIGenome-x.x.x.zip` from [releases](https://github.com/dariober/Java-cafe/releases), unzip it and execute the jar file with

    java -jar /path/to/ASCIIGenome.jar --help

To avoid typing `java -jar ...` every time, you can put both the helper 
script `ASCIIGenome` and the jar file ```ASCIIGenome.jar``` in the same directory in your `PATH` and execute with:

  ASCIIGenome [options]

Note the helper is a bash script. To set the amount of memory available to java use the `-Xmx` option as e.g. `java -Xmx1500m -jar ...`.

If for some reason the text formatting misbehaves, disable it with the `-nf` option. 

Credits
=======

* Bam processing is mostly done with the [samtools/htsjdk](https://github.com/samtools/htsjdk) library.
* Bigwig and tdf are processed with classes from [IGV](https://github.com/igvteam/igv) source code.
* Block compression and indexing done using [jvarkit](https://github.com/lindenb/jvarkit)


DEPRECATED
==========

ASCIIGenome: A text based no GUI genome viewer
================================================

```ASCIIGenome``` is a command line genome viewer useful to browse and visualize sequence alignment and annotation files
on **console screen**. It aims to be similar to ```samtools tview``` but with the flexibility of
GUI genome viewers like IGV.

Features that attempt to combine text based viewers (```tview```) with GUI viewers:

* Command line input and interaction, no graphical interface.
* Can load multiple files in various formats.
* Support for BS-Seq alignment files.
* Navigation and search options

![ex3](https://github.com/dariober/Java-cafe/blob/master/ASCIIGenome/screenshots/ex3.png)

Usage
=====

Quick start
-----------

These are some examples, for brevity using the helper `ASCIIGenome`

Display a bam file together with a gtf annotation file, go straight to position chr7:5566640-5569055 (atcb gene). This is RNA-Seq data:

    ASCIIGenome -r chr7:5566640-5569055 ds051.actb.bam hg19.gencode_genes_v19.gtf.gz

![ex1](https://github.com/dariober/Java-cafe/blob/master/ASCIIGenome/screenshots/ex1.png)


The header line:

```
ds051.actb.bam; Each . = 44.68; max depth: 893.62x; 
```

gives the name of the file, scale of the read depth track (in this example one dot corresponds to 44.68 reads) and the maximum read depth
in the current view (893.62).

The line at the bottom of the tracks

```
chr7:5567688-5567847; 160 bp; 1.0 bp/char; Filters: -q 0 -f 0 -F 4; Mem: 552 MB;
```

shows the current position, the width and scale of the view, filters applied to the bam files and the memory usage.

For visualizing BS-Seq data add the `-bs` flag and provide a reference fasta file:

    ASCIIGenome -bs fa chr7.fa ds051.actb.bam 

<img src="screenshots/exBSmode.png" width="450">

<img src="screenshots/exBSmode-2.png" width="450">


After starting `ASCIIGenome` you can navigate the genome with the following interactive commands. 

<img src="screenshots/bedCluster.png" width="450">

## Moving around the genome

Use the following keys and commands to navigate the genome

```
f / b 
      Small step forward/backward 1/10 window
ff / bb
      Large step forward/backward 1/2 window
zi / zo
      Zoom in / zoom out
p / n
      Go to previous/next visited position
<from>-[to]
      Go to position <from> or to region <from>-[to] on current chromosome. E.g. '10' or '10-1000'
+/-<int>[k,m]
      Move forward/backward by <int> bases. Suffixes k and m allowed. E.g. -2m or +10k
```

Display options
---------------

These options can be used to set the limits of the y axis and their height.

```
visibile [show regex] [hide regex] [track regex]
        In annotation tracks, only include rows captured by [show regex] and exclude [hide regex].
        Apply to annotation tracks captured by [track regex]. With no optional arguments reset to default: "'.*' '^$' '.*'"
        Use '.*' to match everything and '^$' to hide nothing. Ex "visible .*exon.* .*CDS.* .*gtf#.*"
ylim <min> <max> [regex]
        Set limits of y axis for all track IDs captured by regex. Default regex: '.*'
dataCol <idx> [regex]
        Select data column for all bedgraph tracks captured by regex. <idx>: 1-based column index
```

also see:

```
print
        Turn on/off the printing of bed/gtf features in current interval
rNameOn / rNameOff
        Show/Hide read names
history
        Show visited positions
```

These options apply to bam files:

```
-m
    Maximum number of lines to print for read tracks. No limit If < 0
-rpm
    Toggle on/off the normalization of Reads Per Million for bam input. Default off
-d
    Maximum number of lines to print for coverage tracks
-ml
    Maximum number of lines to print for each methylation track
```

Filtering reads
---------------

Reads in bam files can be filtered in and out using the -f, -F and -q options as in `samtools`:

```
-f
    Required sam flags. Use 4096 for reads on top strand
-F
    Filtering sam flags. Use 4096 for reads on top strand
-q
    Minumum mapping quality for a read to be considered

```

## Searching features in annotation files

These options apply to bed, gff and gtf files:

```
next <trackId>
        Move to the next feature in <trackId> on *current* chromosome
find <regex> [trackId]
        Find the next record in trackId matching regex. Use single quotes for strings containing spaces.
        For case insensitive matching prepend (?i) to regex e.g. '(?i).*actb.*'
```

Note that `find` searches the entire lines for matches to the given regular expression. So to find the 'ACTB' gene name use `find .*ACTB.*`, using just `find ACTB` as you would with `grep` will return no matches as no line in a bed file can match just that. 

## Genome option

An optional genome file can be passed to option `-g/--genome` to give a set of allowed sequences and their sizes so that browsing is constrained to the real genomic space. 
The genome file is also used to represent the position of the current window on the chromosome, which is handy to navigate around.

There are three options to pass a genome file:

* A tag identifying a built-in genome, e.g. hg19. See [genomes](http://github.com/dariober/Java-cafe/ASCIIGenome/resources/genomes) for available genomes

* A local file, tab separated with columns chromosome name and length. See [genomes](http://github.com/dariober/Java-cafe/ASCIIGenome/resources/genomes) for examples.

* A bam file with suitable header.

Note that if the input list of files contains a bam file, the `--genome` option is effectively ignored as the genome dictionary is extracted from the bam header.


## Formatting of reads and features

When aligned reads are show at single base resolution, read bases follow the same convention as samtools: 
Upper case letters and `.` for read align to forward strand, lower case and `,` otherwise; second-in-pair reads are underlined;
grey-shaded reads have mapping quality of <=5. In bisulfite mode the characters M, U, m, u are used for methylated and unmethylated bases on forward and reverse strands.
