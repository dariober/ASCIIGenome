Text Only Genome Viewer!
========================

<!-- MarkdownTOC -->

  - [Description](#description)
  - [Key Features](#key-features)
- [Usage](#usage)
  - [Quick start](#quick-start)
  - [General remarks](#general-remarks)
  - [Navigation](#navigation)
  - [Find](#find)
  - [Display](#display)
  - [Alignments](#alignments)
  - [Genome option](#genome-option)
  - [Formatting of reads and features](#formatting-of-reads-and-features)
- [Supported input](#supported-input)
- [Tips gotchas and miscellanea](#tips-gotchas-and-miscellanea)
- [Requirements and Installation](#requirements-and-installation)
  - [Installation quick start](#installation-quick-start)
  - [A little more detail](#a-little-more-detail)
- [Credits](#credits)
- [TODO](#todo)

<!-- /MarkdownTOC -->

<!-- 
MEMO: Compile, package and upload to github releases
- Write-out jar from Eclipse
cd ~/svn_git/ASCIIGenome/trunk
mkdir ASCIIGenome-0.1.0 # This should match the version in ArgParse
cp ASCIIGenome ASCIIGenome-0.1.0/
cp /Users/berald01/Dropbox/Public/ASCIIGenome.jar ASCIIGenome-0.1.0/
zip -r ASCIIGenome-0.1.0.zip ASCIIGenome-0.1.0
rm -r ASCIIGenome-0.1.0

// Upload ASCIIGenome-0.1.0.zip to github releases and delete

 -->

Description
-----------

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
* Can access remote files via URL or ftp address
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

Then move to visualize chromosome 36: `goto 36:1-2682151`, then make visible only the 'transcript' features: `visible '\ttranscript\t'`

General remarks
---------------

As there is no GUI, everything is handled thorough command line. Once `ASCIIGenome` is started enter
one of the commands below and press ENTER to execute, e.g.:

```
[h] for help: ff <ENTER>
```

will move the window forward by half its size. `h <ENTER>` will show help.

Some features of Unix console are enabled, the arrow keys UP and DOWN scroll previous commands and TAB autocompletes commands.
Just pressing ENTER will repeat the previous command and this is handy to quickly scroll along the genome. For example:

```
[h] for help: ff <ENTER> ## Move forward
[h] for help: <ENTER>    ## Move forward again...
[h] for help: <ENTER>    ## ... and again
```

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
from [to]
      Go to position <from> or to region "from to" on current chromosome. E.g. 10 or "10 1000" or "10-1000"
+/-<int>[k,m]
      Move forward/backward by <int> bases. Suffixes k (kilo) and M (mega) allowed. E.g. -2m or +10k
p / n
      Go to previous/next visited position
next / next_start [trackId]
      Move to the next feature in trackId on *current* chromosome
      'next' centers the window on the found feature while 'next_start' sets the window at the start of the feature.
```

The `next` command does exactly that, it moves to the next feature. If there are no more features after the current position it
doesn't rewind to the beginning (use `1` for that) and it doesn't move to another chromosome (use `got chrom`). 


Find
----

```
find_first <regex> [trackId]
      Find the first (next) record in trackId containing regex. Use single quotes for strings with spaces.
find_all <regex> [trackId]
      Find all records on chromosome containing regex. The search stops at the first chromosome
      returning hits starting with the current one. Useful to get all gtf records of a gene
seqRegex <regex>
      Find regex in reference sequence and show matches as and additional track.
      Useful to show restriction enzyme sites, CpGs etc.
```

Display
-------

```
visible [show regex] [hide regex] [track regex]
      In annotation tracks, only include rows containing [show regex] and exclude [hide regex].
      Apply to tracks containing  [track regex]. With no optional arguments reset to default: "'.*' '^$' '.*'"
      Use '.*' to match everything and '^$' to hide nothing. E.g. "visible exon CDS gtf"
trackHeight <int> [track regex]
      Set track height to int lines for all tracks containing regex. Default regex: '.*'
trackColour <colour> [track regex]
      Set colour for tracks containing regex. All colours except white, black, default
      accept the prefix 'light_'. Available colours:
      red green yellow blue magenta cyan grey white black default
      E.g. trackColour light_blue ts.*gtf
ylim <min> <max> [track regex]
      Set limits of y axis for all track IDs containing regex. Use na to autoscale to min and/or max.
      E.g. ylim 0 na. If regex is omitted all tracks will be captured. Default: "ylim na na .*"
dataCol <idx> [regex]
      Select data column for all bedgraph tracks containing regex. <idx>: 1-based column index.
print     [track regex] 
printFull [track regex] 
      Print the lines of the annotation tracks containing  [track regex]. print clips lines to
      fit the screen. printFull wraps long lines. With no arguments all tracks are printed.
showGenome
      Print the genome file
addTracks [file or url]...
      Add tracks
orderTracks [track#1] [track#2]...
      Reorder tracks
history
      Show visited positions
```


Alignments
----------

```
rpm [track regex]
      Toggle display of read coverage from raw count to reads per million
      for alignment files containing [track regex]
-f INT 
-F INT 
      Include (-f) and exclude (-F) reads with INT bits set
mapq INT
      Include reads with mapq >= INT
maxLines INT
      Maximum number of lines to print for alignment tracks
BSseq [track regex]
      Toggle bisulfite mode for read tracks containing regex. Ignored without reference fasta sequence
```

Genome option
-------------

An optional genome file can be passed to option `-g/--genome` to give a set of allowed sequences and their sizes so that browsing is constrained to the real genomic space. 
The genome file is also used to represent the position of the current window on the chromosome, which is handy to navigate around.

There are three ways to pass a genome file:

* A tag identifying a built-in genome, e.g. hg19. See [genomes](http://github.com/dariober/Java-cafe/ASCIIGenome/resources/genomes) for available genomes

* A local file, tab separated with columns chromosome name and length. See [genomes](http://github.com/dariober/Java-cafe/ASCIIGenome/resources/genomes) for examples.

* A bam file with suitable header.

Note that if the input list of files contains a bam file, the `--genome` option is effectively ignored as the genome dictionary is extracted from the bam header.


Formatting of reads and features
--------------------------------

When aligned reads are show at single base resolution, read bases follow the same convention as samtools: 
Upper case letters and `.` for read align to forward strand, lower case and `,` otherwise; second-in-pair reads are underlined;
grey-shaded reads have mapping quality of <=5. 

In bisulfite mode, the characters M, U, m, u are used for methylated and unmethylated bases on forward and reverse strands. For example:

<img src="screenshots/exBSmode-2.png" width="450">


Supported input
===============

File name extensions matter as file types are usually recognized by their extension in case insensitive mode.

* **bam** files should be sorted and indexed, e.g. with `samtools sort` and `samtools index`. 
  Paths to remote URLs are supported but painfully slow (*IGV seems to suffer of the same issue*).
* **bigWig** recognized by extension `.bw` or `.bigWig`. Remote URLs supported.
* **bedGraph** recognized by extension `.bedGraph` or `.bedgraph`
* **bed**, **gtf**, **gff** recognized by respective extensions. Remote URLs supported. 
* **tdf** This is very useful for quickly displaying very large intervals like tens of megabases or entire chromosomes see [tdf](https://www.broadinstitute.org/igv/TDF)
* **vcf** Supported but not too sophisticated representation. URL should be supported but it appears ftp from 1000genomes doesn't work (same for IGV).
* All other extensions (e.g. txt, narrowPeak) will be treated as bed files, provided the format is actually bed!

Notable formats currently **not** supported:  cram, bigBed.

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

* **Performance** Alignment files are typically accessed very quickly but `ASCIIGenome` becomes slow when the window size grows
above a few hundreds of kilobases. Annotation files (bed, gff, gtf) are loaded in memory unless they are indexed with `tabix`. 

* **Regular expression** Use the `(?i)` modifier to match in case insensitve mode, e.g. '(?i).*actb.*'


Requirements and Installation
=============================

Installation quick start 
------------------------

In the commands below replace version number with the latest from [releases](https://github.com/dariober/Java-cafe/releases):

```
wget https://github.com/dariober/ASCIIGenome/releases/download/v0.1.0/ASCIIGenome-0.1.0.zip
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

TODO
====

* Allow a string of options to be set at prompt. E.g. `-F 16; mapq 10; ylim 0 na`
* Add a `bookmark` command create on the fly an IntervalFeatureTrack with bookmarked regions. 
API: `bookmark`: Add current region to bookmarks; `bookmark show` print current bookmarked features; `bookmark clear`, etc.
* Command to "go to other end" of feature?
* Enable (some) options to be set at start 

