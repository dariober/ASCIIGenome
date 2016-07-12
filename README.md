Text Only Genome Viewer!
========================

<!-- MarkdownTOC -->

- [Description](#description)
- [Key Features](#key-features)
- [Usage](#usage)
  - [Quick start](#quick-start)
- [Interactive commands](#interactive-commands)
  - [Navigation](#navigation)
    - [**f** and **b**](#f-and-b)
    - [**ff** and **bb**](#ff-and-bb)
    - [**zi** *x* and **zo** *x*](#zi-x-and-zo-x)
    - [**goto** chrom:**from-to**](#goto-chromfrom-to)
    - [**from INT** *to INT*](#from-int-to-int)
    - [**+/- INT** *k,m*](#--int-km)
    - [**p** and **n**](#p-and-n)
    - [**next** *trackId* and **next_start** *trackId*](#next-trackid-and-next_start-trackid)
  - [Find](#find)
    - [**find_first** regex *trackId*](#find_first-regex-trackid)
    - [**find_all** regex *trackId*](#find_all-regex-trackid)
    - [**seqRegex** *regex*](#seqregex-regex)
  - [Display](#display)
    - [**visible** *show-regex* *hide-regex* *track-regex*](#visible-show-regex-hide-regex-track-regex)
    - [**trackHeight** INT *track regex*](#trackheight-int-track-regex)
    - [**colorTrack** color *track regex*](#colortrack-color-track-regex)
    - [**ylim** min max *track regex*](#ylim-min-max-track-regex)
    - [**dataCol** idx *regex*](#datacol-idx-regex)
    - [print *track regex*](#print-track-regex)
    - [**printFull** *track regex*](#printfull-track-regex)
    - [**showGenome**](#showgenome)
    - [**addTracks** *files or urls*](#addtracks-files-or-urls)
    - [**orderTracks** *track#1 track#2...*](#ordertracks-track1-track2)
    - [**history**](#history)
  - [Alignments](#alignments)
    - [**rpm** *track regex*](#rpm-track-regex)
    - [**-f** INT and  **-F** INT](#-f-int-and---f-int)
    - [**mapq** INT](#mapq-int)
    - [**maxLines** INT](#maxlines-int)
    - [**BSseq** *track regex*](#bsseq-track-regex)
- [Genome option](#genome-option)
- [Formatting of reads and features](#formatting-of-reads-and-features)
- [Saving screenshots](#saving-screenshots)
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
===========

```ASCIIGenome``` is a command-line genome browser running from terminal window and solely based on ASCII characters.
Since ```ASCIIGenome``` does not require a graphical interface it is particularly useful for 
quickly visualizing genomic data on remote servers. With some imagination ```ASCIIGenome``` is the Vim 
of genome viewers.

The closest program to ```ASCIIGenome``` is [samtools tview](http://samtools.sourceforge.net/tview.shtml) but 
```ASCIIGenome``` offers much more flexibility, similar to popular GUI viewers like [IGV](https://www.broadinstitute.org/igv/).


<img src="screenshots/ex3.png" width="600">

Key Features
============

* Command line input and interaction, no graphical interface, minimal [installation and requirements](#requirements-and-installation)
* Can load multiple files in various [formats](#supported-input)
* Can access remote files via URL or ftp address
* Easy [navigation](#moving-around-the-genome), [search](#searching-features-in-annotation-files) and find features and sequences by regular expression, filtering options
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


Interactive commands
====================

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

In the documentation below, commands are in bold, mandatory arguments are in plain style and 
optional argument are in italics.

Navigation
----------

#### **f** and **b**

Move by forward or backward by 1/10 of a window

#### **ff** and **bb**

Move by forward or backward by 1/2 of a window

#### **zi** *x* and **zo** *x*

Zoom in / zoom out *x* times. Each zoom halves or doubles the window size. To zoom in or out really
quickly use `zi 10` or more. Default x is 1.

#### **goto** chrom:**from-to**

Go to region *chrom:from-to*. *E.g.* `goto chr1:1-1000`. Also recognized is the format `goto:chr1:10` or just
`goto:chr1`. The character ':' is a shortcut for `goto`, *e.g.* `:chr1`.

#### **from INT** *to INT*

Go to position `from INT` or to region `from INT to INT` on **current chromosome**.   For example
`[h] for help: 10` will jump to position 10 and `[h] for help: 10 1000` to region spanned by
position  10 to 1000. Also allowed is the hyphenated format `[h] for help: 10-1000`. If a list of
integers is given, the first and last are taken as *from* and *to*; e.g. `10 250
500 750 1000` is the same as `10 1000`. This is handy to copy and paste intervals from the ruler
above the prompt.

#### **+/- INT** *k,m*

Move forward/backward by INT bases. Suffixes k (kilo) and M (mega) expanded to x1000 and
x1,000,000. *E.g.* `-2m` or `+10k` or `+10.5k`.

#### **p** and **n**

Go to the **p**revious or **n**ext visited position. Similar to the back and forward arrows of an
Internet browser.

#### **next** *trackId* and **next_start** *trackId*

Move to the next feature in trackId on *current* chromosome. `next` centers the window on the found
feature and zooms out. This is useful for quickly browsing through annotation files of genes or
ChIP-Seq peaks in combination with read coverage tracks (bigwig, tdf, etc.). `next_start` instead
sets the window right at the start of the feature.

The `next` command does exactly that, it moves to the next feature. If there are no more features
after the current position it doesn't rewind to the beginning (use `1` for that) and it doesn't move
to another chromosome, use `goto chrom` for that.

If `trackId` is omitted, the first annotation track is used. If trackId is not a feature track (bed,
gtf, etc) a more or less ugly warning is issued.

Find
----

#### **find_first** regex *trackId*

Find the first record in trackId containing regex. The search starts from the **end** of the current
window (so the current window is not searched) and moves forward on the current chromosome. At the end 
of the current chromosome move to the next chromosomes and then restart at from the start of
the initial one. The search stops at the first match found.

#### **find_all** regex *trackId*

Find the region containing *all* the records on chromosome containing regex. The search stops at the
first chromosome returning hits starting with the current one. Useful to get all gtf records of a
gene.

E.g. `find_all "ACTB" hg19_genes.gtf.gz` will find the entire ACTB gene:

```
hg19_genes.gtf.gz#1
eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee    ccccccccccc   cccccccccccccccccccccccc                     ccccccccccccc      ccccccc                                           eeeee
                              zz           eeeeeeeeeee   eeeeeeeeeeeeeeeeeeeeeeee                     eeeeeeeeeeeee      eeeeeee                                                
                               cccccccc                                                                                        a                                                
5566779   5566976   5567174   5567371   5567568   5567766   5567963   5568160   5568358   5568555   5568752   5568949   5569147   5569344   5569541   5569739   5569936   557013
chr7:5566779-5570232; 3,454 bp; 19.6 bp/char; Filters: -q 0 -f 0 -F 4; Mem: 1101 MB; 
[h] for help: find_all "ACTB" hg19_genes.gtf.gz
```

#### **seqRegex** *regex*

Find regex in reference sequence and show matches as and additional track. Useful to show
restriction enzyme sites, transcription factor binding motifs etc. The tag of this track is
`seqRegex` and it is not displayed. To adjust its height use `trackHeight 10 seqRegex`.

If regex is omitted the matching is disabled. This command is ignored if the reference fasta file is
missing.

Display
-------

#### **visible** *show-regex* *hide-regex* *track-regex*

In feature tracks, only include rows containing `show regex` and exclude rows containing `hide regex`.
Apply to tracks whose name is matched by `track regex`. 

With no arguments reset to default: `visible .* ^$ .*` which means show everything, hide nothing,
apply to all tracks. Use '.*' to match everything and '^$' to hide nothing. 

This command is useful to filter the annotation in GTF files, for example: 

`visible RNA mRNA gtf`

Will show the rows containing "RNA" and will hide those with mRNA in all the tracks whose name
matches gtf.

#### **trackHeight** INT *track regex*

Set track height to int lines of text for all tracks matching regex. Default regex: `.*`. 
*E.g.* `trackHeight 5 bam`.

#### **colorTrack** color *track regex*

Set colour for tracks containing regex. Default regex is `.*` (all tracks captured). Available
colours:

red, green, yellow, blue, magenta, cyan, grey,  light_red, light_green, light_yellow, light_blue,
light_magenta, light_cyan, light_grey, white, black, default

Default reset to system default colour.

E.g. `colorTrack light_blue ts.*gtf`

Colouring is rendered with ANSI codes 8/16, see also [tip colors and
formatting](http://misc.flogisoft.com/bash/tip_colors_and_formatting)

#### **ylim** min max *track regex*

Set the y-axis limit for all tracks containing regex. Use `na` to autoscale to min and/or max.
*E.g.* `ylim 0 na` will set the lower limit to 0 and the upper limit to the max of the current range.
Default `ylim na na .*`. This command applies only for tracks displaying quantitative data on y-axis
(*e.g.* bigwig, tdf).

#### **dataCol** idx *regex*

Select data column for all bedgraph tracks containing regex. idx: 1-based column index. This command
is not very useful and it might be deprecated.

#### print *track regex*

Print the rows of the feature tracks in the current window matched by *track regex*. Long lines are
clipped if they extend beyond the screen (see also [printFull](#printfull-track-regex)).

Useful to show exactly what features are present in the current window. Features are filtered in/out
according to the `visible` command. For example `print gtf`:

```
1------*------------16M-----------------32M-----------------48M-----------------64M-----------------80M-----------------95M-----------------110M----------------130M----------------140M----------------1
hg19_genes.gtf.gz#1
eeeeeeeeeeeee                         ccccccccccccccccccccccccccccccccccccccccccc                     ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccc                         eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee                     eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee
chr7 unknown exon 5566779 5567522 . - . gene_id "ACTB"; transcript_id "NM_001101"; gene_name "ACTB"; p_id "P18102"; tss_id "TSS29223";
chr7 unknown CDS  5567382 5567522 . - 0 gene_id "ACTB"; transcript_id "NM_001101"; gene_name "ACTB"; p_id "P18102"; tss_id "TSS29223";
chr7 unknown CDS  5567635 5567816 . - 2 gene_id "ACTB"; transcript_id "NM_001101"; gene_name "ACTB"; p_id "P18102"; tss_id "TSS29223";
chr7 unknown exon 5567635 5567816 . - . gene_id "ACTB"; transcript_id "NM_001101"; gene_name "ACTB"; p_id "P18102"; tss_id "TSS29223";
chr7 unknown CDS  5567912 5568350 . - 0 gene_id "ACTB"; transcript_id "NM_001101"; gene_name "ACTB"; p_id "P18102"; tss_id "TSS29223";
chr7 unknown exon 5567912 5568350 . - . gene_id "ACTB"; transcript_id "NM_001101"; gene_name "ACTB"; p_id "P18102"; tss_id "TSS29223";
```

#### **printFull** *track regex*

Same as [print](#print-track-regex) but long lines are wrapped instead of clipped.

#### **showGenome**

Print the genome dictionary, if present, with a representation of chromosome sizes. E.g.

```
[h] for help: showGenome 
chrM  16571
chr1  249250621 ||||||||||||||||||||||||||||||
chr2  243199373 |||||||||||||||||||||||||||||
...
chr21 48129895  ||||||
chr22 51304566  ||||||
chrX  155270560 |||||||||||||||||||
chrY  59373566  |||||||
```

#### **addTracks** *files or urls*

Add tracks from local or remote files.

#### **orderTracks** *track#1 track#2...*

Reorder tracks, e.g. `orderTracks #3 #2 #1` not all the tracks need to be listed. The missing ones 
follow the listed ones in unchanged order.

#### **history**

Show the list of visited positions.

Alignments
----------

These commands apply only to bam files.

#### **rpm** *track regex*

Toggle display of read coverage from raw count to reads per million for alignment files matched by 
*track regex*

#### **-f** INT and  **-F** INT 

Include (-f) and exclude (-F) reads with INT bits set, same as in `samtools`. Note that the flag 4096
can be used to filter in or out reads on top strand, this is useful in bisulfite mode.

#### **mapq** INT

Include reads with mapq >= INT, sanme as `samtools view -q`

#### **maxLines** INT

Maximum number of lines to print for alignment tracks *TODO: Shouldn't this cmd be merged with trackHeight?*

#### **BSseq** *track regex*

Toggle bisulfite mode for read tracks matched by regex. Ignored without reference fasta sequence.

In bisulfite mode, the characters M and m mark methylated bases (*i.e.* unconverted C to T) and U
and u are used for unmethylated bases (*i.e.* C converted to T). Upper case is used for reads on 
forward strand, small case for reverse. For example:

<img src="screenshots/exBSmode-2.png" width="450">


Genome option
=============

An optional genome file can be passed to option `-g/--genome` to give a set of allowed sequences and their sizes so that browsing is constrained to the real genomic space. 
The genome file is also used to represent the position of the current window on the chromosome, which is handy to navigate around.

There are three ways to pass a genome file:

* A tag identifying a built-in genome, e.g. hg19. See [genomes](http://github.com/dariober/Java-cafe/ASCIIGenome/resources/genomes) for available genomes

* A local file, tab separated with columns chromosome name and length. See [genomes](http://github.com/dariober/Java-cafe/ASCIIGenome/resources/genomes) for examples.

* A bam file with suitable header.

Note that if the input list of files contains a bam file, the `--genome` option is effectively ignored as the genome dictionary is extracted from the bam header.


Formatting of reads and features
================================

When aligned reads are show at single base resolution, read bases follow the same convention as samtools: 
Upper case letters and `.` for read align to forward strand, lower case and `,` otherwise; second-in-pair reads are underlined;
grey-shaded reads have mapping quality of <=5. 

Saving screenshots
==================

Screenshots can be saved to file with the commands `save`. Output format is either ASCII text or
png, depending on file name extension. For example:

```
[h] for help: save mygene.txt ## Save to mygene.txt as text
[h] for help: save            ## Save to chrom_start-end.txt as text
[h] for help: save .png       ## Save to chrom_start-end.png as png
[h] for help: save mygene.png ## Save to mygene.png as png
```

Without arguments, `save` writes to file named after the current  genomic position e.g.
`chr1_1000-2000.txt`.  The ANSI formatting (*i.e.* colours) is stripped before saving so that files
can be viewed on any text editor (use a monospace font like `courier`).

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

