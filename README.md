*Note: Versions 0.1.0 and 0.2.0a are mostly experimental. Better version(s) coming soon!*

Text Only Genome Viewer!
========================


- [Description](#description)
- [Requirements and Installation](#requirements-and-installation)
- [Installation quick start](#installation-quick-start)
- [Installation through Homebrew](#installation-through-homebrew)
- [A little more detail](#a-little-more-detail)
- [Usage examples](#usage-examples)
    - [Minimal example](#minimal-example)
    - [Open and browse](#open-and-browse)
    - [Finding & filtering stuff](#finding--filtering-stuff)
- [Supported input](#supported-input)
- [Genome option](#genome-option)
- [Formatting of reads and features](#formatting-of-reads-and-features)
- [Saving screenshots](#saving-screenshots)
- [Tips gotchas and miscellanea](#tips-gotchas-and-miscellanea)
- [Interactive commands](#interactive-commands)
    - [Navigation](#navigation)
        - [f and b](#f-and-b)
        - [ff and bb](#ff-and-bb)
        - [zi *x* and zo *x*](#zi-x-and-zo-x)
        - [goto chrom:*from-to*](#goto-chromfrom-to)
        - [INT *INT*](#int-int)
        - [+/- INT *k,m*](#--int-km)
        - [p and n](#p-and-n)
        - [next *trackId* and next_start *trackId*](#next-trackid-and-next_start-trackid)
    - [Find](#find)
        - [find_first regex *trackId*](#find_first-regex-trackid)
        - [find_all regex *trackId*](#find_all-regex-trackid)
        - [seqRegex *regex*](#seqregex-regex)
    - [Display](#display)
        - [visible *show_regex* *hide_regex* *track_regex*](#visible-show_regex-hide_regex-track_regex)
        - [squash *track_regex*](#squash-track_regex)
        - [gffNameAttr attribute_name *track_regex*](#gffnameattr-attribute_name-track_regex)
        - [trackHeight INT *track_regex*](#trackheight-int-track_regex)
        - [ylim min max *track_regex*](#ylim-min-max-track_regex)
        - [colorTrack color *track_regex*](#colortrack-color-track_regex)
        - [dataCol idx *regex*](#datacol-idx-regex)
        - [print *track_regex*](#print-track_regex)
        - [printFull *track_regex*](#printfull-track_regex)
        - [showGenome](#showgenome)
        - [addTracks *files or urls*](#addtracks-files-or-urls)
        - [orderTracks *regex#1 regex#2 ...*](#ordertracks-regex1-regex2-)
        - [history](#history)
    - [Alignments](#alignments)
        - [rpm *track_regex*](#rpm-track_regex)
        - [-f INT and  -F INT](#-f-int-and---f-int)
        - [mapq INT](#mapq-int)
        - [BSseq *track_regex*](#bsseq-track_regex)
- [Credits](#credits)

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

`ASCIIGenome` is a command-line genome browser running from terminal window and solely based on
ASCII characters. Since `ASCIIGenome` does not require a graphical interface it is particularly
useful for  quickly visualizing genomic data on remote servers. The idea is to make `ASCIIGenome`
the Vim  of genome viewers.

As far as I know, the closest program to `ASCIIGenome` is [samtools tview](http://samtools.sourceforge.net/tview.shtml) but 
`ASCIIGenome` offers much more flexibility, similar to popular GUI viewers like [IGV](https://www.broadinstitute.org/igv/).

Some key features:

* Command line input and interaction, no graphical interface, minimal [installation and requirements](#requirements-and-installation)
* Can load multiple files in various [formats](#supported-input)
* Can access remote files via URL or ftp address
* Easy [navigation](#navigation) and [searching](#find) of features and sequence motifs and filtering options
* Support for BS-Seq alignment

<img src="screenshots/ex3.png" width="800">

Requirements and Installation
=============================

Installation quick start 
------------------------

In the commands below replace version number with the latest from [releases](https://github.com/dariober/ASCIIGenome/releases):

```
wget https://github.com/dariober/ASCIIGenome/releases/download/v0.1.0/ASCIIGenome-0.1.0.zip
unzip ASCIIGenome-0.1.0.zip

cd ASCIIGenome-0.1.0/
chmod a+x ASCIIGenome
cp ASCIIGenome.jar /usr/local/bin/ # Or ~/bin/
cp ASCIIGenome /usr/local/bin/     # Or ~/bin/ 
```

Installation through Homebrew
------------------------------

ASCIIGenome can also be installed through [brew](http://brew.sh/) / [Linux Brew](https://github.com/Linuxbrew/brew), although it is still not an official package:

```
brew install https://raw.githubusercontent.com/dariober/ASCIIGenome/master/install/brew/asciigenome.rb
```


A little more detail
--------------------

`ASCIIGenome.jar` requires **Java 1.7+** and this should be the only requirement. There is virtually no installation needed as `ASCIIGenome` is pure Java and should work on most (all?) platforms. Download the zip file `ASCIIGenome-x.x.x.zip` from [releases](https://github.com/dariober/ASCIIGenome/releases), unzip it and execute the jar file with

```
java -jar /path/to/ASCIIGenome.jar --help
```

To avoid typing `java -jar ...` every time, you can put both the helper 
script `ASCIIGenome` and the jar file ```ASCIIGenome.jar``` in the same directory in your `PATH` and execute with:

```
ASCIIGenome [options]
```

Note the helper is a bash script. To set the amount of memory available to java use the `-Xmx` option as e.g. `java -Xmx1500m -jar ...`.

If for some reason the text formatting misbehaves, disable it with the `-nf` option. I have developed 
ASCIIGenome on MacOS, Ubuntu and CentOS with `bash 4.1`, white colour background.

Usage examples
==============

These are just some functionalities to give an idea behind ASCIIGenome.

### Minimal example

Open a bam file, as simple as:

```
ASCIIGenome aln.bam
```

Open with a reference genome (reference must be indexed, see [Supported input](#supported-input)):

```
ASCIIGenome -fa genome.fa aln.bam
```


### Open and browse 

Open some peak and bigWig files from
[ENCODE](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeSydhTfbs/). Note that
opening remote bigwig files is a little slow (IGV seems equally slow) and it might not work 
with some proxy settings (see also [issue#6](https://github.com/dariober/ASCIIGenome/issues/6)):

```
encode=http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeSydhTfbs

ASCIIGenome -g hg19 \
    $encode/wgEncodeSydhTfbsGm10847NfkbTnfaIggrabPk.narrowPeak.gz \
    $encode/wgEncodeSydhTfbsGm10847NfkbTnfaIggrabSig.bigWig \
    $encode/wgEncodeSydhTfbsGm12892Pol2IggmusPk.narrowPeak.gz \
    $encode/wgEncodeSydhTfbsGm12892Pol2IggmusSig.bigWig
```

Find the first feature on the first file, then change colour of one of the tracks. Reset y axes to
span 0 to 50, finally save as png to default file name:

```
[h] for help: next #1
[h] for help: colorTrack magenta wgEncodeSydhTfbsGm12892Pol2IggmusSig
[h] for help: ylim 0 50
[h] for help: save .png
```

Result on terminal screen should look like this:

<img src="screenshots/encode.png" width="800">

Saved file is `chr1_996137-1003137.png` (currently the png output doesn't include colours though).

### Finding & filtering stuff

Once started, ```ASCIIGenome``` makes it easy to browse the genome. The picture below shows the distribution of transcripts on chromosome 36 of *Leishmania major*. It is clearly visible how transcripts in *Leishmania* tend to be grouped in blocks transcribed from the same direction (blue: forward strand, pink: reverse strand). Note how overlapping features are stacked on top of each other.

This screenshot has been produced by first loading the *L. major* GTF file:

```
ASCIIGenome ftp://ftp.ensemblgenomes.org/pub/release-31/protists/gtf/leishmania_major/Leishmania_major.ASM272v2.31.gtf.gz
```

At command prompt issue the following commands:

```
[h] for help: goto 36:1-2682151
[h] for help: visible '\ttranscript\t'
[h] for help: trackHeight 100
```

<img src="screenshots/leishmania_transcripts.png" width="800">

Now return to the start of the chromosome and find the first feature containing *LmjF.36.TRNAGLN.01*,
print it to screen:

```
[h] for help: 1
[h] for help: find_first LmjF.36.TRNAGLN.01
[h] for help: print 
```

Now showing:

<img src="screenshots/leishmania_find.png" width="800">

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

Bed & gtf file are not required to be sorted or index but in this case they are loaded in memory. To
save memory and time for large files you can again index them as above. Loading in memory is
typically fast for files of up to ~1/2 million records.

For input format specs see also [UCSC format](https://genome.ucsc.edu/FAQ/FAQformat.html) and for
guidelines on the choice of format see [IGV
recommendations](https://www.broadinstitute.org/igv/RecommendedFileFormats).

**Fasta reference**: The reference sequence should be uncompressed and indexed, with *e.g.* [samtools faidx](http://www.htslib.org/doc/samtools.html):

```
samtools faidx genome.fa
```

Genome option
=============

An optional genome file can be passed to option `-g/--genome` to give a set of allowed sequences and their sizes so that browsing is constrained to the real genomic space. 
The genome file is also used to represent the position of the current window on the chromosome, which is handy to navigate around.

There are three ways to pass a genome file:

* A tag identifying a built-in genome, e.g. hg19. See [genomes](https://github.com/dariober/ASCIIGenome/tree/master/resources/genomes) for available genomes

* A local file, tab separated with columns chromosome name and length. See [genomes](https://github.com/dariober/ASCIIGenome/tree/master/resources/genomes) for examples.

* A bam file with suitable header.

Note that if the input list of files contains a bam file, the `--genome` option is effectively ignored as the genome dictionary is extracted from the bam header.



Formatting of reads and features
================================

When aligned reads are show at single base resolution, read bases follow the same convention as samtools: 
Upper case letters and `.` for read align to forward strand, lower case and `,` otherwise; second-in-pair reads are underlined;
grey-shaded reads have mapping quality of <=5. 

GTF/GFF features on are coded according to the feature column as below. For forward strand 
features the colour blue and upper case is used, for reverse strand the colour is pink the case is lower. 
Features with no strand information are in grey.

Feature | Symbol
--------|-------
exon | E  
cds | C  
start_codon | A 
stop_codon | Z 
utr | U 
3utr | U 
5utr | W 
gene | G 
transcript | T 
mrna | M 
trna | X 
rrna | R 
mirna | I 
ncrna | L 
lncrna | L 
sirna | S 
pirna | P 
snorna | O 


If available, the feature name is shown on the feature itself. 
The feature name has a trailing underscore to separate it from the rest of the feature representation. The
last character of the feature is always the feature type. For example, the feature named `myGene` appears as:

```
myGene_EEEEEEEEE ## Enough space for the full name
myGenE           ## Not enough space, name truncated and last char is E
```

For BED features, name is taken from column 4, if available. Default for GTF/GFF is to take name from attribute 
`Name`, if absent try: `ID`, `transcript_name`, `transcript_id`, `gene_id`, `gene_name`. 
To choose an attribute see command `gffNameAttr`.

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


Tips gotchas and miscellanea
============================

* **Performance** Alignment files are typically accessed very quickly but `ASCIIGenome` becomes slow
when the window size grows above a few hundreds of kilobases. Annotation files (bed, gff, gtf) are
loaded in memory unless they are indexed with `tabix`.

* **Regular expression** Use the `(?i)` modifier to match in case insensitve mode, e.g. '(?i).*actb.*'

* When displaying bam files, `ASCIGenome` is hardcoded to disable the coverage and read tracks if
the window size is >100,000 bp. This is to prevent the browsing to become horribly slow. To display
such large windows  consider bigWig or tdf file format.

* When opening bam files, the first chromosome is often the mitochondrial chromosome chrM (or chrMT) which
often has very high read depth (say 10,000x). This can make the opening slow. Consider using the `-r`
option in these cases. E.g. `ASCIIGenome -r chr1 file1.bam file2.bam ...`

Interactive commands
====================

As there is no GUI, everything is handled thorough command line. Once `ASCIIGenome` is started enter
a command and press ENTER to execute.

Some features of Unix console are enabled: 

* Arrow keys UP and DOWN scroll previous commands.
* TAB auto-completes commands.
* ENTER without any argument repeats the previous command.

Examples:

```
[h] for help: ff <ENTER>   ## Move forward
[h] for help: <ENTER>      ## Move forward again...
[h] for help: <ENTER>      ## ... and again
[h] for help: col <TAB>    ## Is expanded to colorTrack
[h] for help: <ARROW UP>   ## Shows previous command
[h] for help: h <ENTER>    ## Show help.
```

In the documentation below, mandatory arguments are in plain style while  optional argument are in
*italics*. When track names are passed as arguments, it is not necessary to  give the full name as
partial matching is enabled. This is handy since track names have an ID appended as suffix which can
be used in place of the full name, e.g. `next myLongfileName.bed#1` can be also typed as `next #1`.

Navigation
----------

#### f and b

Move forward or backward by 1/10 of a window

#### ff and bb

Move forward or backward by 1/2 of a window

#### zi *x* and zo *x*

Zoom in / zoom out *x* times. Each zoom halves or doubles the window size. To zoom in or out really
quickly use x= 5 or 10 like `zi 10` or more. Default x is 1.

#### goto chrom:*from-to*

Go to region *chrom:from-to* or to position *chrom:from* or to start of chromosome *chrom*.
The character ':' is a shortcut for `goto`. Examples:

```
goto chr8:1-1000
goto chr8:1
goto chr8
## The same:
:chr8:1-1000
:chr8:1
:chr8
```

#### INT *INT*

Go to position `INT` or to region `INT INT` on **current chromosome**. Also allowed is the
hyphenated format  separating the two positions. If a list of integers is given, the first and last
are taken as *from* and *to*. This is handy to copy and paste intervals from the ruler above the
prompt. Examples:

```
[h] for help: 10                   ## Will jump to position 10
[h] for help: 10 1000              ## Go to region 10-1000
[h] for help: 10-1000              ## Same as above
[h] for help: 10 250 500 750 1000  ## Same as above again
```

#### +/- INT *k,m*

Move forward/backward by INT bases. Suffixes k (kilo) and M (mega) expanded to x1000 and
x1,000,000. *E.g.* `-2m` or `+10k` or `+10.5k`.

#### p and n

Go to the **p**revious or **n**ext visited position. Similar to the back and forward arrows of an
Internet browser.

#### next *trackId* and next_start *trackId*

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

#### find_first regex *trackId*

Find the first record in trackId containing regex. The search starts from the **end** of the current
window (so the current window is not searched) and moves forward on the current chromosome. At the end 
of the current chromosome move to the next chromosomes and then restart at from the start of
the initial one. The search stops at the first match found.

#### find_all regex *trackId*

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

See below for the encoding of GTF features.


#### seqRegex *regex*

Find regex in reference sequence and show matches as and additional track. Useful to show
restriction enzyme sites, transcription factor motifs, etc. The tag of this track is
`seqRegex` and it is not displayed. To adjust its height use `trackHeight 10 seqRegex`. Matching
is case sensitive, to ignore case use the regex syntax `(?i)` as in `seqReg (?i)ACTG`.

If regex is omitted the matching is disabled. This command is ignored if the reference fasta file is
missing.

For example, find all instances of the motif `GG..T`:

```
1------*------------16M-----------------32M-----------------48M-----------------64M-----------------80M-----------------95M-----------------110M----------------130M----------------140M----------------1
                     <<<<<   >>>>>  >>>>>   >>>>>                      >>>>>      <<<<<    >>>>>                                     >>>>>                                      <<<<<      >>>>>         
                                              <<<<<                                                                                                                                                      
ATACAGCTTGCAGCTGAGCAGAGGCCCTTGGATTTGGGAGTTTGGGAATCCTCTTTTCCCTTGATTCTACTGGATTGCATTCAGTCCAGATGGAATGAGGTTGTGGGAAGGGGAGCCAGGATGAACTCTTAATGGCTTGGAAACAGATCTGATCTCCTAAGAGGACAGGAAGTCAGACTCCTTTTGTGGGTTACAGGGGAA
5567802   5567812   5567822   5567832   5567842   5567852   5567862   5567872   5567882   5567892   5567902   5567912   5567922   5567932   5567942   5567952   5567962   5567972   5567982   5567992   5
chr7:5567802-5568002; 201 bp; 1.0 bp/char; Filters: -q 0 -f 0 -F 4; Mem: 3124 MB; 
[h] for help: seqRegex GG..T
```

Display
-------

#### visible *show_regex* *hide_regex* *track_regex*

In feature tracks, only include rows containing `show regex` and exclude rows containing `hide regex`.
Apply to tracks whose name is matched by `track regex`. 

With no arguments reset to default: `visible .* ^$ .*` which means show everything, hide nothing,
apply to all tracks. Use '.*' to match everything and '^$' to hide nothing. 

This command is useful to filter the annotation in GTF files, for example: 

`visible RNA mRNA gtf`

Will show the rows containing "RNA" and will hide those containing "mRNA", applies to tracks whose name
matches "gtf".

#### squash *track_regex*

For the feature tracks captured by *track_regex*, toggle the squashing of features with the same coordinates (same
start, end, and strand) to allow more compact representation. When multiple features are squashed, only the first one
is shown. Default regex is `.*`. This is an example if features squashed and unsquashed:

```
hg19_genes_head.gtf#1; Show '.*' Hide '' ; squashed
                    NR_026820_eeeeeeee
hg19_genes_head.gtf#2; Show '.*' Hide '' 
                    NR_026820_eeeeeeee
                    NR_026818_eeeeeeee
33967     34296     34625     34954     35283     35612
chr1:33967-39823; 5,857 bp; 32.7 bp/char;
```

#### gffNameAttr attribute_name *track_regex*

For GTF/GFF tracks, choose the attribute to get the feature name from. Use attribute NULL to
reset to default choice of attribute. Applies to all tracks captured by `track_regex`. 
Ignored by non GFF/GTF features. Default `gffNameAttr NULL .*`. 
For examples given the feature:

```
chr1 unknown exon 11874 12227 . + . gene_id "DDX11L1"; Name "NR_046018_1"; gene_name "DDX11L1";
```

Default:

```
gffNameAttr NULL 
E_NR_046018_1_EEEEEE
```

Use `gene_id` attribute:

```
gffNameAttr gene_id 
E_DDX11L1_EEEEEEEEEE
```



#### trackHeight INT *track_regex*

Set track height to int lines of text for all tracks matching regex. Default regex: `.*`. 
*E.g.* `trackHeight 5 bam`.

#### ylim min max *track_regex*

Set the y-axis limit for all tracks containing regex. Use `na` to autoscale to min and/or max.
This command applies only to tracks displaying quantitative data on y-axis (*e.g.* bigwig, tdf), the other
tracks are unaffected. Examples:

```
ylim 0 50      ## Set min= 0 and max= 50 in all tracks.
ylim 0 na      ## Set min to 0 and autoscale the max. Apply to all tracks
ylim na na tdf ## Autoscale min and max. Apply to all tracks matching "tdf"
```

#### colorTrack color *track_regex*

Set colour for tracks containing regex. Default regex is `.*` (all tracks captured). Available
colours:

red, green, yellow, blue, magenta, cyan, grey,  light_red, light_green, light_yellow, light_blue,
light_magenta, light_cyan, light_grey, white, black, default

The "default" colour reset to the system default colour.

E.g. `colorTrack light_blue ts.*gtf`

Colouring is rendered with ANSI codes 8/16, see also [tip colors and
formatting](http://misc.flogisoft.com/bash/tip_colors_and_formatting)

#### dataCol idx *regex*

Select data column for all bedgraph tracks containing regex. idx: 1-based column index. This command
is not very useful and it might be deprecated.

#### print *track_regex*

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
5567470   5567513   5567556   5567600   5567643   5567686   5567729   5567772   5567816   5567859   5567902   5567945   5567988   5568032   5568075   5568118   5568161   5568204   5568248   5568291   5
chr7:5567470-5568334; 865 bp; 4.3 bp/char; Filters: -q 0 -f 0 -F 4; Mem: 3047 MB; 
[h] for help: print gtf
```
To turn off the printing give a regex that does not capture the tracks to be turned off (`print ^$` will
turn off all the tracks).

#### printFull *track_regex*

Same as [print](#print-track-regex) but long lines are wrapped instead of clipped.

#### showGenome

Print the genome dictionary with a representation of chromosome sizes. E.g.

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

#### addTracks *files or urls*

Add tracks from local or remote files.

#### orderTracks *regex#1 regex#2 ...*

Reorder tracks according to the list of regexes. Not all the tracks need to be listed, the missing ones 
follow the listed ones in unchanged order. 

For example, given the track list: `[hela.bam#1, hela.bed#2, hek.bam#3, hek.bed#4]`:

```
orderTracks #2 #1   -> [hela.bed#2, hela.bam#1, hek.bam#3, hek.bed#4]
orderTracks bam bed -> [hela.bam#1, hek.bam#3, hela.bed#2, hek.bed#4]
```

#### history

Show the list of visited positions.

Alignments
----------

These commands apply only to bam files.

#### rpm *track_regex*

Toggle display of read coverage from raw count to reads per million for alignment files matched by 
*track regex*

#### -f INT and  -F INT 

Include (-f) and exclude (-F) reads with INT bits set, same as in `samtools`. Note that the flag 4096
can be used to filter in or out reads on top strand, this is useful in bisulfite mode.

#### mapq INT

Include reads with mapq >= INT, same as `samtools view -q`

#### BSseq *track_regex*

Toggle bisulfite mode for read tracks matched by regex. Ignored without reference fasta sequence.

In bisulfite mode, the characters M and m mark methylated bases (*i.e.* unconverted C to T) and U
and u are used for unmethylated bases (*i.e.* C converted to T). Upper case is used for reads on 
forward strand, small case for reverse. For example:

<img src="screenshots/exBSmode-2.png" width="450">


Credits
=======

* Bam processing is mostly done with the [samtools/htsjdk](https://github.com/samtools/htsjdk) library.
* Bigwig and tdf are processed with classes from [IGV](https://github.com/igvteam/igv) source code.
* Block compression and indexing done using [jvarkit](https://github.com/lindenb/jvarkit).
* Brew installation thanks to [dalloliogm](https://github.com/dalloliogm).
