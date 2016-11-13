.. _Supported_input_and_output:

Input and output
================

Input file formats
------------------

File name extensions matter as file types are usually recognized by their extension in case
insensitive mode. Reading remote files might require starting java with the option
:code:`-Djava.net.useSystemProxies=true`  (see `issue#6 <https://github.com/dariober/ASCIIGenome/issues/6>`_).

Unless noted otherwise remote URL links are supported. 

* **bam** files must be sorted and indexed, *e.g.* with :code:`samtools sort` and :code:`samtools index`. 
  Paths to remote URLs are supported but painfully slow (*IGV seems to suffer of the same issue*). Loaded
  bam files generate two tracks: One for read coverage profile and one for the aligned reads.

* **bigWig** recognized by extension :code:`.bw` or :code:`.bigWig`.

* **bedGraph** recognized by extension :code:`.bedGraph` or :code:`.bedgraph`

* **bed**, **gtf**, **gff** are recognized by respective extensions.

* **tdf** This is very useful for quickly displaying very large intervals like tens of megabases or entire chromosomes (see `tdf here <https://www.broadinstitute.org/igv/TDF>`_).

* **vcf** Supported but their representation is not too sophisticated representation. URL should be supported but it appears ftp from 1000genomes doesn't work (same for IGV).

* All other file extensions (e.g. txt, narrowPeak) will be assumed to be bed formatted files.

All plain text formats (bed, bedgraph, etc) can be read as gzipped and there is no need to decompress them.

For input format specs see also `UCSC format <https://genome.ucsc.edu/FAQ/FAQformat.html>`_ and for
guidelines on the choice of format see the `IGV recommendations <https://www.broadinstitute.org/igv/RecommendedFileFormats>`_.

Some notable formats currently **not** supported are cram and bigBed. bigBed files can be converted to bgzip format with :code:`bigBedToBed` from 
`UCSC utilities <http://hgdownload.soe.ucsc.edu/admin/exe/>`_ and then indexed with tabix (see :ref:`handling_large_files`).

.. _handling_large_files:

Handling large files
++++++++++++++++++++

ASCIIGenome always makes use of indexing to access data so that the memory usage stays low even for large
files. Tabular data files (bed, gtf/gff, bedgraph, etc) which are not block compressed and indexed
are first sorted, compressed and indexed to temporary working files. This is usually quick for files of
up to 1/2 million rows. For larger files consider compressing them with external utilities such as 
`tabix <http://www.htslib.org/doc/tabix.html>`_. For example, to sort, compress and index a bed
file::

    sort -k1,1 -k2,2n my.bed \
    | bgzip > my.bed.gz
    tabix -p bed my.bed.gz

Alternatively, the temporary block compressed files and their indexes created by ASCIIGenome can be saved
by copying them from their temporary location. Use the command `infoTracks <link here>`_ to see where the
tmp working files live.

Setting a genome
----------------

An optional genome file can be passed to option :code:`-g/--genome` or set with the
:code:`setGenome` command to give a set of allowed sequences and their sizes so that browsing is
constrained to the real genomic space.  The genome file is also used to represent the position of
the current window on the chromosome, which is handy to navigate around.

There are different ways to set a genome:

* A tag identifying a built-in genome, e.g. hg19. 
  See `genomes <https://github.com/dariober/ASCIIGenome/tree/master/resources/genomes>`_ for available genomes.

* A local file, tab separated with columns chromosome name and length. 
  See `genomes <https://github.com/dariober/ASCIIGenome/tree/master/resources/genomes>`_ for examples.

* A bam file with suitable header.

* A fasta reference sequence (see :ref:`Fasta-reference-sequence`).

.. _Fasta-reference-sequence:

Fasta reference sequence
++++++++++++++++++++++++

The reference sequence, optional, should be uncompressed and indexed, with *e.g.* `samtools faidx <http://www.htslib.org/doc/samtools.html>`_::

    samtools faidx genome.fa


Output: Formatting of reads and features
----------------------------------------

When aligned reads are show at single base resolution, read bases follow the same convention as
samtools:  Upper case letters and :code:`.` for read align to forward strand, lower case and
:code:`,` otherwise; second-in-pair reads are underlined; grey-shaded reads have mapping quality of <=5.

GTF/GFF features on are coded according to the feature column as below. For forward strand  features
the colour blue and upper case is used, for reverse strand the colour is pink and the case is lower.
Features with no strand information are in grey.

===========  ======
Feature      Symbol
===========  ======
exon         E  
cds          C  
start_codon  A 
stop_codon   Z 
utr          U 
3utr         U 
5utr         W 
gene         G 
transcript   T 
mrna         M 
trna         X 
rrna         R 
mirna        I 
ncrna        L 
lncrna       L   
sirna        S 
pirna        P 
snorna       O 
===========  ======

This is an example of a BAM file and a GTF file. The top track shows the read coverage, the middle
track the aligned reads and the bottom track the GTF features in this genomic window.

.. image:: screenshots/actb_bam_gtf.png

If available, the feature name is shown on the feature itself.  The feature name has a trailing
underscore to separate it from the rest of the feature representation. The last character of the
feature is always the feature type. For example, the feature named :code:`myGene` appears as::

    myGene_EEEEEEEEE ## Enough space for the full name
    myGenE           ## Not enough space, name truncated and last char is E

For BED features, name is taken from column 4, if available. Default for GTF/GFF is to take name
from attribute  :code:`Name`, if absent try: :code:`ID`, :code:`transcript_name`,
:code:`transcript_id`, :code:`gene_id`, :code:`gene_name`.  To choose an attribute see command
:code:`gffNameAttr`.

Read coverage tracks at single base resolution show the consensus sequence obtained from the
underlying reads. If the reference fasta file is present the :code:`=` symbol is used to denote a
match. Heterozygote bases or variants are shown  using the [iupac ambiguity
codes](http://www.bioinformatics.org/sms/iupac.html) for up to two variants (N otherwise). Variants
are called with a not-too-sophisticated heuristics: Only base qualities >= 20 are considered, an
alternative allele is called if supported by at least 3 reads and makes up at least 1% of the total
reads. The first and second allele must make at least  98% of the total reads otherwise the base is
N (see :code:`PileupLocus.getConsensus()` for exact implementation). Insertion/deletions are
currently not considered.

Saving screenshots
------------------

Screenshots can be saved to file with the commands :code:`save`. Output format is either ASCII text or
png, depending on file name extension. For example::

    [h] for help: save mygene.txt ## Save to mygene.txt as text
    [h] for help: save            ## Save to chrom_start-end.txt as text
    [h] for help: save .png       ## Save to chrom_start-end.png as png
    [h] for help: save mygene.png ## Save to mygene.png as png

Without arguments, :code:`save` writes to file named after the current  genomic position e.g.
`chr1_1000-2000.txt`.  The ANSI formatting (*i.e.* colours) is stripped before saving so that files
can be viewed on any text editor (use a monospace font like :code:`courier`). For convenience the 
variable :code:`%r` in the file name is expanded to the current genomic coordinates, for example 
`save mygene.%r.png` is expanded to *e.g.* :code:`mygene.chr1_1000_2000.png`. 

See also :ref:`Batch-processing` for saving screenshots in batch by iterating through a list of
positions.

This is a screenshots of bisulfite-seq data. The `BSseq` mode was set and methylated cytosines are shown in red while unmethylated cytosines in blue. 

.. image:: screenshots/bs.chr7_5560313-5560467.png


