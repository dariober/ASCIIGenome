Examples
========

Open and browse 
---------------

Open some peak and bigWig files from
`ENCODE <http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeSydhTfbs/>`_. Note that
opening remote bigwig files is a little slow (IGV seems equally slow). You might also need to 
start Java with the option ` -Djava.net.useSystemProxies=true` (see also `issue#6 <https://github.com/dariober/ASCIIGenome/issues/6>`_)::

    encode=http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeSydhTfbs

    ASCIIGenome -g hg19 \
        $encode/wgEncodeSydhTfbsGm10847NfkbTnfaIggrabPk.narrowPeak.gz \
        $encode/wgEncodeSydhTfbsGm10847NfkbTnfaIggrabSig.bigWig \
        $encode/wgEncodeSydhTfbsGm12892Pol2IggmusPk.narrowPeak.gz \
        $encode/wgEncodeSydhTfbsGm12892Pol2IggmusSig.bigWig


Find the first feature on the first file, then change colour of one of the tracks. Reset y axes to
span 0 to 50, finally save as png to default file name::

    [h] for help: next #1
    [h] for help: colorTrack red wgEncodeSydhTfbsGm12892Pol2IggmusSig
    [h] for help: ylim 0 50
    [h] for help: save %r.png

Result on terminal screen should look like this:

.. image:: screenshots/chr1_996137-1003137.png

The file is to *chr1_996137-1003137.png*, note that the variable :code:`%r` is expanded to the genomic coordinates.

Finding & filtering stuff
-------------------------

Once started, :code:`ASCIIGenome` makes it easy to browse the genome. The picture below shows the distribution of transcripts on chromosome 36 of *Leishmania major*. It is clearly visible how transcripts in *Leishmania* tend to be grouped in blocks transcribed from the same direction (blue: forward strand, pink: reverse strand). Note how overlapping features are stacked on top of each other.

This screenshot has been produced by first loading the *L. major* GTF file::

    ASCIIGenome ftp://ftp.ensemblgenomes.org/pub/release-31/protists/gtf/leishmania_major/Leishmania_major.ASM272v2.31.gtf.gz

At command prompt issue the following commands::

    [h] for help: goto 36:1-2682151
    [h] for help: grep -i \ttranscript\t
    [h] for help: trackHeight 100

.. image:: screenshots/leishmania_transcripts.png

Now return to the start of the chromosome and find the first feature containing *LmjF.36.TRNAGLN.01*,
print it to screen::

    [h] for help: 1
    [h] for help: find LmjF.36.TRNAGLN.01
    [h] for help: print 

Now showing:

.. image:: screenshots/leishmania_find.png

.. _Batch-processing:

Batch processing
----------------

Often you have a list of regions to visualize in batch for a one or more tracks. For example, you
have a list of ChIP-Seq peaks or RNA-Seq genes and you want to see the coverage profiles together
with an annotation file. :code:`ASCIIGenome` allows easy batch processing  via the
:code:`--batchFile` option.

This script iterates through the intervals in *peaks.bed*. For each interval, it displays two
bigWig, a gtf file and the peak file itself.  Each interval is zoomed out 3 times and the screenshot
saved as png to :code:`/tmp/peak.%r.png`, where `%r` is a special variable  expanded to the current
coordinates as `chrom_start-end`.::

    ASCIIGenome -b peaks.bed \
        -x 'zo 3 && save /tmp/peak.%r.png' \
        chipseq.bigwig \
        input.bigwig \
        gencode_genes.gtf \
        peaks.bed > /dev/null


`convert <http://www.imagemagick.org/script/convert.php>`_ tools from ImageMagick is handy to concatenate png files and create 
a gallery of screenshots in a single file::

    convert -append /tmp/peak.*.png myPeaks.png

A similar task may be achieved by wrapping ASCIIGenome in a for-loop but it would much slower and complicated since each iteration would
require restarting the JVM and re-loading the tracks.

Finding sequence motifs
-----------------------

The reference fasta sequence can be searched for sequence motifs specified via regular expressions 
or via `IUPAC notation <https://en.wikipedia.org/wiki/Nucleic_acid_notation#IUPAC_notation>`_. 

This example is from `Biostars <https://www.biostars.org/p/221325/>`_. We want to find matches of
the motif TATAWAA near gene ENSG00000168487.

First load the reference sequence and a (remote) annotation file. Note that the fasta file must
be indexed)::

    ASCIIGenome -fa Homo_sapiens.GRCh38.dna.chromosome.8.fa \
        ftp://ftp.ensembl.org/pub/release-86/gff3/homo_sapiens/Homo_sapiens.GRCh38.86.chromosome.8.gff3.gz

Then at the command prompt issue these commands::

    find ENSG00000168487
    grep -i \tgene\t.*ENSG00000168487 gff3
    seqRegex -iupac TATAWAA
    zo 8
    print seqRegex
    seqRegex > matches.bed
    save matches.png

Explained: Find the gene ENSG00000168487, for clarity only show the "gene" feature (:code:`grep...`). 
Then search the motif TATAWAA interpreted as iupac notation; zoom out *x* times (e.g. 8 times) to see some
matches in the sequence.

The matches here are shown on screen with :code:`print seqRegex` and then saved to file with :code:`seqRegex > matches.bed`. Finally save a picture as png, shown here:

.. image:: screenshots/matches.png