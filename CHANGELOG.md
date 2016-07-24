New in 0.2.0
=============

* **Performance** Fixed important issue: Reads aligned to reference are much faster. 
Displaying reads on reference sequence was dead slow because  the
reference sequence was retrieved from file for *each* base of *each* read. Now
the sequence is  extracted once only. (Still not as fast `samtools tview`).

* **Commands can be chained!** Instead of executing commands one by one at the prompt you can do:
`mapq 10 && -F 16 && colorTrack red` etc.  This is more convenient and much
faster because the current window is refreshed only at the end of the chain.

* Chained commands can be passed to at start via `--exec` arg. Useful to set right at the start 
the configuration you want e.g. `ASCIIGenome -x 'trackHeight 5 bam && colorTrack red bigwig'`

* **Features have a name** BED and GFF feature names are shown. Use option `gffNameAttr` to choose which GFF attribute to use to
get name from. Now the display of feature tracks is more informative.

* **Merge & squash** features. Option `squash` collapses feature with the same
* **coordinates and `merge` collapses overlapping features.

* Feature track titles display the `showRegex` and `hideRegex` (so you know what you are filtering). 

* `-f -F mapq` can be individually applied to one or more bam tracks.

* **Title line** of read tracks shows the number of alignments in the current window. 
Read coverage tracks show the library size and the filters applied. Useful for RNA-Seq, ChIP-Seq stuff.

* Friendlier behaviour of print/printFull: Now toggle between ON and OFF.

* **List of regexes** For options that take a track regex, a list of regexes can be used instead.
So multiple tracks need not to be captured by a single regex.  E.g. you can do
`print genes.gtf trascr.gff .*bed`

* **Help!** Show help for individual commands using `myCommand -h`. E.g. `visible -h`.

* Read tracks have ID ending in `@INT` (e.g. `aln.bam@3`) so it's easier to
separate them for the corresponding coverage tracks using regexes. E.g. Set
`trackHeight 10 bam# && trackHeight 1 bam@` will set height 10 for coverage
tracks and height 1 for read tracks.

* **New commands** hideTracks, infoTracks