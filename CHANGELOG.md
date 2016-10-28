New in 0.4.0
=============

* **Batch processing** with `--batchFile` option: ASCIIGenome will iterate through each interval in a batch file in bed or gff format. This is super useful to generate a gallery of screenshots in target regions (e.g. genes, peaks, etc). Much much faster then 
iterating ASCIIGenome in a for-loop since the JVM and files are loaded only once!
This should address almost completely [issue 18](https://github.com/dariober/ASCIIGenome/issues/18) (remember to close on merge)

* **Colours** Now the png output has colours.

* **Colours II** ASCIIGenome sets the background colour of the terminal to white, unless started with `--noFormat`. In this way the visual look of ASCIIGenome should be independent of the user's colour scheme of the terminal.

* **Memory** Memory footprint is now even smaller since files are never fully read in memory now. Bed or gff files without tabix index are sorted, block compressed and indexed as needed to temporary files. 

* Fixed [issue 20](https://github.com/dariober/ASCIIGenome/issues/20) (remember to close it on merge with trunk). 

* Fixed [issue 23](https://github.com/dariober/ASCIIGenome/issues/23) (remember to close it on merge with trunk).

* `ASCIIGenome` can start without any input file.

* TDF files can switch to be normalized to reads per million. 

* Major refactoring should make further development easier.

New in 0.2.0
=============

New commands
------------

* **Commands can be chained!** Instead of executing commands one by one at the prompt you can do:
`mapq 10 && -F 16 && colorTrack red` etc.  This is more convenient and much
faster because the current window is refreshed only at the end of the chain.

* Chained commands can also be passed at the start via the `--exec/-x` arg. Useful to set right at the start 
the configuration you want e.g. `ASCIIGenome -x 'trackHeight 5 bam && colorTrack red bigwig' aln.bam genes.gtf ...`

* **Merge & squash** features for more compact view of annotation tracks. Option `squash` collapses feature with the same coordinates and `merge` collapses overlapping features.

* **gap** command to set whether bookended features should be next to each other (more compact) or on different lines (more informative).

* **infoTracks** shows some characteristics of the loaded files.

Features
--------

* **Refit window size** automatically to fit the terminal width. Useful if you change the font size of your terminal while using ASCIIGenome.
For example make the font smaller to fit larger regions and more tracks. The view is resized at the next command executed. 

* **Features have a name!** BED and GFF feature names are shown. Use option `gffNameAttr` to choose which GFF attribute to use to
get name from. Now the display of feature tracks is more informative.

* **Help!** help for individual commands can be called with the popular syntax `myCommand -h`. E.g. `ylim -h`. Help for individual
commands is usually more verbose then the one from `h`. Also, the help invoked with `h` is better formatted and includes the default settings
of each argument

* **Title lines are more informative.** Titles of read tracks shows the number of alignments in the current window. 
Read coverage tracks show the library size and the filters applied. Useful for RNA-Seq, ChIP-Seq stuff.
Titles of Feature track titles display the `include` and `exclude` regexes (so you know what you are filtering).

* **List of regexes** For options that take a track regex, a list of regexes can be used instead.
So multiple tracks need not to be captured by a single regex.  E.g. you can do
`print genes.gtf trascr.gff .*bed`

* **Call consensus sequence.** Alignment tracks show the sequence inferred from the alignments.   

Others
------

* **Performance** Fixed important issue: Reads aligned to reference are much faster. 
Displaying reads on reference sequence was dead slow because  the
reference sequence was retrieved from file for *each* base of *each* read. Now
the sequence is  extracted once only. (Still not as fast `samtools tview`).

* Fixed bug where paired-end flags applied to non-paired reads caused the program to crash (Fixed by disabling `SAMRecord.requireReadPaired()`)

* `-f -F mapq` can be individually applied to one or more bam tracks.

* Friendlier behaviour of print/printFull: Now toggle between ON and OFF.

* Read tracks have ID ending in `@INT` (e.g. `aln.bam@3`) so it's easier to
separate them for the corresponding coverage tracks using regexes. E.g. Set
`trackHeight 10 bam# && trackHeight 1 bam@` will set height 10 for coverage
tracks and height 1 for read tracks.

* Command `visible` renamed to `filter` since this is what it actually does.

* `CG_profile` track is hidden by default (*i.e.* it has `trackHeight 0`).

* Invalid input is generally handled a bit more gently then just by making ASCIIGenome crash or throwing terrifying Java stack traces.
