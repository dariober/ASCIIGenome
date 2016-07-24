New in 0.2.0
=============

* Fixed important performance issue: Speed up of read bases shown on
reference.  Displaying reads on reference sequence was dead slow because  the
reference sequence was retrieved from file for *each* base of *each* read. Now
the sequence is  extracted once only. (Still not as fast `samtools tview`).

* Commands can be chained! Instead of executing them one by one you can do: `mapq 10 && -F 16 && colorTrack red` etc. 
This is more convenient and much faster because the current window is refreshed only at the end of the chain.

* Feature name shown on features. Use option `gffNameAttr` to choose which GFF attribute to use to
get name from. Now the display of feature tracks is more informative.

* Option `squash` to collapse feature with the same coordinates and `merge` to collapse features with overlapping coordinates.

* Feature track titles display the `showRegex` and `hideRegex` (so you know what you are filtering). 

* `-f -F mapq` can be applied individually to one or more bam.

* Title line of read tracks shows the number of alignments in the current window. 
Read coverage tracks show the library size and the filters applied.

* Friendlier behaviour of print/printFull: Now toggle between ON and OFF.

* For options that take a track regex, a list of regexes can be used instead. So multiple tracks need not to be captured by a single regex. 
E.g. you can do `print genes.gtf trascr.gff .*bed`

* Show help for individual commands using `myCommand -h`. E.g. `visible -h`.

