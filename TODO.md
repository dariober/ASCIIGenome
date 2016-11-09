TODO
====

Features
--------

* `export` command to write to file the data in the current interval? API:

```
export [track_regex = .*] > [outdir/]
```

Each track selected by regex will be written to outdir as `chrom_start-end.track_tag.ext` with appropriate `ext`ension.
Coverage and wiggle files will be bedgraph with the score column taken from the Track. Annotation files will bed or gtf
with data taken from the raw. Reference sequence if available will be in fasta format.

* Command to "go to other end" of feature? 

* `next` should say something when no feature track is available or there are no more features to move to.

* Add a `nucPrint` command to print nucleotide counts at each position, a bit like pysamstats. API:
`nucPrint INT [INT] [stranded] [pct]`

Refactor
--------

* ~Tracks CG_profile and seqRegex should be part of TrackSet like all the others. Currently they are treated 
separately and this is a mess.~

* Each command in CommandHelp should be a class on its own extending CommandHelp. The super class command help
should implement a `validate()` method which is customized to each sub class. `.validate()` should check e.g.
that the command exists, regex are fine, number of args is ok, etc.

* Tidy up handling of exceptions when parsing interactive command line options (it's already better now then in 0.1.0). 
    