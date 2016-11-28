TODO
====

Features
--------

* `next` Should move to "next" chromosome. should say something when no feature track is available.

* `print` command with `-head` option to limit the number of rows returned.

* All options that toggle a setting (print, rpm, BSseq, hideTitle, etc) should
  take two or more explicit arguments to override the toggling, e.g. `print -off *.bed`
  turn off printing for all *.bed track, `print -clip *.bed` set clip printing, etc.

* `export` command to write to file the data in the current interval? API:

```
export [track_regex = .*] > [outdir/]
```

Each track selected by regex will be written to outdir as `chrom_start-end.track_tag.ext` with appropriate `ext`ension.
Coverage and wiggle files will be bedgraph with the score column taken from the Track. Annotation files will bed or gtf
with data taken from the raw. Reference sequence if available will be in fasta format.

* Command to "go to other end" of feature? 


* Add a `nucPrint` command to print nucleotide counts at each position, a bit like pysamstats. API:
`nucPrint INT [INT] [stranded] [pct]`

Refactor
--------

* Each command in CommandHelp should be a class on its own extending CommandHelp. The super class command help
should implement a `validate()` method which is customized to each sub class. `.validate()` should check e.g.
that the command exists, regex are fine, number of args is ok, etc.
