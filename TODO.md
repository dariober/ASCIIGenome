TODO
====

Features
--------

Need to turn off DEBUG when reading fasta files!

```
DEBUG   2016-07-31 22:23:19 BlockCompressedOutputStream Using deflater: Deflater
DEBUG   2016-07-31 22:23:19 BlockCompressedOutputStream Using deflater: Deflater
```

* Add a `bookmark` command create on the fly an IntervalFeatureTrack with bookmarked regions.  API:
`bookmark`: Add current region to bookmarks; `bookmark show` print current bookmarked features;
`bookmark clear`, etc.

* Command to "go to other end" of feature? 

* Add `hideTitle/showTitle trackRegex#1 trackRegex#2 ...` to hide the title line from selected
 tracks to allow more compact visualization?

* `next` and `next_start` should say something when no feature track is available or there are no more
features to move to.

* Add a `nucPrint` command to print nucleotide counts at each position, a bit like pysamstats. API:
`nucPrint INT [INT] [stranded] [pct]`

Refactor
--------

* Tracks CG_profile and seqRegex should be part of TrackSet like all the others. Currently they are treated 
separately and this is a mess.

* Each command in CommandHelp should be a class on its own extending CommandHelp. The super class command help
should implement a `validate()` method which is customized to each sub class. `.validate()` should check e.g.
that the command exists, regex are fine, number of args is ok, etc.

* Tidy up handling of exceptions when parsing interactive command line options (it's already better now then in 0.1.0). 
    