TODO
====

Features
--------

* Apply `print` to bam files.

* Add `-awk '<script>'` to `print` command to parse lines before printing. Useful
  to get rid of unnecessary stuff and make lines more readable.

* Indel processing. Check variants like `A,CT` processed ok.

* Converter for ensembl to ucsc chromosome names. This functionality could be included
  in `addTracks` by adding options `-ens2ucsc` and `-ucsc2ens`. E.g. `addTracks -ens2ucsc genes.gtf http://foo/bar/data.bigwig`.
  The option `-ens2ucsc|-ucsc2ens` should iterate through the target file and make an edited copy 
  where each line is converted as appropriate. Then load the edited copy by re-sending it to `addTracks`. 
  For the conversion use the conversion files also used by IGV, put them in resources.
  There are quite a few limitations: Large files, especially bams, will take ages; TDF and bigwig
  cannot easily converted, especially TDF.
  Option 2: Add a field in the `Track` class to indicate that the input file needs to
  be converted. When you fetch a region send to the track the command with the chromosome
  name changed. 

* Check implementation of picard/DownsampleSam, it might be much faster than the one you have.

* Deprecate `-fa` option? If the list of files includes a fasta file use that as reference

* Time out for `updateCheck()`. Something like:

```
long now= System.currentTimeInMiilis();
while(System.currentTimeInMiilis() < (now + 10000)){
  updateCheck();
} 
```

* Check `gffNameAttr` for gtf features

* Allow the user to set the `tmp` directory. By default use a `tmp` other then 
  the system default. On Linux the `tmp` is often pretty small.

* `colorTrack`, `trackHeight` and all commands that are applied to a list of tracks
  should accept a `-v` option to invert track selection, *i.e.* apply *command* to 
  all tracks not matched by the list of track regexes. E.g: `trackHeight -v 10 bam`
  translates to *set height to 10 for all tracks NOT matching bam*.

* Add `-I` flag to `grep` to make it case insensitive.

* Option to log transform quantitative data.

* Implement comments using e.g. // or # 

* Sampling feature: Command to display only a random sample of features in current window. 
  Useful to see the density of features in large intervals with lots of features. E.g.

```
sample INT|percent [track_regex = .*] ...
```

This cmd should operate on the feature set after grep selection. Probably `TrackIntervalFeature.getFeaturesInInterval()` is where
to work for this.

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

* Test runner `test/run_tests.sh` is effectively not used and clunky. Replace with a *e.g.* snakemake script to run tests properly.

* Populate the enum class `Command` and use `Command.cmdName` in place of the hardcoded command names.

* Each command in CommandHelp should be a class on its own extending CommandHelp. The super class command help
should implement a `validate()` method which is customized to each sub class. `.validate()` should check e.g.
that the command exists, regex are fine, number of args is ok, etc.
