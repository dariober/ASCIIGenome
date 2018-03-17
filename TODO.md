TODO
====

* Move files from svn_gitignore.txt to .gitignore

* Move commands to prepare test files from test_data/README.md to build.gradle 

* Read from `stdin`, useful for quick look ups like `intersectBed -a x.vcf -b y.bed | ASCIIGenome -`.
In practice, this would simply read from stdin and write to a temp file which is then sorted, compressed and indexed. 
Ok for small files, fewer then 1/2M records (?), inpractical otherwise. 
You need to guess the format of the input. Unless you enforce syntax like `ASCIIGenome -.vcf`, which is awkward.
ON HOLD - see post on [jline-users](https://groups.google.com/forum/#!topic/jline-users/jCNz0M2ioDc).
 
* Display soft clipped bases.  

FIXME
-----

* `TrackMethylation`: method `getRecordsAsStrings()` hard-coded returns empty list.
  Implement it properly or extend `TrackMethylation` to `TrackReads`

Features
--------

* Use `samtools depth` to compute coverage in `TrackPileup` over large regions. Use of `samtools depth` is avoided when filters except those in `samtools` command 
are enabled.

* Shade base quality to optionally shade a range of qualities with different shades of grey. 
  `shadeBaseQuality -below_baseq x | -range` where `-below_baseq` shades below bseq *x* as it is now.
  `-range` sets a level of grey according to quality.

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

* Deprecate `-fa` option? If the list of files includes a fasta file use that as reference

* Time out for `updateCheck()`. Something like:

```
long now= System.currentTimeInMiilis();
while(System.currentTimeInMiilis() < (now + 10000)){
  updateCheck();
} 
```

* Allow the user to set the `tmp` directory.

* Option to log transform quantitative data.

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
