TODO
====

Features
--------

* Add color to png output. It might not be too difficult: 1) Join the list of strings (lines) in a single 
string; 2) split again the big string at the end-of-formatting escape code (now each string in the least 
starts with the color of that string); 3) loop through the list and write each string using the color 
code at the beginning of the string.

The methods to look for are : `Utils.ansiColourToGraphicsColor()` to convert ansi color number to
graphics color; `Utils.convertTextFileToGraphic()` currently it strips the ansi escapes, it shouldn't if you 
want color; `Utils.convertTextToGraphic()` which loops through the list of strings.

* Accept file of regions to create images from? See https://github.com/dariober/ASCIIGenome/issues/18 and 
https://github.com/dariober/ASCIIGenome/issues/16

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
    