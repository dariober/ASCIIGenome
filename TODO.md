TODO
====

Features
--------

* Alogn with 

* Allow a string of options to be set at prompt. E.g. `-F 16; mapq 10; ylim 0 na`

* Add a `bookmark` command create on the fly an IntervalFeatureTrack with bookmarked regions.  API:
`bookmark`: Add current region to bookmarks; `bookmark show` print current bookmarked features;
`bookmark clear`, etc.

* Command to "go to other end" of feature? 

* Enable (some) options to be set at start

* Add `hideTitle/showTitle trackRegex#1 trackRegex#2 ...` to hide the title line from selected tracks to allow more compact visualization? 

* ~~Remove "exon" feature from ideogram if there is also a CDS feature with the same coords? This makes the view less crowded but it could hide information?~~ See squash

Refactor
--------

* Tidy up handling of exceptions when parsing interactive command line options. 

* Create an enumerate class for all the commands. Use that instead of hardcoding strings.