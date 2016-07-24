TODO
====

Features
--------

* Add a `bookmark` command create on the fly an IntervalFeatureTrack with bookmarked regions.  API:
`bookmark`: Add current region to bookmarks; `bookmark show` print current bookmarked features;
`bookmark clear`, etc.

* Command to "go to other end" of feature? 

* Enable (some) options to be set at start

* Add `hideTitle/showTitle trackRegex#1 trackRegex#2 ...` to hide the title line from selected
 tracks to allow more compact visualization?

* `next` and `next_start` should say something when no feature track is available or you are no more
features to move to.

Refactor
--------

* Each command in CommandHelp should be a class on its own extending CommandHelp. The super class command help
implements a `validate()` method which is customized to each sub class. `.validate()` should check e.g.
that the command exists, regex are fine, number of args is ok, etc.

* Tidy up handling of exceptions when parsing interactive command line options. 

