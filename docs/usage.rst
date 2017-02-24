Usage
=====

Quick start
-----------

The interface to *ASCIIGenome* should look familiar to those used to command line programs.  To see
the command line options use the usual syntax :code:`-h/--help` as :code:`ASCIIGenome -h`.

Open an indexed a bam file, as simple as::

    ASCIIGenome aln.bam

Open with a reference genome::

    ASCIIGenome -fa genome.fa aln.bam

In interactive mode use :code:`-h` to browse available commands in interactive mode::

    ASCIIGenome <options>
    [h] for help: -h

And to see the help for a given command use :code:`cmd_name -h`, for example::

    [h] for help: next -h

    next [-start] [track_id]
          Move to the next feature on track_id on *current* chromosome. `next` centers
          the window on the found feature and zooms out. This is useful for quickly browsing
          through annotation files of genes or ChIP-Seq peaks in combination with read
          coverage tracks (bigwig, tdf, etc.). The `-start` flag sets the window right
          at the start of the feature, without centering and zooming out.
          ...      
    
See also :ref:`Supported_input_and_output`.

Interactive commands
--------------------

As there is no GUI, everything is handled thorough command line. Once *ASCIIGenome* is started type
a command and press ENTER to execute.

Some features of Unix console are enabled: 

* Arrow keys UP and DOWN scroll previous commands.
* TAB auto-completes commands.
* ENTER without any argument repeats the previous command.

Examples::

    [h] for help: ff <ENTER>   # Move forward
    [h] for help: <ENTER>      # Move forward again...
    [h] for help: <ENTER>      # ... and again
    [h] for help: col <TAB>    # Is expanded to colorTrack
    [h] for help: <ARROW UP>   # Shows previous command
    [h] for help: h <ENTER>    # Show help.

When track names are passed as arguments, it is not necessary to give the full name as
partial matching is enabled. This is handy since track names have an ID appended as suffix which can
be used in place of the full name, e.g. `next myLongfileName.bed#1` can be also typed as `next #1`.

These are just some functionalities to give an idea behind *ASCIIGenome*. See :ref:`Supported_input_and_output2` for 
the individual commands available.

Tips gotchas and miscellanea
----------------------------

* **Regular expression** *ASCIIGenome* makes extensive use of regular expressions. 
  Most commands use regular expression in *case sensitive* mode. 
  Use the :code:`(?i)` modifier to match in case insensitve mode, e.g. '(?i)bam' to capture 
  'foo.bam' and 'foo.BAM'. Note that the command :code:`seqRegex` by default is case insensitive,
  unless the flag :code:`-c` is set.

* When displaying bam files, *ASCIGenome* is hardcoded to disable the coverage and read tracks if
  the window size is >100,000 bp. This is to prevent the browsing to become horribly slow. To display
  such large windows  consider bigWig or tdf file format. Also consider hiding the 
  read track if not necessary with :code:`trackHeight 0 @*.bam` (or other suitable regex).

* When opening bam files, the first chromosome is often the mitochondrial chromosome chrM (or chrMT) which
  often has very high read depth (say 10,000x). This can make the opening slow. Consider using the :code:`-r`
  option in these cases. E.g. :code:`ASCIIGenome -r chr1 file1.bam file2.bam ...`
  
* If the background colour of your terminal is not white, after exiting ASCIIGenome you will have a rather awkward mix 
  of white background and user's colours. To get rid of this mix you can issue the Unix command `clear` after exiting ASCIIGenome.
