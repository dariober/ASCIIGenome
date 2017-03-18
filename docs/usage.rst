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