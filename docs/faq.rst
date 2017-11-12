FAQ and miscellanea
===================

Can I change colour theme / Can I use my configuration?
-------------------------------------------------------

To change colour theme use the ``-c/--config`` command line option or the
interactive command ``setConfig``. To make your own theme and set it as default,
use as template one of the files in the repository directory `config
<https://github.com/dariober/ASCIIGenome/blob/master/resources/config/>`_, edit
it as desired and save it as ``~/.asciigenome_config``.

Examples::

    ASCIIGenome -c metal ...    <- Use the "metal" built-in them
    ASCIIGenome -c mytheme.conf <- Read configuration from this file
    ASCIIGenome ...             <- No args to -c: Read file ~/.asciigenome_config or use default theme

For available colour names see the help in ``colorTrack`` or this `cheat sheet <http://jonasjacek.github.io/colors/>`_.

Can I turn off case sensitivity?
--------------------------------

For command that do not explicitly enable turning on or off case sensitivity,
you can prepend ``(?i)`` to your regex to match in case insensitve
mode, e.g. ``(?i)bam`` will capture  ``foo.bam`` and ``foo.BAM``. This is standard regular expression
syntax unrelated to ASCIIGenome.

Note that the command :code:`seqRegex` by default is case insensitive, unless
the flag :code:`-c` is set.

Why reads and coverage in bam tracks sometimes disappear?
---------------------------------------------------------

When displaying bam files, *ASCIIGenome* is configured to disable the coverage and read tracks if
the window size is >100,000 bp. This is to prevent the browsing to become too slow. To display
such large windows  consider bigWig or tdf file format. 

Can I execute multiple commands inside ``-x/--exec`` or at the command prompt?
--------------------------------------------------------------------------------------

Use the ``&&`` to concatenate commands (similar to Unix syntax).  E.g.
``colorTrack red && goto chr1:150000 && zo``.

How can I print the header of a VCF file?
-----------------------------------------

Assuming :code:`bcftools` is available, use the :code:`sys` command, for example::

	sys bcftools view -H my_variants.vcf.gz | less

Of course you can further parse the header by piping to standard Unix tools. For 
example, to exclude the ``contig`` lines use::
	
	sys bcftools view -h mutect/WW00282.vcf.gz | grep -v '##contig' | less
 