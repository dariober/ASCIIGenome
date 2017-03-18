FAQ and miscellanea
===================

* **The default colour theme**, *black on white*, **is annoying on my dark background 
  terminal. How do I changed it?** 

  To change colour theme use the ``-c/--config`` command line option or the 
  interactive command ``setConfig``. To make your own theme and set it as default,
  use as template one of the files in the repository directory `config <https://github.com/dariober/ASCIIGenome/blob/master/resources/config/>`_, edit it as desired and save it
  as ``~/.asciigenome_config``.

  For available colours execute see the help in ``trackColor``.

|

* **Can I turn off case sensitivity in ASCIIGenome's commands?**

  For command that do not explicitly enable turning on or off case sensitivity,
  you can prepend ``(?i)`` to your regex to match in case insensitve
  mode, e.g. '(?i)bam' will capture  'foo.bam' and 'foo.BAM'. This is standard regular expression
  syntax unrelated to ASCIIGenome.

  Note that the command :code:`seqRegex` by default is case insensitive, unless
  the flag :code:`-c` is set.

|

* **Why read and coverage tracks in bam files sometimes disappear?**

  When displaying bam files, *ASCIGenome* is hardcoded to disable the coverage and read tracks if
  the window size is >100,000 bp. This is to prevent the browsing to become horribly slow. To display
  such large windows  consider bigWig or tdf file format. Also consider hiding the 
  read track if not necessary with :code:`trackHeight 0 @*.bam` (or other suitable regex).

| 

* **How do I execute multiple commands inside the** ``-x/--exec`` **command line option
  or at the command prompt?**

  Use the ``&&`` to concatenate commands (similar to Unix syntax). 
  E.g. ``colorTrack red && goto chr1:150000 && zo``.

|

* When opening bam files, the first chromosome is often the mitochondrial chromosome chrM (or chrMT) which
  often has very high read depth (say 10,000x). This can make the opening slow. Consider using the :code:`-r`
  option in these cases. E.g. :code:`ASCIIGenome -r chr1 file1.bam file2.bam ...`
