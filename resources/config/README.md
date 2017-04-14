Configuration file  
===================

The configuration file contains the settings for (mostly) colours and other parameters. 

This directory contains some examples which are shipped with the 
compiled jar file `ASCIIGenome.jar`. At the time of this writing the default configuration
is the one in `metal.conf`. 

You can create you own configuration file using one of these as template. The 
structure is quite simple. There are two columns separated by one or more whitespaces.

* First column is the parameter name. Do not change, add or remove parameters, changing their order is fine.

* Second column is the parameter value. For parameters that take a colour see 
this nice [cheat sheet](https://jonasjacek.github.io/colors/) or execute `ASCIIGenome -ni -x 'colorTrack -h'`.

For actually making use of your own configuration see [FAQ](http://asciigenome.readthedocs.io/en/latest/faq.html)
or the command [setconfig](http://asciigenome.readthedocs.io/en/latest/commandHelp.html#setconfig)