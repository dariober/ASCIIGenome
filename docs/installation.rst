Installation
============

Quick start 
------------------------

In the commands below replace version number with the latest from `releases <https://github.com/dariober/ASCIIGenome/releases>`_::

    wget https://github.com/dariober/ASCIIGenome/releases/download/v0.x.0/ASCIIGenome-x.y.z.zip
    unzip ASCIIGenome-x.y.z.zip

    cd ASCIIGenome-x.y.z/
    chmod a+x ASCIIGenome
    cp ASCIIGenome.jar /usr/local/bin/ # Or else in your PATH e.g. ~/bin/
    cp ASCIIGenome /usr/local/bin/     # Or else in your PATH e.g. ~/bin/

With Homebrew
------------------------------

ASCIIGenome can also be installed through `brew <http://brew.sh>`_ / `Linux Brew <https://github.com/Linuxbrew/brew>`_, although it is still not an official package::

    brew install https://raw.githubusercontent.com/dariober/ASCIIGenome/master/install/brew/asciigenome.rb

A little more detail
--------------------

:code:`ASCIIGenome.jar` requires **Java 1.7+** and this should be the only requirement. There is virtually no installation needed as :code:`ASCIIGenome` is pure Java and should work on most (all?) platforms. Download the zip file :code:`ASCIIGenome-x.x.x.zip` from `releases <https://github.com/dariober/ASCIIGenome/releases>`_, unzip it and execute the jar file with::

    java -jar /path/to/ASCIIGenome.jar --help

To avoid typing :code:`java -jar ...` every time, you can put both the helper 
script :code:`ASCIIGenome` and the jar file :code:`ASCIIGenome.jar` in the same directory in your :code:`PATH` and execute with::

    ASCIIGenome [options]

Note the helper is a bash script. To set the amount of memory available to java use the :code:`-Xmx` option as e.g. :code:`java -Xmx1500m -jar ...`.

If for some reason the text formatting misbehaves, disable it with the :code:`-nf` option.

