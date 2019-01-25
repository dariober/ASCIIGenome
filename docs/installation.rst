Installation
============

Quick start 
------------------------

Basically, download :code:`ASCIIGenome-X.Y.Z.zip` from `releases <https://github.com/dariober/ASCIIGenome/releases>`_, 
unzip, copy :code:`ASCIIGenome` and :code:`ASCIIGenome.jar` to a directory of your liking and that's it.

For example, in the commands below replace version number with the latest from `releases <https://github.com/dariober/ASCIIGenome/releases>`_::

    wget https://github.com/dariober/ASCIIGenome/releases/download/vX.Y.Z/ASCIIGenome-x.y.z.zip
    unzip ASCIIGenome-x.y.z.zip

    cd ASCIIGenome-x.y.z/
    chmod a+x ASCIIGenome
    cp ASCIIGenome.jar /usr/local/bin/ # Or else in your PATH e.g. ~/bin/
    cp ASCIIGenome /usr/local/bin/     # Or else in your PATH e.g. ~/bin/

With conda
----------

ASCIIGenome is available as a `bioconda package <https://bioconda.github.io/recipes/asciigenome/README.html>`_
and it can be installed with the conda package manager::

    conda install -c bioconda asciigenome

With Homebrew
------------------------------

ASCIIGenome can also be installed through `brew <http://brew.sh>`_ / `Linux Brew <https://github.com/Linuxbrew/brew>`_, although it is still not an official package::

    brew install https://raw.githubusercontent.com/dariober/ASCIIGenome/master/install/brew/asciigenome.rb

A little more detail
--------------------

:code:`ASCIIGenome.jar` requires **Java 1.8+** and this should be the only requirement. There is virtually no installation needed as :code:`ASCIIGenome` is pure Java and should work on most (all?) platforms. Download the zip file :code:`ASCIIGenome-x.x.x.zip` from `releases <https://github.com/dariober/ASCIIGenome/releases>`_, unzip it and execute the jar file with::

    java -jar /path/to/ASCIIGenome.jar --help

To avoid typing :code:`java -jar ...` every time, you can put both the helper 
script :code:`ASCIIGenome` and the jar file :code:`ASCIIGenome.jar` in the same directory in your :code:`PATH` and execute with::

    ASCIIGenome [options]

.. note:: This part below has nothing to do with ASCIIGenome specifically. These are just general instructions to add executable files to your PATH. 

For Unix users: If you have administrator rights and you want to make ASCIIGenome available to all users,
a popular choice of installation directory is :code:`/usr/local/bin/`, *e.g.*::

    cp ASCIIGenome.jar /usr/local/bin/
    cp ASCIIGenome /usr/local/bin/

If you don't have administrator rights (*i.e.* you get a :code:`Permission denied` error) you can instead copy to a directory that you have on your 
PATH and where you have permission to write. A popular user directory for executable files is :code:`$HOME/bin` (*e.g.* :code:`/home/myName/bin` or :code:`/Users/myName/bin` or short :code:`~/bin`). *e.g.*::

    cp ASCIIGenome.jar ~/bin/
    cp ASCIIGenome ~/bin/

If :code:`~/bin` does not exist or is not on your PATH create it with::

    mkdir ~/bin/

And to add to it to your PATH edit your profile file to add the new 
directory. *E.g.* edit :code:`~/.bash_profile` to::
    
    PATH=/home/myName/bin:$PATH

Then reload the profile file or log off and log back in to make the changes effective.

Note the helper is a bash script. To set the amount of memory available to java use the :code:`-Xmx` option as e.g. :code:`java -Xmx1500m -jar ...`.

If for some reason the text formatting misbehave, disable it with the :code:`-nf` option.

Compiling the source code
-------------------------

From version 1.13 ASCIIGenome is built using the `gradle build tool <https://gradle.org/>`_. If you want to edit the source code and 
re-compile it to an executable jar, all you need to do is::

    # Get source
    git clone https://github.com/dariober/ASCIIGenome.git
    cd ASCIIGenome
    
    ./gradlew clean
    ./gradlew build -x test

The executable jar file will be in :code:`build/libs/ASCIIGenome.jar`. 
The :code:`-x test` option builds the code without running the tests. 
