ASCIIGenome: Text only genome viewer!
=====================================

For full documentation, including installation please see
http://asciigenome.readthedocs.io/en/latest/index.html

Installation
------------

Download the latest version of ASCIIGenome as a zip file from
https://github.com/dariober/ASCIIGenome/releases and unzip the file.  (you may
already have done that). 

Than simply copy the `ASCIIGenome` script and the `ASCIIGenome.jar` file to a
convenient directory on your `PATH`. Remember to make the `ASCIIGenome` script
executable.

Fo example, replace here x.y.z with the latest version:

``` 
wget https://github.com/dariober/ASCIIGenome/releases/download/vX.Y.Z/ASCIIGenome-x.y.z.zip
unzip ASCIIGenome-x.y.z.zip

cd ASCIIGenome-x.y.z/ 
chmod a+x ASCIIGenome # Make script executable

cp ASCIIGenome.jar /usr/local/bin/ # Or else in your PATH e.g. ~/bin/ 
cp ASCIIGenome /usr/local/bin/     # Or else in your PATH e.g. ~/bin/
```

Minimal example
---------------

Now you should be all set to go, e.g.:

```
ASCIIGenome ftp://ftp.ensembl.org/pub/release-89/gff3/homo_sapiens/Homo_sapiens.GRCh38.89.chromosome.Y.gff3.gz
```

At the ASCIIGenome command prompt type `h` for general help or `<cmd> -h` for help on command `<cmd>`.
