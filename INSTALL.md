ASCIIGenome: Text only genome viewer!
=====================================

For full documentation, including installation please see
http://asciigenome.readthedocs.io/en/latest/index.html

Installation
------------

Download the latest version of ASCIIGenome as a zip file from
https://github.com/dariober/ASCIIGenome/releases and unzip it.

Then, copy the `ASCIIGenome` script and the `ASCIIGenome.jar` file to a
convenient directory on your `PATH`. Remember to make the `ASCIIGenome` script
executable.

For example (replace x.y.z with the latest version):

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
ASCIIGenome http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/functional_annotation/filtered/ALL.chr21.phase3_shapeit2_mvncall_integrated_v5.20130502.sites.annotation.vcf.gz
```

At the ASCIIGenome command prompt type `h` for general help or `<cmd> -h` for
help on command `<cmd>` (e.g. `print -h`).
