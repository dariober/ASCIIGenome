Text Only Genome Viewer!
========================


<!-- 
MEMO: Compile, package and upload to github releases
- Write-out jar from Eclipse
cd ~/svn_git/ASCIIGenome/trunk
mkdir ASCIIGenome-1.0.0 # This should match the version in ArgParse
cp ASCIIGenome ASCIIGenome-1.0.0/
cp /Users/berald01/Dropbox/Public/ASCIIGenome.jar ASCIIGenome-1.0.0/
zip -r ASCIIGenome-1.0.0.zip ASCIIGenome-1.0.0
rm -r ASCIIGenome-1.0.0

// Upload ASCIIGenome-x.y.z.zip to github releases and delete

// Update brew formula 
// Edit install/brew/asciigenome.rb to change release version and sha sum.
shasum -a 256 ASCIIGenome-1.0.0.zip
-->

Description
-----------

`ASCIIGenome` is a genome browser based on command line interface and designed for running from console terminals.

Since `ASCIIGenome` does not require a graphical interface it is particularly
useful for  quickly visualizing genomic data on remote servers while offering flexibility similar to popular GUI viewers like [IGV](https://www.broadinstitute.org/igv/).

**Documentation** is at [readthedocs/asciigenome](http://asciigenome.readthedocs.io/en/latest/).

**Support**: Bugs, comments and issues can be reported here on [GitHub](https://github.com/dariober/ASCIIGenome/issues) or on [Biostars.org](https://www.biostars.org/).

Some key features:

* Command line input and interaction, no graphical interface, minimal installation and requirements.
* Can load multiple files in various formats.
* Can access remote files via URL or ftp address.
* Easy navigation and searching for features and sequence motifs and filtering options
* Support for BS-Seq alignment

<img src="docs/screenshots/composite.png" width="800">

