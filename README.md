Text Only Genome Viewer!
========================

- [Description](#description)
- [Requirements and Installation](#requirements-and-installation)
    - [Installation quick start](#installation-quick-start)
    - [Installation through Homebrew](#installation-through-homebrew)
    - [A little more detail](#a-little-more-detail)
- [Usage examples](#usage-examples)
    - [Minimal example](#minimal-example)
    - [Open and browse](#open-and-browse)
    - [Finding & filtering stuff](#finding--filtering-stuff)
    - [Chaining commands](#chaining-commands)
    - [Batch processing](#batch-processing)
- [Supported input](#supported-input)
- [Genome option](#genome-option)
- [Formatting of reads and features](#formatting-of-reads-and-features)
- [Saving screenshots](#saving-screenshots)
- [Tips gotchas and miscellanea](#tips-gotchas-and-miscellanea)
- [Interactive commands](#interactive-commands)
- [Credits](#credits)


<!-- 
MEMO: Compile, package and upload to github releases
- Write-out jar from Eclipse
cd ~/svn_git/ASCIIGenome/trunk
mkdir ASCIIGenome-0.4.0 # This should match the version in ArgParse
cp ASCIIGenome ASCIIGenome-0.4.0/
cp /Users/berald01/Dropbox/Public/ASCIIGenome.jar ASCIIGenome-0.4.0/
zip -r ASCIIGenome-0.4.0.zip ASCIIGenome-0.4.0
rm -r ASCIIGenome-0.4.0

// Upload ASCIIGenome-0.4.0.zip to github releases and delete

// Update brew formula 
// Edit install/brew/asciigenome.rb to change release version and sha sum.
shasum -a 256 ASCIIGenome-0.4.0.zip
-->

Description
===========

`ASCIIGenome` is a command-line genome browser running from terminal window and solely based on
ASCII characters. Since `ASCIIGenome` does not require a graphical interface it is particularly
useful for  quickly visualizing genomic data on remote servers. The idea is to make `ASCIIGenome`
the Vim  of genome viewers.

As far as I know, the closest program to `ASCIIGenome` is [samtools tview](http://samtools.sourceforge.net/tview.shtml) but 
`ASCIIGenome` offers much more flexibility, similar to popular GUI viewers like [IGV](https://www.broadinstitute.org/igv/).

Some key features:

* Command line input and interaction, no graphical interface, minimal installation and requirements.
* Can load multiple files in various formats.
* Can access remote files via URL or ftp address.
* Easy navigation and searching for features and sequence motifs and filtering options
* Support for BS-Seq alignment

<img src="screenshots/ex3.png" width="800">

