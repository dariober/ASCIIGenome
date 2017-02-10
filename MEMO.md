<!-- MarkdownTOC -->

- [Notes on maintenance and development](#notes-on-maintenance-and-development)
    - [Fix bugs and add features](#fix-bugs-and-add-features)
    - [Release new version](#release-new-version)
        - [Upload to github](#upload-to-github)
        - [Upload ASCIIGenome-x.y.z.zip to github](#upload-asciigenome-xyzzip-to-github)
        - [Update brew formula](#update-brew-formula)
        - [Update bioconda](#update-bioconda)

<!-- /MarkdownTOC -->


Notes on maintenance and development
====================================

These notes are howto's and they are not part of the documentation. They are
here just for pro memoria. A lot of references are not general but depend on the
system.

Fix bugs and add features
-------------------------

* 

Release new version
-------------------

Steps to follow once you are happy with the changes to a development branch and
you want to release it as a new version.

* Ensure all tests succeed. In Eclipse, in panel *Package explorer*, right-click on the test
  source directory, then launch `Run As -> JUnit Test`. This will run all the tests
  under `test`, which should mean all of them.

* Make sure the version set in `ArgParse.VERSION` is correct.

* Merge branch to master: 

### Upload to github

We need to create a zip file containing the jar and the helper bash script. This 
is what users will download and use.

* From eclipse, write out the definitive jar file. 

* Prepare zip

```
cd ~/svn_git/ASCIIGenome/trunk ## Or wherever the latest local dir is

mkdir ASCIIGenome-1.0.0        ## The x.y.z tag should match the version in ArgParse.VERSION

## Copy helper script and jar file to future zip dir
cp ASCIIGenome ASCIIGenome-1.0.0/
cp /Users/berald01/Dropbox/Public/ASCIIGenome.jar ASCIIGenome-1.0.0/

## Zip up
zip -r ASCIIGenome-1.0.0.zip ASCIIGenome-1.0.0
rm -r ASCIIGenome-1.0.0
```

### Upload ASCIIGenome-x.y.z.zip to github 

* Create a new release (*Draft new release*). Format of the name must be v*X.Y.Z*
  e.g. *v1.2.3*. As always, X.Y.Z must match throughout.

* Upload (e.g. using web interface).

* Add some release notes, typically copied from `CHANGELOG.md` file.

### Update brew formula 

Edit `install/brew/asciigenome.rb` to change **release version** and **sha sum** of the zip file.

Get sha256 sum with:

```
shasum -a 256 ASCIIGenome-x.y.z.zip
```

### Update bioconda

Similar to brew: edit [meta.yaml](https://github.com/bioconda/bioconda-recipes/blob/master/recipes/asciigenome/meta.yaml) 
as appropriate. NB: You should include sha sum here as well!


