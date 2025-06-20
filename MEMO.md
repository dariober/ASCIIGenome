<!-- vim-markdown-toc GFM -->

* [Notes on maintenance and development](#notes-on-maintenance-and-development)
    * [Code style](#code-style)
    * [Release new version](#release-new-version)
        * [Upload release to github](#upload-release-to-github)
        * [Upload ASCIIGenome-x.y.z.zip to github](#upload-asciigenome-xyzzip-to-github)
        * [Update brew formula](#update-brew-formula)
    * [Test brew installation](#test-brew-installation)
        * [Update bioconda](#update-bioconda)
    * [Start new development branch](#start-new-development-branch)
        * [Create a new branch:](#create-a-new-branch)
* [Install or update gradle](#install-or-update-gradle)

<!-- vim-markdown-toc -->

Notes on maintenance and development
====================================

These notes are howto's and they are not part of the documentation. They are
here just for pro memoria. A lot of references are not general but depend on the
system. For authoritative references see github, svn, eclipse, ..., docs. 

Code style
----------

Reformat code *in place*:

```
curl -O -L https://github.com/google/google-java-format/releases/download/v1.22.0/google-java-format-1.22.0-all-deps.jar
java -jar google-java-format-1.22.0-all-deps.jar -i `find src/ -name '*.java'`
```

Release new version
-------------------

Steps to follow once you are happy with the changes to a development branch and
you want to release it as a new version.

* Make sure the version set in `ArgParse.VERSION` is bumped as appropriate.

* Test and build jar file. All tests should PASS.

* Reformat code (see above)

```
./gradlew clean
./gradlew build 
``` 

### Upload release to github

We need to create a zip file containing the jar and the helper bash script. This 
is what users will download and use.

* Prepare zip

```
cd ~/git_repos/ASCIIGenome ## Or wherever the latest local dir is

VERSION='1.12.0' # To match ArgParse.VERSION

mkdir ASCIIGenome-${VERSION}

## Copy helper script and jar file to future zip dir
cp ASCIIGenome ASCIIGenome-${VERSION}/
cp build/libs//ASCIIGenome.jar ASCIIGenome-${VERSION}/
cp INSTALL.md ASCIIGenome-${VERSION}/

## Zip up
zip -r ASCIIGenome-${VERSION}.zip ASCIIGenome-${VERSION}
rm -r ASCIIGenome-${VERSION}
```

### Upload ASCIIGenome-x.y.z.zip to github 

* Create a new release (*Draft new release*). Format of the name must be v*X.Y.Z*
  e.g. *v1.2.3*. As always, X.Y.Z must match throughout.

* Upload (e.g. using web interface).

* Add some release notes, typically copied from `CHANGELOG.md` file.

### Update brew formula 

* Edit `install/brew/asciigenome.rb` to change **release version** and **sha sum** of the zip file.

Get sha256 sum with:

```
shasum -a 256 ASCIIGenome-${VERSION}.zip
vi install/brew/asciigenome.rb ## Edit version and sha
```

* Add, commit and push edits:

```
git status
git add
git commit -m 'Update brew'
git push
```

* Merge branch to master: 

```
git merge master              # Resolve conflict
git checkout master           # Switch to master  
git merge --no-ff <my-branch> # Merge branch into master
git commit -m 'Merge to master'
git push
```

Test brew installation 
----------------------

```
brew install https://raw.githubusercontent.com/dariober/ASCIIGenome/master/install/brew/asciigenome.rb
# Or more likely:
brew upgrade https://raw.githubusercontent.com/dariober/ASCIIGenome/master/install/brew/asciigenome.rb

/usr/local/Cellar/asciigenome/${VERSION}/bin/ASCIIGenome # Full path to make sure you use the brew version
```

### Update bioconda

Similar to brew: edit [meta.yaml](https://github.com/bioconda/bioconda-recipes/blob/master/recipes/asciigenome/meta.yaml) 
as appropriate. NB: You should include sha sum here as well!

Start new development branch
----------------------------

You have uploaded a new release and now you want to develop new features. 
For development, create a branch from the master. Work on it and once happy return to
the steps to [release a new version](#release-new-version), in an iterative way.

### Create a new branch:

```
cd ~/git_repos/ASCIIGenome # Or wherever you the repo
git checkout -b v1.18.0    # Create new branch
git checkout v1.18.0       # Switch to new branch 
git push -u origin v1.18.0 # Add branch to remote     
```

# Install or update gradle

Assuming you don't have gradle already installed, in a conda env:

```
mamba install gradle
```

Create or update the gradlew wrapper:

```
gradle wrapper
```
