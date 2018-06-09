Notes on maintenance and development
====================================

These notes are howto's and they are not part of the documentation. They are
here just for pro memoria. A lot of references are not general but depend on the
system. For authoritative references see github, svn, eclipse, ..., docs. 

Release new version
-------------------

Steps to follow once you are happy with the changes to a development branch and
you want to release it as a new version.

* Make sure the version set in `ArgParse.VERSION` is bumped as appropriate.

* Test and build jar file. All tests should PASS.

```
./gradlew clean
./gradlew build 
``` 

### Upload release to github

We need to create a zip file containing the jar and the helper bash script. This 
is what users will download and use.

* From eclipse, write out the definitive jar file. 

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
git checkout -b v1.13.0    # Create new branch
git checkout v1.13.0       # Switch to new branch 
git push -u origin v1.13.0 # Add branch to remote     
```

### Set up Eclipse project

* In Eclipse: *File* `->` *New* `->` *Java project*

* *Project name*: something meaningful, it doesn't really matter. Use the branch name maybe.

* Uncheck *Use default location* and browse instead to the new branch directory, for 
example it may be `/Users/berald01/git_repos/ASCIIGenome/branches/argparse`.

* If not selected, choose execution environment Java 1.7. Other options should be fine as default.

* You may need to edit the classpath if your directory structure changed since last time. 
Use *Build Path* `->` *Configure Build Path* to set up the appropriate packages.
If Eclipse doesn't complain about missing imports, than the classpath is fine. 

That's it, start developing in Eclipse.
