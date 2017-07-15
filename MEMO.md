Notes on maintenance and development
====================================

These notes are howto's and they are not part of the documentation. They are
here just for pro memoria. A lot of references are not general but depend on the
system. For authoritative references see github, svn, eclipse, ..., docs. 

Release new version
-------------------

Steps to follow once you are happy with the changes to a development branch and
you want to release it as a new version.

* Ensure all tests succeed. In Eclipse, in panel *Package explorer*, right-click on the test
  source directory, then launch `Run As -> JUnit Test`. This will run all the tests
  under `test`, which should mean all of them.

* Make sure the version set in `ArgParse.VERSION` is correct.

* Commit changes to repository.

* Merge branch to master: 

### Upload to github

We need to create a zip file containing the jar and the helper bash script. This 
is what users will download and use.

* From eclipse, write out the definitive jar file. 

* Prepare zip

```
cd ~/git_repos/ASCIIGenome ## Or wherever the latest local dir is

mkdir ASCIIGenome-1.9.0        ## The 1.9.0 tag should match the version in ArgParse.VERSION

## Copy helper script and jar file to future zip dir
cp ASCIIGenome ASCIIGenome-1.9.0/
cp ~/Dropbox/Public/ASCIIGenome.jar ASCIIGenome-1.9.0/
cp INSTALL.md ASCIIGenome-1.9.0/

## Zip up
zip -r ASCIIGenome-1.9.0.zip ASCIIGenome-1.9.0
rm -r ASCIIGenome-1.9.0
```

### Upload ASCIIGenome-1.9.0.zip to github 

* Create a new release (*Draft new release*). Format of the name must be v*X.Y.Z*
  e.g. *v1.2.3*. As always, X.Y.Z must match throughout.

* Upload (e.g. using web interface).

* Add some release notes, typically copied from `CHANGELOG.md` file.

### Update brew formula 

Edit `install/brew/asciigenome.rb` to change **release version** and **sha sum** of the zip file.

Get sha256 sum with:

```
shasum -a 256 ASCIIGenome-1.9.0.zip
vi install/brew/asciigenome.rb ## Edit version and sha
```

### Update bioconda

Similar to brew: edit [meta.yaml](https://github.com/bioconda/bioconda-recipes/blob/master/recipes/asciigenome/meta.yaml) 
as appropriate. NB: You should include sha sum here as well!

### Merge branch to master

* Go to the master page at https://github.com/dariober/ASCIIGenome

* Select `New pull request`

* From the scroll down menu `compare to: ...` choose the branch you want to merge 
(typically the one you used to produce the zip file above).

* Hopefully github tells you that the branch can be autmatically merged. If so, 
just follow the `Create pull request` link. (If it cannot merged, deal with it to fugure out why...!)

Start new development branch
----------------------------

You have uploaded a new release and now you want to develop new features. 
For development, create a branch from the master. Work on it and once happy return to
the steps to [release a new version](#release-new-version), in an iterative way.

### Create a new branch:

* Go to the master page https://github.com/dariober/ASCIIGenome

* From scroll down menu `Branch: master` type in a new for the new branch and create 
the new branch.

* Check out the new branch:

```
cd /Users/berald01/git_repos/
svn co https://github.com/dariober/ASCIIGenome
```

This command effectively checks out the entire repository. If
`/Users/berald01/git_repos/` already contains dir `ASCIIGenome` with the trunk and
old branches, svn will not download them again. That's good because downloading
the entire repository takes a while.

* Get hold of test data that is not in the repository. See `README.md` in `test_data`
to download and prepare them or just copy it from one of the old branches.

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
