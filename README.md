# Nek5000 Documentation

This is the repository of the [Nek5000](http://nek5000.mcs.anl.gov/) documentation written using
the [Sphinx](http://www.sphinx-doc.org/) documentation framework.

## How to build

* 'make html' builds the user documentation as a set of interlinked HTML 
  and image files.  The top-level webpage is build/html/index.html.
  Supporting HTML and image files are also in the build directory.
  
Note: This requires the [Sphinx](https://pypi.python.org/pypi/Sphinx) and
[sphinx_rtd_theme](https://pypi.python.org/pypi/sphinx_rtd_theme) Python packages.  Both are
available from the [Python Package Index](http://www.sphinx-doc.org://pypi.python.org/pypi).  

## How to contribute

Please create a fork of the repository and make pull/merge requests. Keep in 
mind that the number of binary files should be kept minimal. The Makefile should be 
adapted to any special build requirements.

New issues or requests are welcome to be reported.

## How to publish on GitHub Pages

To update the GitHub Page, a contributor must have write permissions to the main [NekDoc][main-repo] repository.  
The GitHub Page should *not* contain any edits that are newer than the [master branch][master] 
of the main repository.  

Workflow:

1. Checkout the latest master
2. run `make gh-pages`

[gh-page]:   https://nek5000.github.io/NekDoc/Nek_users.html "Nek5000 user documentation on GitHub Pages"
[main-repo]: https://github.com/Nek5000/NekDoc "NekDoc repository"
[master]:    https://github.com/Nek5000/NekDoc/tree/master "NekDoc master branch"
