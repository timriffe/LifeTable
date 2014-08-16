LifeTable
=========

This package includes a lifetable function, ```LT()```, with optional mx and ax smoothing, different methods 
of estimating ax, and options for single and abridged lifetables.

*this package is soon to undergo an overhaul.

Installation
============

To download the most recent version of LifeTable:

1. make sure you have the most recent version of R
2. look at the notes below if you're on Windows or Mac.

Download the [zip ball](https://github.com/timriffe/LifeTable/zipball/master) or [tar ball](https://github.com/timriffe/LifeTable/tarball/master), decompress and run `R CMD INSTALL` on it in the terminal command line, or use the **devtools** package to install the development version:

```r
# install.packages("devtools")

library(devtools)
install_github("timriffe/LifeTable", subdir = "LifeTable")
```

**Note**: Windows users need [Rtools](http://cran.r-project.org/bin/windows/Rtools/) to install from github as shown above. Get the most recent version of [R for Windows](http://cran.r-project.org/bin/windows/base/) and download and install the version of Rtools that corresponds to it.

**Note**: Mac users might be required to install the appropriate version [XTools](https://developer.apple.com/xcode/) from the [Apple Developer site](https://developer.apple.com/) in order to install the development version.  You may need to [register as an Apple developer](https://developer.apple.com/programs/register/).  An older version of XTools may also be required.

Help
===============
All functions are documented in the standard way, which means that once you load the package using ```library(LifeTable)```
you can just type for example ```?LT``` to see the help file. Help files are also available in html format, 
[here](http://timriffe.github.io/LifeTable/help/).

To report a bug
===============
Just go to the [main repository page](https://github.com/timriffe/LifeTable) and click on the ```Issues``` 
button on the right side. That's a convenient way to track bugs. Otherwise, just email the maintainer. Feature 
requests can also be made to the maintainer. Motivated individuals are also free to offer assistance by collaborating 
via the ```git``` version control system and github.

