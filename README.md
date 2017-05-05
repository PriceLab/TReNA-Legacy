# TReNA
 Fit transcriptional regulatory networks using gene expression, priors, machine learning

## Getting Started 

To build and test:

 - clone this repository
 - install R 3.2.3 or later; RUnit 0.4.31 or later (see below)
 - install the following solver packages:
   - glmnet R package 2.0.3 or later
   - randomForest
   - vbsr
   - flare
   - lassopv
 - cd TReNA
 - R CMD INSTALL .
 
The most reliable way to install package dependencies (and other of their dependencies):

````
source("http://bioconductor.org/biocLite.R")
biocLite(c("glmnet", "RUnit"))
````

## Using TReNA

 - open an R session
 - source("inst/unitTests/test_TReNA.R")
 - runTests()

The unitTests perform double duty: they ensure the package performs as (currently) expected;
they introduce the package to the user and developer.
Thus [test_TReNA.R](https://github.com/PriceLab/TReNA/blob/master/inst/unitTests/test_TReNA.R)
is one entry point into this project.

[Jupyter Notebook](http://nbviewer.jupyter.org/github/PriceLab/TReNA/blob/master/inst/demos/Assess_Distributions.ipynb)
