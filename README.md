# TReNA
 Fit transcriptional regulatory networks using gene expression, priors, machine learning

To build and test:

 - clone this repository
 - install R 3.2.3 or later
 - install glmnet R package 2.0.3 or later
 - cd TReNA
 - R CMD BUILD .
 - open an R session
 - source("inst/unitTests/test_TReNA.R")
 - runTests()

The unitTests perform double duty: they ensure the package performs as (currently) expected; they introduce the package to the user and developer.  Thus [test_TReNA.R](https://github.com/PriceLab/TReNA/blob/master/inst/unitTests/test_TReNA.R) is your entry point into this project.
