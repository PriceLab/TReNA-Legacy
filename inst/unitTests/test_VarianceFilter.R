library(TReNA)
library(RUnit)
#----------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
    test_VarianceFilter()

} # runTests
#----------------------------------------------------------------------------------------------------
test_VarianceFilter <- function()
{
    printf("--- test_VarianceFilter")

    # Create dummy data and filter it based on variance
    x <- test_developAndFitDummyTestData(quiet=TRUE)
    var.filter <- VarianceFilter(mtx.assay = x$assay)
    tf.list <- getCandidates(var.filter, extraArgs = list("target.gene" = x$correlated.target,
                                                          "var.size" = 0.5))
    
    checkTrue(length(tf.list$tfs) > 0)

} # test_VarianceFilter
#----------------------------------------------------------------------------------------------------
if(!interactive()) runTests()
