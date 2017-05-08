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
    load(system.file(package="TReNA", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
    var.filter <- VarianceFilter(mtx.assay = mtx.sub)
    tf.list <- getCandidates(var.filter, extraArgs = list("target.gene" = "MEF2C",
                                                          "var.size" = 0.5))
    
    checkTrue(length(tf.list$tfs) == 1)

} # test_VarianceFilter
#----------------------------------------------------------------------------------------------------
if(!interactive()) runTests()
