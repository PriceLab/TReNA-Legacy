library(TReNA)
library(RUnit)
#----------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_NullFilter()
   
} # runTests
#----------------------------------------------------------------------------------------------------
test_NullFilter <- function()
{
    printf("--- test_NullFilter")

    # Check that dummy data returns all gene names
    load(system.file(package="TReNA", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
    null.filter <- NullFilter(mtx.assay = mtx.sub)
    checkEquals(getCandidates(null.filter),rownames(mtx.sub))

} # test_NullFilter
#----------------------------------------------------------------------------------------------------
if(!interactive()) runTests()
