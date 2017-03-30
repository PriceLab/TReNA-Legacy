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
    x <- test_developAndFitDummyTestData(quiet=TRUE)
    null.filter <- NullFilter(mtx.assay = x$assay)
    checkEquals(getCandidates(null.filter),rownames(x$assay))

} # test_NullFilter
#----------------------------------------------------------------------------------------------------
if(!interactive()) runTests()
