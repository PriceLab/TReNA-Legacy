library(TReNA)
library(RUnit)
#----------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_CandidateFilter()
   
} # runTests
#----------------------------------------------------------------------------------------------------
test_CandidateFilter <- function()
{
    printf("--- test_CandidateFilter")

    # Create a filter with an empty matrix
    candidate.filter <- CandidateFilter(matrix())
    checkEquals(is(candidate.filter), "CandidateFilter")

} #test_CandidateFilter
#----------------------------------------------------------------------------------------------------
if(!interactive()) runTests()
