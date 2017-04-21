library(TReNA)
library(RUnit)
#----------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
    test_directMode()

} # runTests
#----------------------------------------------------------------------------------------------------
test_directMode <- function()
{
    printf("--- test_directMode")

    require(org.Hs.eg.db)
    goFilter <- GeneOntologyFilter(org.Hs.eg.db)
    candidates <- getCandidates(goFilter, list(mode="direct"))
    checkTrue(all(c("tbl", "tfs") %in% names(candidates)))
    first.genes <- head(candidates$tfs)
       # simple test: they should be sorted, and all start with "A"
    checkEquals(length(grep("^A", first.genes)), 6)

} # test_VarianceFilter
#----------------------------------------------------------------------------------------------------
if(!interactive()) runTests()
