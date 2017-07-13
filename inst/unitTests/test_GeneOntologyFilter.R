library(TReNA)
library(RUnit)
#----------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
    test_basicUse()

} # runTests
#----------------------------------------------------------------------------------------------------
test_basicUse <- function()
{
    printf("--- test_directMode")

    require(org.Hs.eg.db)

    #  GO:0010467                        gene expression
    #  GO:0097659   nucleic acid-templated transcription
    #  GO:0001172           transcription, RNA-templated
    #  GO:0006351           transcription, DNA-templated

    goFilter <- GeneOntologyFilter(org.Hs.eg.db, GOTerm="GO:0006351")
    candidates <- getCandidates(goFilter)
    checkTrue(all(c("tbl", "tfs") %in% names(candidates)))
    first.genes <- head(candidates$tfs)
       # simple test: they should be sorted, and all start with "A"
    checkEquals(length(grep("^A", first.genes)), 6)

} # test_basicUse
#----------------------------------------------------------------------------------------------------
if(!interactive()) runTests()
