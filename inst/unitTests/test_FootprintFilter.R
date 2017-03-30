library(TReNA)
library(RUnit)
#----------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_FootprintFilter()
   
} # runTests
#----------------------------------------------------------------------------------------------------
test_FootprintFilter <- function()
{
    printf("--- test_FootprintFilter")

    # Load ampAD data and filter it based on footprints
    load(system.file(package="TReNA", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
    footprint.filter <- FootprintFilter(mtx.assay = mtx.sub)

    db.address <- system.file(package="TReNA", "extdata")
    genome.db.uri <- paste("sqlite:/",db.address,"genome.sub.db", sep = "/")
    project.db.uri <- paste("sqlite:/",db.address,"project.sub.db", sep = "/")
    target.gene <- "MEF2C"
    list.out <- getCandidates(footprint.filter, extraArgs = list("target.gene" = target.gene,
                                                            "genome.db.uri" = genome.db.uri,
                                                            "project.db.uri" = project.db.uri,
                                                            "size.upstream" = 1000,
                                                            "size.downstream" = 1000))

    # Make sure it grabs the right number of genes
    checkEquals(length(list.out[["tfs"]]), 64)

} # test_FootprintFilter
#----------------------------------------------------------------------------------------------------
if(!interactive()) runTests()
