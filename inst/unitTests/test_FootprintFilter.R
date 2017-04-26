library(TReNA)
library(RUnit)
#----------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_FootprintFilter.geneCentered()

} # runTests
#----------------------------------------------------------------------------------------------------
test_FootprintFilter.geneCentered <- function()
{
    printf("--- test_FootprintFilter.geneCentered")

    # Load ampAD data and filter it based on footprints
    #load(system.file(package="TReNA", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))

    db.address <- system.file(package="TReNA", "extdata")
    genome.db.uri    <- paste("sqlite:/", db.address, "genome.sub.db",  sep = "/")
    footprint.db.uri <- paste("sqlite:/", db.address, "project.sub.db", sep = "/")
    target.gene <- "MEF2C"

    recipe <- list(genomeDB=genome.db.uri,
                   footprintDB=footprint.db.uri,
                   geneCenteredSpec=list(targetGene=target.gene,
                                         tssUpstream=1000,
                                         tssDownstream=1000),
                   regionsSpec=list())

    filter <- FootprintFilter(recipe$genomeDB,
                              recipe$footprintDB,
                              recipe$geneCenteredSpec,
                              recipe$regionsSpec,
                              quiet=TRUE)

    list.out <- getCandidates(filter)
    mef2c.tss <- 88904257   # minus strand
    min.fp.start <- min(list.out$tbl$start)
    max.fp.end   <- max(list.out$tbl$end)
    checkTrue(min.fp.start > (mef2c.tss - 1000))
    checkTrue(max.fp.end <  (mef2c.tss + 1000))
    span <- max.fp.end - min.fp.start   # less than 2000, more than nothing
    checkTrue(span > 100)
    checkTrue(span <= 20001)

    checkEquals(length(list.out$tfs), 185)

} # test_FootprintFilter.geneCentered
#----------------------------------------------------------------------------------------------------
if(!interactive()) runTests()
