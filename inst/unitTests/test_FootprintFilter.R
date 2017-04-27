library(TReNA)
library(RUnit)
#----------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_FootprintFilter.geneCentered()
   test_FootprintFilter.byRegion()
   test_FootprintFilter.byTwoRegions()
   test_FootprintFilter.geneCenteredPlusRegion()

} # runTests
#----------------------------------------------------------------------------------------------------
test_FootprintFilter.geneCentered <- function()
{
   printf("--- test_FootprintFilter.geneCentered")

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
# should exactly replicate the geneCentered results above
test_FootprintFilter.byRegion <- function()
{
   printf("--- test_FootprintFilter.byRegion")

   db.address <- system.file(package="TReNA", "extdata")
   genome.db.uri    <- paste("sqlite:/", db.address, "genome.sub.db",  sep = "/")
   footprint.db.uri <- paste("sqlite:/", db.address, "project.sub.db", sep = "/")
   #target.gene <- "MEF2C"

   mef2c.tss <- 88904257   # minus strand
   chrom <- "chr5"
   start <- mef2c.tss - 1000
   end   <- mef2c.tss + 1000
   chromLoc <- sprintf("%s:%d-%d", chrom, start, end)

   recipe <- list(genomeDB=genome.db.uri,
                  footprintDB=footprint.db.uri,
                  geneCenteredSpec=list(),
                  regionsSpec=c(chromLoc))

   filter <- FootprintFilter(recipe$genomeDB,
                             recipe$footprintDB,
                             recipe$geneCenteredSpec,
                             recipe$regionsSpec,
                             quiet=TRUE)

   list.out <- getCandidates(filter)
   min.fp.start <- min(list.out$tbl$start)
   max.fp.end   <- max(list.out$tbl$end)
   checkTrue(min.fp.start > (mef2c.tss - 1000))
   checkTrue(max.fp.end <  (mef2c.tss + 1000))
   span <- max.fp.end - min.fp.start   # less than 2000, more than nothing
   checkTrue(span > 100)
   checkTrue(span <= 20001)

   checkEquals(length(list.out$tfs), 185)

} # test_FootprintFilter.byRegion
#----------------------------------------------------------------------------------------------------
# get footprints for tss +/- 1000, by specifying two contiguous regions, one 500 bp, one 1500 bp
test_FootprintFilter.byTwoRegions <- function()
{
   printf("--- test_FootprintFilter.byTwoRegions")

   db.address <- system.file(package="TReNA", "extdata")
   genome.db.uri    <- paste("sqlite:/", db.address, "genome.sub.db",  sep = "/")
   footprint.db.uri <- paste("sqlite:/", db.address, "project.sub.db", sep = "/")

   mef2c.tss <- 88904257   # minus strand
   chrom <- "chr5"
   start <- mef2c.tss - 1000
   end   <- mef2c.tss - 500
   chromLoc.1 <- sprintf("%s:%d-%d", chrom, start, end)

   start <- mef2c.tss - 500
   end   <- mef2c.tss + 1000
   chromLoc.2 <- sprintf("%s:%d-%d", chrom, start, end)

   recipe <- list(genomeDB=genome.db.uri,
                  footprintDB=footprint.db.uri,
                  geneCenteredSpec=list(),
                  regionsSpec=c(chromLoc.1, chromLoc.2))

   filter <- FootprintFilter(recipe$genomeDB,
                             recipe$footprintDB,
                             recipe$geneCenteredSpec,
                             recipe$regionsSpec,
                             quiet=TRUE)

   list.out <- getCandidates(filter)
   min.fp.start <- min(list.out$tbl$start)
   max.fp.end   <- max(list.out$tbl$end)
   checkTrue(min.fp.start > (mef2c.tss - 1000))
   checkTrue(max.fp.end <  (mef2c.tss + 1000))
   span <- max.fp.end - min.fp.start   # less than 2000, more than nothing
   checkTrue(span > 100)
   checkTrue(span <= 20001)

   checkEquals(length(list.out$tfs), 185)
   checkTrue(nrow(list.out$tbl) > 20)
   checkTrue(nrow(list.out$tbl) < 25)

} # test_FootprintFilter.byTwoRegions
#----------------------------------------------------------------------------------------------------
test_FootprintFilter.geneCenteredPlusRegion <- function()
{
   printf("--- test_FootprintFilter.geneCenteredPlusRegion")

   db.address <- system.file(package="TReNA", "extdata")
   genome.db.uri    <- paste("sqlite:/", db.address, "genome.sub.db",  sep = "/")
   footprint.db.uri <- paste("sqlite:/", db.address, "project.sub.db", sep = "/")
   target.gene <- "MEF2C"
   mef2c.tss <- 88904257   # minus strand

   extra.region <- sprintf("chr5:%d-%d", mef2c.tss+1000, mef2c.tss+2000)

   recipe <- list(genomeDB=genome.db.uri,
                  footprintDB=footprint.db.uri,
                  geneCenteredSpec=list(targetGene=target.gene,
                                        tssUpstream=1000,
                                        tssDownstream=1000),
                  regionsSpec=c(extra.region))

   filter <- FootprintFilter(recipe$genomeDB,
                             recipe$footprintDB,
                             recipe$geneCenteredSpec,
                             recipe$regionsSpec,
                             quiet=TRUE)

   list.out <- getCandidates(filter)
   min.fp.start <- min(list.out$tbl$start)
   max.fp.end   <- max(list.out$tbl$end)
   checkTrue(min.fp.start > (mef2c.tss - 1000))
   checkTrue(max.fp.end <  (mef2c.tss + 2000))
   span <- max.fp.end - min.fp.start   # less than 2000, more than nothing
   checkTrue(span > 2000)
   checkTrue(span <= 3000)

   checkTrue(nrow(list.out$tbl) > 25)   # we see 22 rows without the extra 1k bp

} # test_FootprintFilter.geneCenteredPlusRegion
#----------------------------------------------------------------------------------------------------
if(!interactive()) runTests()
