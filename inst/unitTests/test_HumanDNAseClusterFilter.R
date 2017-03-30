library(TReNA)
library(RUnit)
#----------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_defaultConstructor()
   test_.getRegions()
   test_.fetchSequence()
   test_.getScoredMotifs()
   test_mef2cPromoter.incrementally()

   test_mef2cPromoter.normalUse()

} # runTests
#----------------------------------------------------------------------------------------------------
test_defaultConstructor <- function()
{
   printf("--- test_defaultConstructor")

   ocf <- HumanDNAseClusterFilter();
   load(system.file(package="TReNA", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))

} # test_defaultConstructor
#----------------------------------------------------------------------------------------------------
test_.getRegions <- function()
{
   printf("--- test_.getRegions");

   chrom <- "chr5"
   start <- 88819630
   end   <- 88835936

   tbl.regions <- TReNA:::.getRegions(chrom, start, end)
   checkTrue(nrow(tbl.regions) > 20)
   checkEquals(colnames(tbl.regions), c("chrom", "chromStart", "chromEnd", "score", "sourceCount"))
   checkTrue(all(tbl.regions$chrom == chrom))
   checkTrue(all(tbl.regions$chromStart >= start))
   checkTrue(all(tbl.regions$chromStart <= end))
   checkTrue(all(tbl.regions$chromEnd >= start))
   checkTrue(all(tbl.regions$chromEnd <= end))

      # a very small region, revealed by previous call, known to have no hits
   start <-  88820860
   end   <-  88820880
   tbl.regions <- TReNA:::.getRegions(chrom, start, end)
   checkEquals(nrow(tbl.regions), 0)

} # test_.getRegions
#----------------------------------------------------------------------------------------------------
test_.fetchSequence <- function()
{
   printf("--- test_.fetchSequence")
   chroms <- rep("chr5", 3)
   starts <- c(88819700, 88820700, 88820980)
   ends   <- c(88819910, 88820850, 88821130)

   tbl.regions <- data.frame(chrom=chroms, chromStart=starts, chromEnd=ends, stringsAsFactors=FALSE)
   seqs <- TReNA:::.getSequence(tbl.regions)
   expected.lengths <- 1 + ends - starts
   checkEquals(unlist(lapply(seqs, nchar)), expected.lengths)
   invisible(seqs)

} # test_.fetchSequence
#----------------------------------------------------------------------------------------------------
test_.getScoredMotifs <- function()
{
   printf("--- test_.getScoredMotifs")
   seqs <- test_.fetchSequence()
   motifs <- TReNA:::.getScoredMotifs(seqs, min.match.percentage=95)
   checkEquals(unlist(lapply(motifs, nrow)), c(2, 4, 3))

} # test_.getScoredMotifs
#----------------------------------------------------------------------------------------------------
test_mef2cPromoter.incrementally <- function()
{
   printf("--- test_mef2cPromoter.incrementally")

    # chr5:88,813,245-88,832,344
   chrom <- "chr5"
   start <- 88824500
   end   <- 88832344

   tbl.regions <- TReNA:::.getRegions(chrom, start, end, score.threshold=0)   # 31 open chromatin regions
   checkEquals(dim(tbl.regions), c(18, 5))
   checkTrue(all(tbl.regions$chrom == chrom))
   checkTrue(all(tbl.regions$chromStart >= start))
   checkTrue(all(tbl.regions$chromEnd <= end))

   tbl.regions <- TReNA:::.getRegions(chrom, start, end, score.threshold=200)   # 31 open chromatin regions
   checkEquals(dim(tbl.regions), c(13, 5))

   tbl.regions <- TReNA:::.getRegions(chrom, start, end, score.threshold=700)   # 31 open chromatin regions
   checkEquals(dim(tbl.regions), c(2, 5))
   checkTrue(all(tbl.regions$score >= 700))

   seqs <- TReNA:::.getSequence(tbl.regions)
   motifs <- TReNA:::.getScoredMotifs(seqs, min.match.percentage=95)

} # test_mef2cPromoter.incrementally
#----------------------------------------------------------------------------------------------------
test_mef2cPromoter.normalUse <- function()
{
   printf("--- test_mef2cPromoter.normalUse")

   load(system.file(package="TReNA", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
   checkTrue(exists("mtx.sub"))
   hdcf <- HumanDNAseClusterFilter(mtx.sub)


    # chr5:88,813,245-88,832,344: has just a few high scoring clusters
   chrom <- "chr5"
   start <- 88824500
   end   <- 88832344

   args <- list(chrom=chrom, start=start, end=end,
                region.score.threshold=700,
                motif.min.match.percentage=95)

   x <- getCandidates(hdcf, args)
   checkEquals(sort(names(x)), c("tbl", "tfs"))
   checkTrue(all(c("chrom", "regionStart", "regionEnd", "regionScore", "sourceCount", "motif", "match", "motif.start", "motif.end", "motif.width", "motif.score", "strand", "tf") %in% colnames(x$tbl)))
   checkTrue(nrow(x$tbl) > 15)     # 16 on (29 mar 2017)
   checkTrue(length(x$tfs) > 200)  # 217 on (29 mar 2017)
   checkEquals(length(which(duplicated(x$tfs))), 0)

} # test_mef2cPromoter.normalUse
#----------------------------------------------------------------------------------------------------
