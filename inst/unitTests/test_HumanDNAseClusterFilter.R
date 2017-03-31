library(TReNA)
library(MotifDb)
library(RUnit)
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_defaultConstructor()
   test_.getRegions()
   test_.fetchSequence()
   test_.matchForwardAndReverse()
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

   tbl.regions <- TReNA:::.getRegions(chrom, start, end, score.threshold=200, quiet=FALSE)
   checkTrue(nrow(tbl.regions) > 20)
   checkEquals(colnames(tbl.regions), c("chrom", "chromStart", "chromEnd", "score", "sourceCount"))
   checkTrue(all(tbl.regions$chrom == chrom))
   checkTrue(all(tbl.regions$chromStart >= start))
   checkTrue(all(tbl.regions$chromStart <= end))
   checkTrue(all(tbl.regions$chromEnd >= start))
   checkTrue(all(tbl.regions$chromEnd <= end))

      # some DHS regions, good for testing small region overlap handling
      #    chrom chromStart chromEnd score sourceCount
      # 1   chr1   88801880 88802150   573          64
      # 2   chr1   88802480 88802930   287          13
      # 3   chr1   88803000 88803270   541          60
      # 4   chr1   88811140 88811290   100           1
      # 5   chr1   88811400 88811550    68           1

     # a small region, entirely within a DHS region
   chrom = "chr1"
   start <-  88802520
   end   <-  88802530
   tbl.regions <- TReNA:::.getRegions(chrom, start, end, quiet=FALSE)
   checkEquals(nrow(tbl.regions), 1)
   checkEquals(tbl.regions$chromStart, start)
   checkEquals(tbl.regions$chromEnd, end)

     # a small region, overhanging that DHS region at its upperbound, on the right
   chrom = "chr1"
   start <-  88802145
   end   <-  88802155
   tbl.regions <- TReNA:::.getRegions(chrom, start, end, quiet=FALSE)
   checkEquals(nrow(tbl.regions), 1)
   checkEquals(tbl.regions$chromStart, start)
   checkEquals(tbl.regions$chromEnd, 88802150)

     # a small region, overhanging that DHS region at its lowerbound, on the left
   chrom = "chr1"
   start <-  88801878
   end   <-  88801883
   tbl.regions <- TReNA:::.getRegions(chrom, start, end, quiet=FALSE)
   checkEquals(nrow(tbl.regions), 1)
   checkEquals(tbl.regions$chromStart, 88801880)
   checkEquals(tbl.regions$chromEnd, end)

      # another completely contained-in-DHS-region, small area of interest
   chrom <- "chr1"
   start <- 167830160
   end   <- 167830180
   tbl.regions <- TReNA:::.getRegions(chrom, start, end, quiet=FALSE)
   checkEquals(nrow(tbl.regions), 1)
   checkEquals(tbl.regions$chromStart, start)
   checkEquals(tbl.regions$chromEnd, end)

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
test_.matchForwardAndReverse <- function()
{
   printf("--- test_.matchForwardAndReverse")

   chrom <- "chr1"
   start <- 167829960
   end   <- 167830230

   motifName <- "MA0476.1"
   mtx <- query(MotifDb, motifName)[[1]];
   sequence <- as.character(getSeq(reference.genome, chrom, start, end))

   tbl <- TReNA:::.matchForwardAndReverse(sequence, mtx, motifName, min.match.percentage=90, quiet=TRUE)

   checkEquals(nrow(tbl), 1)
   checkEquals(tbl$start, 57)
   checkEquals(tbl$end, 67)
   checkEquals(tbl$width, 11)
   checkEqualsNumeric(tbl$score, 7.98, tol=0.1)
   checkEquals(tbl$motif, "MA0476.1")
   checkEquals(tbl$strand, "+")
   checkEquals(tbl$match, substring(sequence, tbl$start, tbl$end))

   motifName <- "MA0478.1"
   mtx <- query(MotifDb, motifName)[[1]];

   tbl <- TReNA:::.matchForwardAndReverse(sequence, mtx, motifName, min.match.percentage=90, quiet=FALSE)
      # fimo finds:
      #  X.pattern.name sequence.name start stop strand   score  p.value q.value matched.sequence
      #        MA0478.1        ma0478    58   68      - 14.5455 1.21e-05   0.006      GCATGACTCAG
      #   library(FimoClient)
      #   FIMO_HOST <- "whovian"
      #   FIMO_PORT <- 5558
      #   fc <<- FimoClient(FIMO_HOST, FIMO_PORT, quiet=TRUE)
      #   tbl.fc <- requestMatch(fc, list(ma0478=sequence))
      #   subset(tbl.fc, X.pattern.name=="MA0478.1")
      #
      #   X.pattern.name sequence.name start stop strand   score  p.value q.value matched.sequence
      #         MA0478.1        ma0478    58   68      - 14.5455 1.21e-05   0.006      GCATGACTCAG

   checkEquals(nrow(tbl), 1)
   tbl$score <- round(tbl$score, 2)  # for easy comparison
   checkEquals(as.list(tbl), list(start=58, end=68, width=11, score=7.94, motif="MA0478.1",
                                  match="GCATGACTCAG", strand="-"))

} # test_.matchForwardAndReverse
#----------------------------------------------------------------------------------------------------
checkFimo <- function()
{

} # checkFimo
#----------------------------------------------------------------------------------------------------
test_.findMofits <- function()
{
   printf("--- test_.findMotifs")
   x <- .findMotifs("ACTATTCCCCT", pfms, 90)
   seqs <- test_.fetchSequence()
   motifs <- TReNA:::.getScoredMotifs(seqs, min.match.percentage=95)
   checkEquals(unlist(lapply(motifs, nrow)), c(2, 4, 3))

} # test_.findMotifs
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
# one of the igap gwas alzheimer's snps, hand-verified to fall within a motif in a dhs region
test_rs34423320 <- function()
{
   printf("--- test_rs34423320")
     # chr1 167830170 167830170  rs34423320-C-T

   load(system.file(package="TReNA", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
   checkTrue(exists("mtx.sub"))
   hdcf <- HumanDNAseClusterFilter(mtx.sub)

     # chr1 167830170 167830170  rs34423320-C-T

   chrom <- "chr1"
   snp.loc.hg38 <- 167830170
   checkEquals(as.character(getSeq(reference.genome, chrom, snp.loc.hg38, snp.loc.hg38)), "C")

   start <- snp.loc.hg38 - 10
   end   <- snp.loc.hg38 + 10

   args <- list(chrom=chrom, start=start, end=end,
                region.score.threshold=200,
                motif.min.match.percentage=90)

   x <- getCandidates(hdcf, args)
   checkEquals(sort(names(x)), c("tbl", "tfs"))
   checkTrue(all(c("chrom", "regionStart", "regionEnd", "regionScore", "sourceCount", "motif", "match", "motif.start", "motif.end", "motif.width", "motif.score", "strand", "tf") %in% colnames(x$tbl)))
   checkTrue(nrow(x$tbl) == 2)
   checkTrue(length(x$tfs) == 25
   checkEquals(length(which(duplicated(x$tfs))), 0)
   checkEquals(x$tbl$motif, c("MA0081.1", "MA0056.1"))

   seq.wt <-  as.character(getSeq(reference.genome, "chr1", 167830170-10, 167830170+10))
   seq.mut <- sprintf("%s%s%s", substr(seq.wt, 1, 10), "T", substr(seq.wt, 12, 21))

   TReNA:::.findMotifs(seq.wt, hdcf@pfms["MA0081.1"], 90)
   TReNA:::.findMotifs(seq.wt, hdcf@pfms["MA0056.1"], 90)
   TReNA:::.findMotifs(seq.mut, hdcf@pfms["MA0081.1"], 90)
   TReNA:::.findMotifs(seq.mut, hdcf@pfms["MA0056.1"], 90)

   TReNA:::.getScoredMotifs(list(seq.wt, seq.mut), min.match.percentage=90, quiet=TRUE)

} # test_rs34423320
#----------------------------------------------------------------------------------------------------
