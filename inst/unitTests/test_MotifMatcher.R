library(TReNA)
library(MotifDb)
library(RUnit)
#----------------------------------------------------------------------------------------------------
# the vrk2 promoter snp
# chr2:57907313-57907333
sequence <- "ACCAGCATGCAAATTAGACAA"
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_basicConstructor()
   test_getSequence()
   test_.matchPwmForwardAndReverse()
   test_.getScoredMotifs()
   test_findMatchesByChromosomalRegion()

} # runTests
#------------------------------------------------------------------------------------------------------------------------
test_basicConstructor <- function(reuse=FALSE)
{
   if(!reuse)  printf("--- test_basicConstructor")

   mm <- MotifMatcher()
   mm <- MotifMatcher(name="rs13384219.neighborhood",
                      genomeName="hg38")

} # test_basicConstructor
#----------------------------------------------------------------------------------------------------
test_.matchPwmForwardAndReverse <- function()
{
   printf("--- test_.matchPwmForwardAndReverse")

   motifName.1 <- "MA0507.1"
   mtx.1 <- query(MotifDb, motifName.1)[[1]];
   motifName.2 <- "MA0063.1" #"MA0627.1"
   mtx.2 <- query(MotifDb, motifName.2)[[1]];

   sequence <- "TTGTCTAATTTGCATGCTGGT"


   tbl.1 <- TReNA:::.matchPwmForwardAndReverse(sequence, mtx.1, motifName.1, min.match.percentage=90, quiet=TRUE)
   checkEquals(nrow(tbl.1), 1)
   tbl.2 <- TReNA:::.matchPwmForwardAndReverse(sequence, mtx.2, motifName.2, min.match.percentage=60, quiet=TRUE)
   checkEquals(nrow(tbl.2), 4)


} # test_.matchForwardAndReverse
#----------------------------------------------------------------------------------------------------
test_.getScoredMotifs <- function()
{
   printf("--- test_.getScoredMotifs")
   seqs <- test_getSequence(indirect=TRUE)

   motifs <- TReNA:::.getScoredMotifs(seqs, min.match.percentage=95)  # very high threshold
   checkEquals(unlist(lapply(motifs, nrow)), c(2, 4, 3))

   motifs <- TReNA:::.getScoredMotifs(seqs, min.match.percentage=80)  # relatexed threshold
   checkEquals(unlist(lapply(motifs, nrow)), c(220, 290, 129))

   motifs <- TReNA:::.getScoredMotifs(seqs, min.match.percentage=100)  # nigh impossible threshold
   checkEquals(unlist(lapply(motifs, nrow)), c(0, 0, 0))

} # test_.getScoredMotifs
#----------------------------------------------------------------------------------------------------
test_getSequence <- function(indirect=FALSE)
{
   if(!indirect)
      printf("--- test_getSequence")

   mm <- MotifMatcher(genomeName="hg38")
   chroms <- rep("chr5", 3)
   starts <- c(88819700, 88820700, 88820980)
   ends   <- c(88819910, 88820850, 88821130)
   tbl.regions <- data.frame(chrom=chroms, start=starts, end=ends, stringsAsFactors=FALSE)

   seqs <- getSequence(mm, tbl.regions)
   checkEquals(length(seqs), 3)
   checkEquals(unlist(lapply(seqs, nchar)), 1 + tbl.regions$end - tbl.regions$start)

   invisible(seqs)

} # test_getSequence
#----------------------------------------------------------------------------------------------------
test_findMatchesByChromosomalRegion <- function()
{
   printf("--- test_findMatchesByChromosomalRegion")
     # the vrk2 promoter snp,  chr2:57907313-57907333
   mm <- MotifMatcher(name="rs13384219.neighborhood",
                      genomeName="hg38")

   tbl.regions <- data.frame(chrom="chr2", start=57907313, end=57907333, stringsAsFactors=FALSE)
   x <- findMatchesByChromosomalRegion(mm, tbl.regions, pwmMatchMinimumAsPercentage=90)
   checkEquals(sort(names(x)), c("tbl", "tfs"))
   checkEquals(dim(x$tbl), c(20, 11))
   checkTrue(length(x$tfs) > 150)

} # test_findMatchesByChromosomalRegion
#----------------------------------------------------------------------------------------------------
