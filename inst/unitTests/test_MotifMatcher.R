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
   test_.matchPwmForwardAndReverse()
   #test_.getScoredMotifs()

} # runTests
#------------------------------------------------------------------------------------------------------------------------
test_basicConstructor <- function(reuse=FALSE)
{
   if(!reuse)  printf("--- test_basicConstructor")

   mm <- MotifMatcher()
   mm <- MotifMatcher(name="rs13384219.neighborhood",
                      sequence=sequence,
                      pwmMatchMinimumAsPercentage=85)

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
test_.findMofits <- function()
{
   printf("--- test_.findMotifs")
   x <- .findMotifs("ACTATTCCCCT", pfms, 90)
   seqs <- test_getSequence()
   motifs <- TReNA:::.getScoredMotifs(seqs, min.match.percentage=95)
   checkEquals(unlist(lapply(motifs, nrow)), c(2, 4, 3))

} # test_.findMotifs
#----------------------------------------------------------------------------------------------------
test_.getScoredMotifs <- function()
{
   printf("--- test_.getScoredMotifs")
   seqs <- test_getSequence()
   motifs <- TReNA:::.getScoredMotifs(seqs, min.match.percentage=95)
   checkEquals(unlist(lapply(motifs, nrow)), c(2, 4, 3))

} # test_.getScoredMotifs
#----------------------------------------------------------------------------------------------------
