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
   test_parseVariantString()
   test_getSequenceWithVariants()
   test_.matchPwmForwardAndReverse()
   test_.getScoredMotifs()
   test_findMatchesByChromosomalRegion()
   test_findMatchesByChromosomalRegion.twoAlternateAlleles()

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
   ends   <- c(88819710, 88820705, 88821000)
   tbl.regions <- data.frame(chrom=chroms, start=starts, end=ends, stringsAsFactors=FALSE)

   tbl.regions <- getSequence(mm, tbl.regions)
   checkEquals(colnames(tbl.regions), c("chrom", "start", "end", "seq", "variant"))
   seqs <- tbl.regions$seq
   checkEquals(length(seqs), 3)
   checkEquals(unlist(lapply(seqs, nchar)), 1 + tbl.regions$end - tbl.regions$start)

   invisible(seqs)

} # test_getSequence
#----------------------------------------------------------------------------------------------------
test_.injectSnp <- function()
{
   printf("--- test_.injectSnp")
   tbl.regions <- data.frame(chrom="chr2", start=57907318, end=57907328, seq="CATGCAAATTA",
                             stringsAsFactors=FALSE)
   variant.info <- list(chrom="chr2", loc=57907323, wt="A", mut="G")
   tbl.overlaps <- data.frame(queryHits=1,  subjectHits=1, stringsAsFactors=FALSE)
   tbl.regionsInjected <- TReNA:::.injectSnp(tbl.regions, queryIndex=1, subjectIndex=1,
                                             variant.info)
   checkEquals(tbl.regionsInjected$seq[1], "CATGCGAATTA")

     #---- a snp with two alternate alleles
   tbl.regions <- data.frame(chrom="chr18", start=26864390, end=26864430,
                             seq="GGAAGAGCTGGCTCCACAGGGGGGTGGCCAGCCACATCCCA",
                             stringsAsFactors=FALSE)
   variant.info <- list(chrom="chr18", loc=26864410, wt="G", mut=c("A", "T"))
   tbl.regionsInjected <- TReNA:::.injectSnp(tbl.regions, queryIndex=1, subjectIndex=1,
                                             variant.info)

   checkEquals(tbl.regionsInjected$seq,
               c("GGAAGAGCTGGCTCCACAGGAGGGTGGCCAGCCACATCCCA",
                 "GGAAGAGCTGGCTCCACAGGTGGGTGGCCAGCCACATCCCA"))

 } # test_.injectSnp
#----------------------------------------------------------------------------------------------------
#  rs13384219  A->G
#  gtcagtagtggtggaaccagcatgc[A/G]aattagacaatgtgacttcatagcc
#  Chromosome: 2:57907323
#  chr2:57907323:A:G
test_getSequenceWithVariants <- function()
{
   printf("--- test_getSequenceWithVariants")

   mm <- MotifMatcher(genomeName="hg38")

   chroms <- "chr2"
   starts <- 57907318
   ends   <- 57907328
   tbl.regions <- data.frame(chrom=chroms, start=starts, end=ends, stringsAsFactors=FALSE)

   tbl.wt <- getSequence(mm, tbl.regions)
   tbl.mut <- getSequence(mm, tbl.regions, "rs13384219")

   checkEquals(tbl.wt$variant, "")
   checkEquals(tbl.mut$variant, "chr2:57907323:G")
   checkEquals(tbl.wt$seq,  "CATGCAAATTA")
   checkEquals(tbl.mut$seq, "CATGCGAATTA")

    # now a snp with two variants
   chromosome <- "chr18"
   loc <- 26864410
   tbl.regions <- data.frame(chrom=chromosome, start=loc-10, end=loc+10, stringsAsFactors=FALSE)
   rsid <- "rs3763040"
   tbl.mut <- getSequence(mm, tbl.regions, rsid)

   mm <- MotifMatcher(name=", two alternate alleles", genomeName="hg38")

} # test_getSequenceWithVariants
#----------------------------------------------------------------------------------------------------
test_findMatchesByChromosomalRegion <- function()
{
   printf("--- test_findMatchesByChromosomalRegion")
     # the vrk2 promoter snp,  chr2:57907313-57907333
   mm <- MotifMatcher(name="rs13384219.neighborhood", genomeName="hg38")

   tbl.regions <- data.frame(chrom="chr2", start=57907313, end=57907333, stringsAsFactors=FALSE)
   x <- findMatchesByChromosomalRegion(mm, tbl.regions, pwmMatchMinimumAsPercentage=90)
   checkEquals(sort(names(x)), c("tbl", "tfs"))
   checkEquals(dim(x$tbl), c(20, 13))
   checkEquals(unique(x$tbl$variant), "")
   checkTrue(length(x$tfs) > 150)

} # test_findMatchesByChromosomalRegion
#----------------------------------------------------------------------------------------------------
test_findMatchesByMultipleChromosomalRegions <- function()
{
   printf("--- test_findMatchesByMultipleChromosomalRegions")
     # the vrk2 promoter snp,  chr2:57907313-57907333
   mm <- MotifMatcher(name="rs13384219.neighborhood", genomeName="hg38")

   tbl.regions <- data.frame(chrom=c("chr2", "chr18"),
                             start=c(57907313, 26864400),
                             end=c(57907333, 26864420),
                             stringsAsFactors=FALSE)

   x <- findMatchesByChromosomalRegion(mm, tbl.regions, pwmMatchMinimumAsPercentage=90)
   checkEquals(sort(names(x)), c("tbl", "tfs"))
   checkEquals(dim(x$tbl), c(22, 13))
   checkTrue(length(x$tfs) > 150)
      # two motifs from chr18, the rest from chr2
   checkEquals(nrow(subset(x$tbl, chrom=="chr18")), 2)
   checkEquals(nrow(subset(x$tbl, chrom=="chr2")), 20)

} # test_findMatchesByMultipleChromosomalRegions
#----------------------------------------------------------------------------------------------------
test_findMatchesByChromosomalRegion.twoAlternateAlleles <- function()
{
   printf("--- test_findMatchesByChromosomalRegion.twoAlternateAlleles")

     # rsid <- "rs3763040"  # 18:26864410 A/C/T     (T/G/A)
     # getSequence(mm, data.frame(chrom="chr18", start=26864410, end=26864410))[1] "G"
     #                                |
     # chr18:26864405-26864415   ACAGGGGGGTG
     # thus variants must be (wrt to + strand) T and A
     # target.gene <- "APQR"
     # genome <- "hg38"

   chromosome <- "chr18"
   loc <- 26864410

   mm <- MotifMatcher(name="rs3763040, two alternate alleles", genomeName="hg38")
   tbl.regions <- data.frame(chrom=chromosome, start=loc-6, end=loc+3, stringsAsFactors=FALSE)
   x0 <- findMatchesByChromosomalRegion(mm, tbl.regions, pwmMatchMinimumAsPercentage=80L, "chr18:26864410:G:T")
      # two low-scoring ~80% matches to the sequence with C substituted
      #              |
      #        CACAGGTGGG
   checkEquals(length(grep("CAGGTG", x0$tbl$match)), nrow(x0$tbl))
   x1 <- findMatchesByChromosomalRegion(mm, tbl.regions, pwmMatchMinimumAsPercentage=80L, "chr18:26864410:G:A")
      # 10 often higher-scoring matches to the core of the sequence with T substituted
      #              |
      #          CAGGTG
   checkEquals(length(grep("CAGGAG", x1$tbl$match)), nrow(x1$tbl))

     # now find both results from one call.  MotifMatcher looksup the rsid,
     # learns of three alleles expressed in an ambiguity code (D: AGT) then
     # subtracts out the reference (G) returning a table which combines motifs matching
     # the A and T injected motifs
   x2 <- findMatchesByChromosomalRegion(mm, tbl.regions, pwmMatchMinimumAsPercentage=80L, "rs3763040")
   checkTrue(all(x2$tbl$chrom == chromosome))
   checkEquals(sort(names(x2)), c("tbl", "tfs"))
   checkEquals(length(grep("CAGGTG", x2$tbl$match)), nrow(x0$tbl))
   checkEquals(length(grep("CAGGAG", x2$tbl$match)), nrow(x1$tbl))

} # test_findMatchesByChromosomalRegion.twoAlternateAlleles
#----------------------------------------------------------------------------------------------------
test_parseVariantString <- function()
{
   mm <- MotifMatcher(name="rs13384219.neighborhood", genomeName="hg38")

   tbl.variant <- TReNA:::parseVariantString(mm, "rs13384219")
   checkEquals(dim(tbl.variant), c(1, 4))
   checkEquals(tbl.variant$chrom, "chr2")
   checkEquals(tbl.variant$loc, 57907323)
   checkEquals(tbl.variant$wt, "A")
   checkEquals(tbl.variant$mut, "G")

       # the same snp, but here expressed as a string (not an rsid)
   tbl.v2 <- TReNA:::parseVariantString(mm, "chr2:57907323:A:G")
   checkEquals(dim(tbl.v2), c(1, 4))
   checkEquals(tbl.v2$chrom, "chr2")
   checkEquals(tbl.v2$loc, 57907323)
   checkEquals(tbl.v2$wt, "A")
   checkEquals(tbl.v2$mut, "G")

   tbl.2vars <- TReNA:::parseVariantString(mm, "rs3763040")  # snp with two variant alleles
   checkEquals(dim(tbl.2vars), c(2,4))
   checkEquals(tbl.2vars$chrom, rep("chr18", 2))
   checkEquals(tbl.2vars$loc, rep(26864410, 2))
   checkEquals(tbl.2vars$wt, rep("G", 2))
   checkEquals(tbl.2vars$mut, c("A", "T"))

} # test_parseVariantString
#----------------------------------------------------------------------------------------------------
