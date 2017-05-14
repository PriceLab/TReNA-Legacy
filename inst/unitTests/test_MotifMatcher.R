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
   test_parseVariantString()
   test_.injectSnp()
   test_getSequenceWithVariants()
   test_.getScoredMotifs()
   test_findMatchesByChromosomalRegion()
   test_findMatchesByChromosomalRegion.twoAlternateAlleles()

} # runTests
#------------------------------------------------------------------------------------------------------------------------
test_basicConstructor <- function(reuse=FALSE)
{
   if(!reuse)  printf("--- test_basicConstructor")

   mm <- MotifMatcher()
   mm <- MotifMatcher(name="rs13384219.neighborhood", genomeName="hg38")

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

   motifs <- TReNA:::.getScoredMotifs(seqs, min.match.percentage=80)  # relatexed threshold
   checkEquals(unlist(lapply(motifs, nrow)), c(2, 0, 20))

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

   for(i in 1:3){
     tbl.regions <- data.frame(chrom=chroms[1], start=starts[1], end=ends[1], stringsAsFactors=FALSE)
     tbl.regions <- getSequence(mm, tbl.regions)
     checkEquals(colnames(tbl.regions), c("chrom", "start", "end", "seq", "alt"))
     seqs <- tbl.regions$seq
     checkEquals(length(seqs), 1)
     checkEquals(unlist(lapply(seqs, nchar)), 1 + tbl.regions$end - tbl.regions$start)
     }

        # --- how to get multiple sequences
    tbl.regions <- data.frame(chrom=chroms, start=starts, end=ends, stringsAsFactors=FALSE)
    x <- lapply(1:nrow(tbl.regions), function(r) getSequence(mm, tbl.regions[r,]))
    tbl.regions <- do.call(rbind, x)
    invisible(tbl.regions$seq)

} # test_getSequence
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
test_.injectSnp <- function()
{
   printf("--- test_.injectSnp")

   tbl.regions <- data.frame(chrom="chr2", start=57907318, end=57907328, seq="CATGCAAATTA", stringsAsFactors=FALSE)
   checkEquals(dim(tbl.regions), c(1, 4))
   tbl.variants <- data.frame(chrom="chr2", loc=57907323, wt="A", mut="G", stringsAsFactors=FALSE)
   tbl.regions.new <- TReNA:::.injectSnp(tbl.regions, tbl.variants)
   checkEquals(dim(tbl.regions.new), c(1, 5))
   checkEquals(tbl.regions.new$seq, "CATGCGAATTA")
   checkEquals(tbl.regions.new$alt, "chr2:57907323(A->G)")


     #---- a snp with two alternate alleles: rs3763040
   tbl.regions <- data.frame(chrom="chr18", start=26864405, end=26864415, seq="ACAGGGGGGTG",stringsAsFactors=FALSE)
   tbl.variants <- data.frame(chrom=rep("chr18",2), loc=rep(26864410,2), wt=rep("G",2) , mut=c("A", "T"),stringsAsFactors=FALSE)
   tbl.regions.new <- TReNA:::.injectSnp(tbl.regions, tbl.variants)
   checkEquals(dim(tbl.regions.new), c(2, 5))
     #                                      |
   checkEquals(tbl.regions.new$seq, c("ACAGGAGGGTG",
                                      "ACAGGTGGGTG"))
   checkEquals(tbl.regions.new$alt, c("chr18:26864410(G->A)",
                                      "chr18:26864410(G->T)"))

      # that same snp but with a region that does not containt i
   tbl.regions <- data.frame(chrom="chr2", start=26864405, end=26864415, seq="ACAGGGGGGTG",stringsAsFactors=FALSE)
   tbl.variants <- data.frame(chrom=rep("chr18",2), loc=rep(26864410,2), wt=rep("G",2) , mut=c("A", "T"),stringsAsFactors=FALSE)
   tbl.regions.new <- TReNA:::.injectSnp(tbl.regions, tbl.variants)
   checkEquals(tbl.regions, tbl.regions.new)

      # now with one overlapping region, two non-overlapping
   tbl.regions <- data.frame(chrom="chr18",
                             start=c(26864305, 26864405, 26864505),
                             end=  c(26864315, 26864415, 26864515),
                             seq=c("CTAGCCCTTAG", "ACAGGGGGGTG","CTCTAGAGGAA"),
                             stringsAsFactors=FALSE)
   tbl.variants <- data.frame(chrom=rep("chr18",2), loc=rep(26864410,2), wt=rep("G",2) , mut=c("A", "T"),stringsAsFactors=FALSE)
   tbl.regions.new <- TReNA:::.injectSnp(tbl.regions, tbl.variants)
   checkEquals(dim(tbl.regions.new), c(4, 5))
   checkEquals(tbl.regions.new$chrom, rep("chr18", 4))
   checkEquals(tbl.regions.new$start, c(26864305, 26864405, 26864405, 26864505))
   checkEquals(tbl.regions.new$end, c(26864315, 26864415, 26864415, 26864515))
   checkEquals(tbl.regions.new$seq, c("CTAGCCCTTAG", "ACAGGAGGGTG", "ACAGGTGGGTG", "CTCTAGAGGAA"))
   checkEquals(tbl.regions.new$alt, c("wt", "chr18:26864410(G->A)", "chr18:26864410(G->T)", "wt"))

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

   checkEquals(tbl.wt$alt, "wt")
   checkEquals(tbl.mut$alt, "chr2:57907323(A->G)")
   checkEquals(tbl.wt$seq,  "CATGCAAATTA")
   checkEquals(tbl.mut$seq, "CATGCGAATTA")

    # now a snp with two variants
   chromosome <- "chr18"
   loc <- 26864410
   tbl.regions <- data.frame(chrom=chromosome, start=loc-5, end=loc+5, stringsAsFactors=FALSE)
   rsid <- "rs3763040"

   tbl.wt <- getSequence(mm, tbl.regions)
   checkEquals(dim(tbl.wt), c(1, 5))
   checkEquals(tbl.wt$seq, "ACAGGGGGGTG")
   checkEquals(tbl.wt$alt, "wt")

   tbl.mut <- getSequence(mm, tbl.regions, rsid)
   checkEquals(dim(tbl.mut), c(2, 5))
   checkEquals(tbl.mut$seq, c("ACAGGAGGGTG", "ACAGGTGGGTG"))
   checkEquals(tbl.mut$alt, c("chr18:26864410(G->A)", "chr18:26864410(G->T)"))

} # test_getSequenceWithVariants
#----------------------------------------------------------------------------------------------------
test_findMatchesByChromosomalRegion <- function()
{
   printf("--- test_findMatchesByChromosomalRegion")
     # the vrk2 promoter snp,  chr2:57907313-57907333

   motifMatcher <- MotifMatcher(name="rs13384219.neighborhood", genomeName="hg38")

   tbl.regions <- data.frame(chrom="chr2", start=57907313, end=57907333, stringsAsFactors=FALSE)
   x <- findMatchesByChromosomalRegion(motifMatcher, tbl.regions, pwmMatchMinimumAsPercentage=90)
   checkEquals(sort(names(x)), c("tbl", "tfs"))
   checkEquals(dim(x$tbl), c(20, 13))
   checkEquals(unique(x$tbl$alt), "wt")
   checkTrue(length(x$tfs) > 150)

} # test_findMatchesByChromosomalRegion
#----------------------------------------------------------------------------------------------------
test_findMatchesByMultipleChromosomalRegions <- function()
{
   printf("--- test_findMatchesByMultipleChromosomalRegions")
     # the vrk2 promoter snp,  chr2:57907313-57907333
   motifMatcher <- MotifMatcher(name="rs13384219.neighborhood", genomeName="hg38")

   tbl.regions <- data.frame(chrom=c("chr2", "chr18"),
                             start=c(57907313, 26864400),
                             end=c(57907333, 26864420),
                             stringsAsFactors=FALSE)

   x <- findMatchesByChromosomalRegion(motifMatcher, tbl.regions, pwmMatchMinimumAsPercentage=90)
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
   checkEquals(sort(unique(x0$tbl$motifname)),
       c("MA0004.1", "MA0104.3", "MA0522.2", "MA0622.1", "MA0626.1", "MA0745.1", "MA0820.1", "MA0824.1","MA0830.1"))
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
