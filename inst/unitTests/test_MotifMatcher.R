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
   test_.parseVariantString()
   test_.injectSnp()
   test_getSequenceWithVariants()
   test_.getScoredMotifs()
   test_findMatchesByChromosomalRegion()
   test_findMatchesByChromosomalRegion_contrastReferenceWithVariant()
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
   checkEquals(colnames(tbl.1),
               c("start","end","width","score","maxScore","relativeScore","motif","match","strand"))
   tbl.2 <- TReNA:::.matchPwmForwardAndReverse(sequence, mtx.2, motifName.2, min.match.percentage=60, quiet=TRUE)
   checkEquals(nrow(tbl.2), 4)
   checkEquals(colnames(tbl.2),
               c("start","end","width","score","maxScore","relativeScore","motif","match","strand"))

} # test_.matchForwardAndReverse
#----------------------------------------------------------------------------------------------------
test_.getScoredMotifs <- function()
{
   printf("--- test_.getScoredMotifs")
   seqs <- test_getSequence(indirect=TRUE)

   motifs <- TReNA:::.getScoredMotifs(seqs, min.match.percentage=80)  # relatexed threshold
   checkEquals(unlist(lapply(motifs, nrow)), c(2, 0, 20))
   checkEquals(colnames(motifs[[1]]),
                        c("start", "end", "width", "score", "maxScore", "relativeScore", "motif", "match", "strand"))
   checkEquals(colnames(motifs[[3]]),
                        c("start", "end", "width", "score", "maxScore", "relativeScore", "motif", "match", "strand"))

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
     tbl.regions.noSeq <- data.frame(chrom=chroms[i], start=starts[i], end=ends[i], stringsAsFactors=FALSE)
     tbl.regions <- getSequence(mm, tbl.regions.noSeq)
     checkEquals(colnames(tbl.regions), c("chrom", "start", "end", "seq", "status"))
     seqs <- tbl.regions$seq
     checkEquals(length(seqs), 1)
     checkEquals(unlist(lapply(seqs, nchar)), 1 + tbl.regions$end - tbl.regions$start)
     }

        # --- how to get multiple sequences
    tbl.regions.noSeq <- data.frame(chrom=chroms, start=starts, end=ends, stringsAsFactors=FALSE)
    x <- lapply(1:nrow(tbl.regions.noSeq), function(r) getSequence(mm, tbl.regions.noSeq[r,]))
    tbl.regions <- do.call(rbind, x)
    invisible(tbl.regions$seq)

} # test_getSequence
#----------------------------------------------------------------------------------------------------
test_.parseVariantString <- function()
{
   mm <- MotifMatcher(name="rs13384219.neighborhood", genomeName="hg38")

   tbl.variant <- TReNA:::.parseVariantString(mm, "rs13384219")
   checkEquals(dim(tbl.variant), c(1, 4))
   checkEquals(tbl.variant$chrom, "chr2")
   checkEquals(tbl.variant$loc, 57907323)
   checkEquals(tbl.variant$wt, "A")
   checkEquals(tbl.variant$mut, "G")

       # the same snp, but here expressed as a string (not an rsid)
   tbl.v2 <- TReNA:::.parseVariantString(mm, "chr2:57907323:A:G")
   checkEquals(dim(tbl.v2), c(1, 4))
   checkEquals(tbl.v2$chrom, "chr2")
   checkEquals(tbl.v2$loc, 57907323)
   checkEquals(tbl.v2$wt, "A")
   checkEquals(tbl.v2$mut, "G")

   tbl.2vars <- TReNA:::.parseVariantString(mm, "rs3763040:A")  # snp with two variant alleles; one must be specified
   checkEquals(dim(tbl.2vars), c(1,4))
   checkEquals(tbl.2vars$chrom, "chr18")
   checkEquals(tbl.2vars$loc, 26864410)
   checkEquals(tbl.2vars$wt, "G")
   checkEquals(tbl.2vars$mut, "A")

   tbl.2vars <- TReNA:::.parseVariantString(mm, "rs3763040:T")  # snp with two variant alleles; one must be specified
   checkEquals(dim(tbl.2vars), c(1,4))
   checkEquals(tbl.2vars$chrom, "chr18")
   checkEquals(tbl.2vars$loc, 26864410)
   checkEquals(tbl.2vars$wt, "G")
   checkEquals(tbl.2vars$mut, "T")

   tbl.2vars <- TReNA:::.parseVariantString(mm, "rs3763040")  # snp with two variant alleles; first is chosen by default
   checkEquals(dim(tbl.2vars), c(1,4))
   checkEquals(tbl.2vars$chrom, "chr18")
   checkEquals(tbl.2vars$loc, 26864410)
   checkEquals(tbl.2vars$wt, "G")
   checkEquals(tbl.2vars$mut, "A")

      # now fail on purpose, explictly specifying an alt allele not include in the snp
   checkException(
      tbl.2vars <- TReNA:::.parseVariantString(mm, "rs3763040:C"),
      "caught!", silent=TRUE)

} # test_.parseVariantString
#----------------------------------------------------------------------------------------------------
test_.injectSnp <- function()
{
   printf("--- test_.injectSnp")

   mm <- MotifMatcher(genomeName="hg38")
   tbl.regions.noSeq <- data.frame(chrom="chr2", start=57907318, end=57907328, stringsAsFactors=FALSE)
   tbl.regions <- getSequence(mm, tbl.regions.noSeq)
   checkEquals(tbl.regions$seq, "CATGCAAATTA")

   tbl.variants <- data.frame(chrom="chr2", loc=57907323, wt="A", mut="G", stringsAsFactors=FALSE)
   tbl.regions.new <- TReNA:::.injectSnp(tbl.regions, tbl.variants)
   checkEquals(dim(tbl.regions.new), c(1, 5))
   checkEquals(tbl.regions.new$seq, "CATGCGAATTA")
   checkEquals(tbl.regions.new$status, "mut")
   #checkEquals(tbl.regions.new$alt, "chr2:57907323(A->G)")


     #---- a snp with two alternate alleles: rs3763040

   #   tbl.regions.noSeq <- data.frame(chrom="chr18", start=26864405, end=26864415,  stringsAsFactors=FALSE)
   #   tbl.regions <- getSequence(mm, tbl.regions.noSeq)
   #   checkEquals(tbl.regions$seq, "ACAGGGGGGTG")
   #   checkEquals(tbl.regions$status, "wt")
   #
   #   tbl.variants <- data.frame(chrom=rep("chr18",2), loc=rep(26864410,2), wt=rep("G",2) , mut=c("A", "T"),stringsAsFactors=FALSE)
   #   tbl.regions.new <- TReNA:::.injectSnp(tbl.regions, tbl.variants)
   #   checkEquals(dim(tbl.regions.new), c(2, 5))
   #     #                                      |
   #   checkEquals(tbl.regions.new$seq, c("ACAGGAGGGTG",
   #                                      "ACAGGTGGGTG"))
   #   checkEquals(tbl.regions.new$alt, c("chr18:26864410(G->A)",
   #                                      "chr18:26864410(G->T)"))
   #
   #      # that same snp but with a region that does not containt it
   #   tbl.regions <- data.frame(chrom="chr2", start=26864405, end=26864415, seq="ACAGGGGGGTG",stringsAsFactors=FALSE)
   #   tbl.variants <- data.frame(chrom=rep("chr18",2), loc=rep(26864410,2), wt=rep("G",2) , mut=c("A", "T"),stringsAsFactors=FALSE)
   #   tbl.regions.new <- TReNA:::.injectSnp(tbl.regions, tbl.variants)
   #   checkEquals(tbl.regions, tbl.regions.new)
   #
   #      # now with one overlapping region, two non-overlapping
   #   tbl.regions <- data.frame(chrom="chr18",
   #                             start=c(26864305, 26864405, 26864505),
   #                             end=  c(26864315, 26864415, 26864515),
   #                             seq=c("CTAGCCCTTAG", "ACAGGGGGGTG","CTCTAGAGGAA"),
   #                             stringsAsFactors=FALSE)
   #   tbl.variants <- data.frame(chrom=rep("chr18",2), loc=rep(26864410,2), wt=rep("G",2) , mut=c("A", "T"),stringsAsFactors=FALSE)
   #   tbl.regions.new <- TReNA:::.injectSnp(tbl.regions, tbl.variants)
   #   checkEquals(dim(tbl.regions.new), c(4, 5))
   #   checkEquals(tbl.regions.new$chrom, rep("chr18", 4))
   #   checkEquals(tbl.regions.new$start, c(26864305, 26864405, 26864405, 26864505))
   #   checkEquals(tbl.regions.new$end, c(26864315, 26864415, 26864415, 26864515))
   #   checkEquals(tbl.regions.new$seq, c("CTAGCCCTTAG", "ACAGGAGGGTG", "ACAGGTGGGTG", "CTCTAGAGGAA"))
   #   checkEquals(tbl.regions.new$alt, c("wt", "chr18:26864410(G->A)", "chr18:26864410(G->T)", "wt"))

      # now try two variants in one region, different locations

   mm <- MotifMatcher(genomeName="hg38")
   tbl.variants <- data.frame(chrom=c("chr18", "chr18"),
                              loc=c(26865466, 26865469),
                              wt=c("G", "T"),
                              mut=c("C", "C"),
                              stringsAsFactors=FALSE)
   tbl.regions.noSeq <- data.frame(chrom="chr18", start=26865463, end=26865472, stringsAsFactors=FALSE)
   tbl.regions <- getSequence(mm, tbl.regions.noSeq)

   tbl.regions.new <- TReNA:::.injectSnp(tbl.regions, tbl.variants)
      #                                 |  |
   checkEquals(tbl.regions$seq,     "AAAGCATCCC")
   checkEquals(tbl.regions.new$seq, "AAACCACCCC")

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
   tbl.regions.noSeq <- data.frame(chrom=chroms, start=starts, end=ends, stringsAsFactors=FALSE)

   tbl.wt <- getSequence(mm, tbl.regions.noSeq)
   tbl.mut <- getSequence(mm, tbl.regions.noSeq, "rs13384219")

   checkEquals(tbl.wt$status, "wt")
   checkEquals(tbl.mut$status, "mut") #"chr2:57907323(A->G)")
   checkEquals(tbl.wt$seq,  "CATGCAAATTA")
   checkEquals(tbl.mut$seq, "CATGCGAATTA")

    # now a snp with two variants
   chromosome <- "chr18"
   loc <- 26864410
   tbl.regions.noSeq <- data.frame(chrom=chromosome, start=loc-5, end=loc+5, stringsAsFactors=FALSE)
   rsid <- "rs3763040"

   tbl.wt <- getSequence(mm, tbl.regions.noSeq)
   checkEquals(dim(tbl.wt), c(1, 5))
   checkEquals(tbl.wt$seq, "ACAGGGGGGTG")
   checkEquals(tbl.wt$status, "wt")

   tbl.mut <- getSequence(mm, tbl.regions.noSeq, rsid)
   checkEquals(dim(tbl.mut), c(1, 5))
   checkEquals(tbl.mut$seq, "ACAGGAGGGTG")
   checkEquals(tbl.mut$status, "mut")    # c("chr18:26864410(G->A)", "chr18:26864410(G->T)"))

   tbl.mut <- getSequence(mm, tbl.regions.noSeq, "rs3763040:A")
   checkEquals(dim(tbl.mut), c(1, 5))
   checkEquals(tbl.mut$seq, "ACAGGAGGGTG")
   checkEquals(tbl.mut$status, "mut")    # c("chr18:26864410(G->A)", "chr18:26864410(G->T)"))

   tbl.mut <- getSequence(mm, tbl.regions.noSeq, "rs3763040:T")
   checkEquals(dim(tbl.mut), c(1, 5))
   checkEquals(tbl.mut$seq, "ACAGGTGGGTG")
   checkEquals(tbl.mut$status, "mut")    # c("chr18:26864410(G->A)", "chr18:26864410(G->T)"))

      # now provide two variants (at two locations), here just three bases apart
      #  rs750694782: chr18:26865466
      #  rs3875089:   chr18:26865469

   variants <- c("rs750694782", "rs3875089")
   tbl.regions.noSeq <- data.frame(chrom="chr18", start=26865463, end=26865472, stringsAsFactors=FALSE)
   tbl.wt <- getSequence(mm, tbl.regions.noSeq)
   checkEquals(as.list(tbl.wt[, c("seq", "status")]),  list(seq="AAAGCATCCC", status="wt"))
   tbl.mut <- getSequence(mm, tbl.regions.noSeq, variants) #        |  |
   checkEquals(as.list(tbl.mut[, c("seq", "status")]), list(seq="AAACCACCCC", status="mut"))

} # test_getSequenceWithVariants
#----------------------------------------------------------------------------------------------------
test_findMatchesByChromosomalRegion <- function()
{
   printf("--- test_findMatchesByChromosomalRegion")
     # the vrk2 promoter snp,  chr2:57907313-57907333

   motifMatcher <- MotifMatcher(name="rs13384219.neighborhood", genomeName="hg38", quiet=FALSE)

   tbl.regions <- data.frame(chrom="chr2", start=57907313, end=57907333, stringsAsFactors=FALSE)
   x <- findMatchesByChromosomalRegion(motifMatcher, tbl.regions, pwmMatchMinimumAsPercentage=92)
   checkEquals(sort(names(x)), c("tbl", "tfs"))
   checkEquals(dim(x$tbl), c(7, 13))
   checkEquals(unique(x$tbl$status), "wt")
   checkTrue(length(x$tfs) > 50)

      # the best match (longest, all bases in agreement with the motif:
      # http://jaspar.genereg.net/cgi-bin/jaspar_db.pl?rm=present&collection=CORE&ID=MA0792.1

   best <- as.list(x$tbl[1,])
   checkEquals(best$motifName, "MA0792.1")
   checkEquals(best$chrom, "chr2")
   checkEquals(best$motifStart, 57907318)
   checkEquals(best$motifEnd, 57907326)
   checkEquals(best$strand, "+")
   checkEqualsNumeric(best$motifScore, 6.707832)
   checkEqualsNumeric(best$motifRelativeScore, 0.9208761)
   checkEquals(best$match, "CATGCAAAT")
   checkEquals(best$chromStart, 57907313)
   checkEquals(best$chromEnd, 57907333)
   checkEquals(best$seq, "ACCAGCATGCAAATTAGACAA")
   checkEquals(best$status, "wt")
   checkEquals(best$tf, "POU5F1B")

} # test_findMatchesByChromosomalRegion
#----------------------------------------------------------------------------------------------------
test_findMatchesByChromosomalRegion_contrastReferenceWithVariant <- function()
{
   printf("--- test_findMatchesByChromosomalRegion_contrastReferenceWithVariant")

     # the vrk2 promoter snp,  chr2:57907313-57907333

   motifMatcher <- MotifMatcher(name="rs13384219.neighborhood", genomeName="hg38", quiet=TRUE)

   tbl.regions <- data.frame(chrom="chr2", start=57907313, end=57907333, stringsAsFactors=FALSE)
   x.wt <- findMatchesByChromosomalRegion(motifMatcher, tbl.regions, pwmMatchMinimumAsPercentage=92)
   x.mut <- findMatchesByChromosomalRegion(motifMatcher, tbl.regions, pwmMatchMinimumAsPercentage=92, "rs13384219")

   tbl <- rbind(x.wt$tbl[, c(1,12, 2,3,4,5,7,8,13)], x.mut$tbl[, c(1,12, 2,3,4,5,7,8,13)])
   tbl <- tbl[order(tbl$motifName, tbl$motifRelativeScore, decreasing=TRUE),]
   motifs.in.both <- tbl$motifName[which(duplicated(tbl$motifName))]
   tbl.summary <- subset(tbl, !motifName %in% motifs.in.both)
      # three wt motifs lost: MA0792.1, MA0701.1, MA0075.2
   checkTrue(all(c("MA0792.1", "MA0701.1", "MA0075.2") %in% tbl.summary$motifName))
   checkEquals(unique(tbl.summary$status), "wt")

      # now try again with a more permissive matching threshold
   x.wt <- findMatchesByChromosomalRegion(motifMatcher, tbl.regions, pwmMatchMinimumAsPercentage=80)
   x.mut <- findMatchesByChromosomalRegion(motifMatcher, tbl.regions, pwmMatchMinimumAsPercentage=80, "rs13384219")
   tbl <- rbind(x.wt$tbl[, c(1,12, 2,3,4,5,7,8,13)], x.mut$tbl[, c(1,12, 2,3,4,5,7,8,13)])
   tbl <- tbl[order(tbl$motifName, tbl$motifRelativeScore, decreasing=TRUE),]

     # the POU5F1B motif MA0792.1 drops from a 92% to an 83% match.
     # assessing the relevance of this will be up to the biologist...

   checkEqualsNumeric(subset(tbl, motifName=="MA0792.1")$motifRelativeScore, c(0.92, 0.83), tol=1e-2)
   motifs.in.both <- tbl$motifName[which(duplicated(tbl$motifName))]
   tbl.summary <- subset(tbl, !motifName %in% motifs.in.both)

   # two possible new binding motifs: MA0891.1, MA0648.1
   tbl.deNovoMut <- subset(tbl.summary, motifRelativeScore > 0.85)
   checkEquals(sort(tbl.deNovoMut$motifName), c("MA0648.1", "MA0891.1"))
   checkEquals(unique(tbl.deNovoMut$status), "mut")

} # test_findMatchesByChromosomalRegion_contrastReferenceWithVariant
#----------------------------------------------------------------------------------------------------
test_findMatchesByMultipleChromosomalRegions <- function()
{
   printf("--- test_findMatchesByMultipleChromosomalRegions")
     # the vrk2 promoter snp,  chr2:57907313-57907333
   motifMatcher <- MotifMatcher(name="rs13384219.neighborhood", genomeName="hg38")

     # 21 bases on chr2, 21 bases on chr18
   tbl.regions <- data.frame(chrom=c("chr2", "chr18"),
                             start=c(57907313, 26864400),
                             end  =c(57907333, 26864420),
                             stringsAsFactors=FALSE)

   x <- findMatchesByChromosomalRegion(motifMatcher, tbl.regions, pwmMatchMinimumAsPercentage=93)
   checkEquals(sort(names(x)), c("tbl", "tfs"))
   tbl <- x$tbl

   m1 <- as.list(tbl[1,])

   checkEquals(m1$motifName, "MA0630.1")
   checkEquals(m1$chrom, "chr2")
   checkEquals(m1$motifStart, 57907317)
   checkEquals(m1$motifEnd, 57907324)
   checkEquals(m1$strand, "-")
   checkEqualsNumeric(m1$motifScore, 5.23, tol=1e-2)
   checkEqualsNumeric(m1$motifRelativeScore, 0.93, tol=1e-2)
   checkEquals(m1$match, "CAAATTAG")
   checkEquals(m1$chromStart, 57907313)
   checkEquals(m1$chromEnd, 57907333)
   checkEquals(m1$seq, "ACCAGCATGCAAATTAGACAA")
   checkEquals(m1$alt, "wt")
   tfs <- strsplit(m1$tf, ";")[[1]]
   checkTrue("SHOX" %in% tfs)     # the one tf associated to this motif by jaspar
   checkEquals(length(tfs), 43)   # do we believe this?

   m2 <- as.list(tbl[2,])
   checkEquals(m2$motifName, "MA0130.1")
   checkEquals(m2$chrom, "chr18")
   checkEquals(m2$motifStart, 26864401)
   checkEquals(m2$motifEnd, 26864406)
   checkEquals(m2$strand, "+")
   checkEquals(m2$motifScore, 5, tol=1e-2)
   checkEquals(m2$motifRelativeScore, 0.98, tol=1e-2)
   checkEquals(m2$match, "CTCCAC")
   checkEquals(m2$chromStart, 26864400)
   checkEquals(m2$chromEnd, 26864420)
   checkEquals(m2$seq, "GCTCCACAGGGGGGTGGCCAG")
   checkEquals(m2$alt, "wt")
   checkEquals(m2$tf,  "ZNF354C")

      # now repeat with a looser match threshold
    x <- findMatchesByChromosomalRegion(motifMatcher, tbl.regions, pwmMatchMinimumAsPercentage=80)
    tbl.freq <- as.data.frame(table(x$tbl$chrom))

      # one match each found for these two regions at the stringent match threshold
   checkEquals(dim(x$tbl), c(2, 13))
   checkTrue(length(x$tfs) > 40)
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
   xwt <- findMatchesByChromosomalRegion(mm, tbl.regions, pwmMatchMinimumAsPercentage=80L)
   xmut <- findMatchesByChromosomalRegion(mm, tbl.regions, pwmMatchMinimumAsPercentage=80L, "chr18:26864410:G:T")

   checkEquals(unique(xwt$tbl$seq), "CACAGGGGGG")
   checkEquals(unique(xmut$tbl$seq), "CACAGGTGGG")

      # two motifs are shared
   motifs.shared <- sort(intersect(xwt$tbl$motifName, xmut$tbl$motifName))
   checkEquals(motifs.shared, c("MA0522.2", "MA0830.1"))

      # the snp increases the scores of these two
   checkEqualsNumeric(subset(xwt$tbl,  motifName %in% motifs.shared[1])$motifScore, 5.79, tol=1e-2)
   checkEqualsNumeric(subset(xmut$tbl, motifName %in% motifs.shared[1])$motifScore, 6.79, tol=1e-2)

       # MA0522.2, TCF3, class: bHLH, family: E2A-related factors, reverse complement, from JASPAR
       #                                                          relative score   motif score
       # wt:    C   A     C    A    G    G    G    G    G    G            0.8016    5.79
       # mut:   C   A     C    A    G    G    T    G    G    G            0.9394    6.79
       # A  [2247 2554    5 7799    0  317   19    0  213 1839 ]
       # C  [2093 1978 7799    7  465 1251    0    0 2842 1014 ]
       # G  [1393 3014   21   35 7799 7799    5 7799 1321 2123 ]
       # T  [2066  253   27   43    0    0 7799    0 3423 2823 ]

   checkEqualsNumeric(subset(xwt$tbl,  motifName %in% motifs.shared[2])$motifScore, 5.89, tol=1e-2)
   checkEqualsNumeric(subset(xmut$tbl, motifName %in% motifs.shared[2])$motifScore, 6.89, tol=1e-2)


       # MA0830.1, TCF4, class: bHLH, family: E2A-related factors, reverse complement, from JASPAR
       #                                                          relative score   motif score
       # wt:     C     A     C     A     G     G     G     G     G     G            0.8047    5.89
       # mut:    C     A     C     A     G     G     T     G     G     G            0.9412    6.89
       # A  [ 5380  6301    25 20335     0  1334     0    13   409  5185 ]
       # C  [ 4838  6324 20335     0   387  2315     0     0  9640  2128 ]
       # G  [ 4812  7336     5     0 20335 20335     0 20335  2465  8284 ]
       # T  [ 5305   374    27   117    56     0 20335     6  7822  4738 ]

     # motifs lost and gained
   motifs.lost   <- sort(setdiff(xwt$tbl$motifName, xmut$tbl$motifName))
   motifs.gained <- sort(setdiff(xmut$tbl$motifName, xwt$tbl$motifName))
   checkEquals(motifs.lost, "MA0056.1")

       # MA0056.1, MZF1, class: C2H@ zinc finger factors, family: more than 3 adjacent zinc finger factors, from JASPAR
       #                                                                    relative score   motif score
       # wt:  G  G  G  G  G  G            0.8039    4.10
       # A  [ 3  0  2  0  0 18 ]
       # C  [ 5  0  0  0  0  0 ]
       # G  [ 4 19 18 19 20  2 ]
       # T  [ 8  1  0  1  0  0 ]

   checkEquals(motifs.gained, c("MA0004.1", "MA0104.3", "MA0622.1", "MA0626.1", "MA0745.1", "MA0820.1", "MA0824.1"))

       # look at (for now( just one of the gain motifs
       # MA0745.1, SNAi2, class C2H2 zinc finger factors, more than 3 adjacent zinc finger factors
       #                                                                     relative score   motif score
       # mut:     C      A      C      A      G      G      T      G      G       0.9516795      6.199955
       # A  [ 37982  56979   5150 110706  16150    361     57  16003   9858 ]
       # C  [ 15746   4756 110706   1908   1458      0    525   5615  37564 ]
       # G  [ 35246  53727   2678    784 110706 110706      0 110706  25335 ]
       # T  [ 21733  16594   4022   3454   3864    323 110706  21327  37949 ]


   checkEquals(sort(unique(xmut$tbl$motifName)),
       c("MA0004.1", "MA0104.3", "MA0522.2", "MA0622.1", "MA0626.1", "MA0745.1", "MA0820.1", "MA0824.1","MA0830.1"))
      # two low-scoring ~80% matches to the sequence with C substituted
      #              |
      #        CACAGGTGGG
   checkEquals(length(grep("CAGGTG", xmut$tbl$match)), nrow(xmut$tbl))

      # now substitute, not G->T, but G->A
      # this results in the loss of the GGGGGG motif from the wildtype, but the consequences
      # for the other two wt motifs are undetectable in the scores: both wt and this mutant
      # provide a poor match to the T expected at that position
   x1 <- findMatchesByChromosomalRegion(mm, tbl.regions, pwmMatchMinimumAsPercentage=80L, "chr18:26864410:G:A")
   checkEqualsNumeric(x1$tbl$motifRelativeScore[1], 0.805, tol=10e-2)
   checkEqualsNumeric(x1$tbl$motifRelativeScore[2], 0.802, tol=10e-2)
      # 10 often higher-scoring matches to the core of the sequence with T substituted
      #              |
      #          CAGGTG
   checkEquals(length(grep("CAGGAG", x1$tbl$match)), nrow(x1$tbl))


} # test_findMatchesByChromosomalRegion.twoAlternateAlleles
#----------------------------------------------------------------------------------------------------
