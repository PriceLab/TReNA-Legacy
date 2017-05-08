library(TReNA)
library(MotifDb)
library(RUnit)
#----------------------------------------------------------------------------------------------------
mef2c.tss <- 88904257
mef2c.promoter.region <- list(chrom="chr5", start=mef2c.tss-100, end=mef2c.tss+100)
mef2c.promoter.string <- with(mef2c.promoter.region, sprintf("%s:%d-%d", chrom, start, end))
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_create.vrk2.candidateFilterSpec()
   test_basicConstructor()
   test_geneSymbolToTSS()
   test_getEncodeRegulatoryTableNames()
   test_checkAllEncodeTables(quiet=FALSE)

   test_getRegulatoryRegions()
   test_getCandidates.vrk2()
   test_getCandidates.vrk2.rs13384219.neighborhood()
   #test_getSequence()
   #test_.matchForwardAndReverse()
   #test_.getScoredMotifs()

   #test_getRegulatoryRegions_hardCase()
   #test_.matchForwardAndReverse()
   #test_.getScoredMotifs()
   #test_mef2cPromoter.incrementally()
   #test_mef2cPromoter.normalUse()

} # runTests
#----------------------------------------------------------------------------------------------------
# reuse this in several tests
create.vrk2.candidateFilterSpec <- function(geneCentered=TRUE, promoter.length=1000)
{
   target.gene <- "VRK2"
   genome <- "hg38"
   chromosome <- "chr2"
   tss <- 57907651
   promoter.length <- 1000

   candidateFilterSpec <- list(filterType="EncodeDNaseClusters",
                               genomeName=genome,
                               encodeTableName="wgEncodeRegDnaseClustered",
                               pwmMatchPercentageThreshold=85L,
                               geneInfoDB="postgres://whovian/gtf",
                               regionsSpec=NA_character_,  # no explicit regions in this recipe
                               geneCenteredSpec=list(targetGene=target.gene,
                                                     tssUpstream=promoter.length,
                                                     tssDownstream=promoter.length))
   if(!geneCentered) {
      candidateFilterSpec$geneCenteredSpec <- list()
      candidateFilterSpec$regionsSpec=sprintf("%s:%d-%d", chromosome, tss-promoter.length, tss+promoter.length)
      }

   return(candidateFilterSpec)

} # create.vrk2.candidateFilterSpec
#----------------------------------------------------------------------------------------------------
# reuse this in several tests
create.vrk2.candidateFilterSpec.twoRegions <- function()
{
   target.gene <- "VRK2"
   genome <- "hg38"
   chromosome <- "chr2"
   tss <- 57907651
   promoter.length <- 1000

   cfSpec <- list(filterType="EncodeDNaseClusters",
                  genomeName=genome,
                  encodeTableName="wgEncodeRegDnaseClustered",
                  pwmMatchPercentageThreshold=85L,
                  geneInfoDB="postgres://whovian/gtf",
                        # 2 dhs regions found by inspection
                  regionsSpec=c("chr2:57906700-57906870",
                                "chr2:57907740-57908150"),
                  geneCenteredSpec=list())
   return(cfSpec)

} # create.vrk2.candidateFilterSpec
#----------------------------------------------------------------------------------------------------
#  rs13384219  A->G
#  rs13384219 at chr2:57907073-57907573
#  gtcagtagtggtggaaccagcatgc[A/G]aattagacaatgtgacttcatagcc
#  Chromosome: 2:57907323
#  vrk2 tss at chr2:57908651
create.vrk2.rs13384219.neighborhood.candidateFilterSpec <- function(shoulder=10)
{
   spec <- create.vrk2.candidateFilterSpec(geneCentered=FALSE)

   rs13384219.loc <- 57907323
   left.loc <- rs13384219.loc - shoulder
   right.loc <- rs13384219.loc + shoulder
   spec$regionsSpec <- sprintf("chr2:%d-%d", left.loc, right.loc)

   spec

} # create.vrk2.rs13384219.neighborhood.candidateFilterSpec
#----------------------------------------------------------------------------------------------------
test_create.vrk2.candidateFilterSpec <- function()
{
   printf("--- test_create.vrk2.candidateFilterSpec")
   spec.0 <- create.vrk2.candidateFilterSpec(geneCentered=TRUE, promoter.length=1000)
   spec.1 <- create.vrk2.candidateFilterSpec(geneCentered=FALSE)

   checkEquals(spec.0$filterType, "EncodeDNaseClusters")
   checkEquals(spec.0$genomeName, "hg38")
   checkEquals(spec.0$encodeTableName, "wgEncodeRegDnaseClustered")
   #checkEquals(spec.0$fimoDB, "postgres://whovian/fimo")
   checkEquals(spec.0$geneInfoDB, "postgres://whovian/gtf")
   checkTrue(is.na(spec.0$regionsSpec))
   checkEquals(spec.0$geneCenteredSpec$targetGene, "VRK2")
   checkEquals(spec.0$geneCenteredSpec$tssUpstream, 1000)
   checkEquals(spec.0$geneCenteredSpec$tssDownstream, 1000)

   checkEquals(spec.1$filterType, "EncodeDNaseClusters")
   checkEquals(spec.1$genomeName, "hg38")
   checkEquals(spec.1$encodeTableName, "wgEncodeRegDnaseClustered")
   #checkEquals(spec.1$fimoDB, "postgres://whovian/fimo")
   checkEquals(spec.1$geneInfoDB, "postgres://whovian/gtf")
   checkEquals(spec.1$regionsSpec, "chr2:57906651-57908651")
   checkEquals(spec.1$geneCenteredSpec, list())

} # test_create.vrk2.candidateFilterSpec
#------------------------------------------------------------------------------------------------------------------------
test_basicConstructor <- function(reuse=FALSE)
{

   if(!reuse)  printf("--- test_basicConstructor")

   candidateFilterSpec <- create.vrk2.candidateFilterSpec()
   checkEquals(candidateFilterSpec$geneCenteredSpec$targetGene, "VRK2")
   checkEquals(candidateFilterSpec$geneCenteredSpec$tssUpstream, 1000)
   checkEquals(candidateFilterSpec$geneCenteredSpec$tssDownstream, 1000)
   checkTrue(is.na(candidateFilterSpec$regionsSpec))

   hdf <- with(candidateFilterSpec,
               HumanDHSFilter(genomeName,
                              encodeTableName=encodeTableName,
                              #fimoDatabase.uri=fimoDB,
                              pwmMatchPercentageThreshold=85L,
                              geneInfoDatabase.uri=geneInfoDB,
                              regionsSpec=regionsSpec,
                              geneCenteredSpec=geneCenteredSpec))

   checkTrue(all(c("HumanDHSFilter", "CandidateFilter") %in% is(hdf)))

      # make sure an unsupported genome triggers an error
   candidateFilterSpec$genomeName <- "intentional error in genome name"

   if(!reuse)
     checkException(hdf <- with(candidateFilterSpec,
                           HumanDHSFilter(genomeName,
                              encodeTableName=encodeTableName,
                              pwmMatchPercentageThreshold=85L,
                              geneInfoDatabase.uri=geneInfoDB,
                              regionsSpec=regionsSpec,
                              geneCenteredSpec=geneCenteredSpec)))
   if(reuse)
      return(hdf)

} # test_basicConstructor
#----------------------------------------------------------------------------------------------------
test_getEncodeRegulatoryTableNames <- function()
{
   printf("--- test_getEncodeRegulatoryTableNames")

   candidateFilterSpec <- create.vrk2.candidateFilterSpec()
   hdf <- with(candidateFilterSpec,
               HumanDHSFilter(genomeName,
                              encodeTableName=encodeTableName,
                              #fimoDatabase.uri=fimoDB,
                              pwmMatchPercentageThreshold=85L,
                              geneInfoDatabase.uri=geneInfoDB,
                              regionsSpec=regionsSpec,
                              geneCenteredSpec=geneCenteredSpec))

    names <- getEncodeRegulatoryTableNames(hdf)
    checkTrue(length(names) > 90)   # 96 on (13 apr 2017)

} # test_getEncodeRegulatoryTableNames
#----------------------------------------------------------------------------------------------------
test_checkAllEncodeTables <- function(quiet=TRUE)
{
   printf("--- test_checkAllEncodeTables")

   candidateFilterSpec <- create.vrk2.candidateFilterSpec()
   #candidateFilterSpec <- create.vrk2.rs13384219.neighborhood.candidateFilterSpec()

   hdf <- with(candidateFilterSpec,
               HumanDHSFilter(genomeName,
                              encodeTableName=encodeTableName,
                              #fimoDatabase.uri=fimoDB,
                              pwmMatchPercentageThreshold=85L,
                              geneInfoDatabase.uri=geneInfoDB,
                              regionsSpec=regionsSpec,
                              geneCenteredSpec=geneCenteredSpec,
                              quiet=TRUE))

   tableNames <- getEncodeRegulatoryTableNames(hdf)

   chrom <- "chr5"
   start <- 8800000
   end   <- 8850000

   #  rs13384219 at chr2:57907073-57907573
   chrom <- "chr2"
   rs13384219.loc <- 57907323
   start <- rs13384219.loc - 10
   end <- rs13384219.loc + 10

   for(tableName in tableNames){
      tbl <-getRegulatoryRegions(hdf, tableName, chrom, start, end)
      if(!quiet) printf("--- %s: %d rows", tableName, nrow(tbl))
      checkTrue(nrow(tbl) >= 0)
      checkEquals(colnames(tbl), c("chrom", "chromStart", "chromEnd",  "count",  "score"))
      #browser(); xyz <- 99
      }

} # test_checkAllEncodeTables
#----------------------------------------------------------------------------------------------------
# use this sample code to poke at the encode data offered by uscs
# note that most of the tables here only serve to list, not regions, but
# metadata:  what the inputs where, where the bb (bigBed) file can be found.
# for instance, the Hotspot table has this single line:
#                                                                     fileName
# 1 /gbdb/hg38/bbi/wgEncodeRegDnase/wgEncodeRegDnaseUwA549Hotspot.broadPeak.bb
explore.ucsc.database <- function()
{
   library(RMySQL)
   driver <- MySQL()
   host <- "genome-mysql.cse.ucsc.edu"
   user <- "genome"
   dbname <- "hg38"

   db <- dbConnect(driver, user = user, host = host, dbname = dbname)
   tables <- c("wgEncodeRegDnaseClustered", "wgEncodeRegDnaseUwA549Hotspot",  "wgEncodeRegDnaseUwA549Peak")
   main.clause <- sprintf("select * from %s where", tables[1]);

   chrom <- "chr5"
   start <- 88819630
   end   <- 88835936

   query <- paste(main.clause,
                  sprintf("chrom = '%s'", chromosome),
                   sprintf("and chromStart >= %d", start),
                   sprintf("and chromEnd <= %d", end),
                   collapse = " ")
   suppressWarnings(dbGetQuery(db, sprintf("select * from %s limit 5", tables[3])))


} # explore.ucsc.database
#----------------------------------------------------------------------------------------------------
test_geneSymbolToTSS <- function()
{
   printf("--- test_geneSymbolToTSS")
   candidateFilterSpec <- create.vrk2.candidateFilterSpec(geneCentered=TRUE)
   hdf <- with(candidateFilterSpec,
               HumanDHSFilter(genomeName,
                              encodeTableName=encodeTableName,
                              #fimoDatabase.uri=fimoDB,
                              pwmMatchPercentageThreshold=85L,
                              geneInfoDatabase.uri=geneInfoDB,
                              geneCenteredSpec=geneCenteredSpec,
                              regionsSpec=regionsSpec))
   x <- geneSymbolToTSS(hdf)
   checkEquals(x$chrom, "chr2")
   checkEquals(x$tss, 57907651)

} # test_geneSymbolToTSS
#----------------------------------------------------------------------------------------------------
test_getRegulatoryRegions <- function()
{
   printf("--- test_getRegulatoryRegions");

   candidateFilterSpec <- create.vrk2.candidateFilterSpec(geneCentered=FALSE)
   hdf <- with(candidateFilterSpec,
               HumanDHSFilter(genomeName,
                              encodeTableName=encodeTableName,
                              pwmMatchPercentageThreshold=85L,
                              geneInfoDatabase.uri=geneInfoDB,
                              geneCenteredSpec=geneCenteredSpec,
                              regionsSpec=regionsSpec))

   tableNames <- getEncodeRegulatoryTableNames(hdf)
   table <- "wgEncodeRegDnaseClustered"
   checkTrue(table %in% tableNames)
   x <- TReNA:::.parseChromLocString(candidateFilterSpec$regionsSpec)

   tbl.regions <-getRegulatoryRegions(hdf, table, x$chrom, x$start, x$end)

   checkTrue(nrow(tbl.regions) >= 5)
   checkEquals(colnames(tbl.regions), c("chrom", "chromStart", "chromEnd", "count", "score"))
   checkTrue(all(tbl.regions$chrom == x$chrom))
   checkTrue(all(tbl.regions$chromStart >= x$start))
   checkTrue(all(tbl.regions$chromStart <= x$end))
   checkTrue(all(tbl.regions$chromEnd >= x$start))
   checkTrue(all(tbl.regions$chromEnd <= x$end))

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
   tbl.regions <- getRegulatoryRegions(hdf, table, chrom, start, end)
   checkEquals(nrow(tbl.regions), 1)
   checkEquals(tbl.regions$chromStart, start)
   checkEquals(tbl.regions$chromEnd, end)

     # a small region, overhanging that DHS region at its upperbound, on the right
   chrom = "chr1"
   start <-  88802145
   end   <-  88802155
   tbl.regions <- getRegulatoryRegions(hdf, table, chrom, start, end)
   checkEquals(nrow(tbl.regions), 1)
   checkEquals(tbl.regions$chromStart, start)
   checkEquals(tbl.regions$chromEnd, 88802150)

     # a small region, overhanging that DHS region at its lowerbound, on the left
   chrom = "chr1"
   start <-  88801878
   end   <-  88801883
   tbl.regions <- getRegulatoryRegions(hdf, table, chrom, start, end)
   checkEquals(nrow(tbl.regions), 1)
   checkEquals(tbl.regions$chromStart, 88801880)
   checkEquals(tbl.regions$chromEnd, end)

      # another completely contained-in-DHS-region, small area of interest
   chrom <- "chr1"
   start <- 167830160
   end   <- 167830180
   tbl.regions <- getRegulatoryRegions(hdf, table, chrom, start, end)
   checkEquals(nrow(tbl.regions), 1)
   checkEquals(tbl.regions$chromStart, start)
   checkEquals(tbl.regions$chromEnd, end)

} # test_getRegulatoryRegions
#----------------------------------------------------------------------------------------------------
# test_getRegulatoryRegions_hardCase <- function()
# {
#    printf("--- test_getRegulatoryRegions_hardCase")
#
#      # a hard case.  the requested region overlaps incompletely with two neighboring
#      # dhs regions.  we want to get those two partial intersecting regions,
#      # one shortened on the left, the other shortened on the right
#
#    dhsFilter <- HumanDHSFilter("hg38", quiet=TRUE);
#
#    chrom <- "chr2"
#    start <-  127107283
#    end   <-  127107637
#
#        # dhs region 1: chr2:127,107,061 - 470
#        # dhs region 2:              620 - 929
#        # we expect back:
#        #      107,283 - 470
#        #      107,620 - 637
#
#     tbl.reg <- getRegulatoryRegions(dhsFilter, "wgEncodeRegDnaseClustered", chrom, start, end)
#     checkEquals(dim(tbl.reg), c(2, 5))
#     checkEquals(tbl.reg$chrom, c(chrom, chrom))
#     checkEquals(tbl.reg$chromStart, c(start, 127107620))
#     checkEquals(tbl.reg$chromEnd, c(127107470, end))
#
#        # now skip the second partial dhs region, get only the first
#     end   <-  127107500
#     tbl.reg <- getRegulatoryRegions(dhsFilter, "wgEncodeRegDnaseClustered", chrom, start, end)
#     checkEquals(nrow(tbl.reg), 1)
#     checkEquals(tbl.reg$chromStart, start)
#     checkEquals(tbl.reg$chromEnd, 127107470)
#
#        # now place start and end BOTH within the dhs region
#     end <- 127107465
#     tbl.reg <- getRegulatoryRegions(dhsFilter, "wgEncodeRegDnaseClustered", chrom, start, end)
#     checkEquals(nrow(tbl.reg), 1)
#     checkEquals(tbl.reg$chromStart, start)
#     checkEquals(tbl.reg$chromEnd, end)
#
# } # test_getRegulatoryRegions_hardCase
#----------------------------------------------------------------------------------------------------
test_getSequence <- function()
{
   printf("--- test_getSequence")
   candidateFilterSpec <- create.vrk2.candidateFilterSpec(geneCentered=FALSE)
   hdf <- with(candidateFilterSpec,
               HumanDHSFilter(genomeName,
                              encodeTableName=encodeTableName,
                              #fimoDatabase.uri=fimoDB,
                              pwmMatchPercentageThreshold=85L,
                              geneInfoDatabase.uri=geneInfoDB,
                              geneCenteredSpec=geneCenteredSpec,
                              regionsSpec=regionsSpec))



   chroms <- rep("chr5", 3)
   starts <- c(88819700, 88820700, 88820980)
   ends   <- c(88819910, 88820850, 88821130)

   tbl.regions <- data.frame(chrom=chroms, chromStart=starts, chromEnd=ends, stringsAsFactors=FALSE)
   seqs <- getSequence(hdf, tbl.regions)

   expected.lengths <- 1 + ends - starts
   checkEquals(unlist(lapply(seqs, nchar)), expected.lengths)
   invisible(seqs)

} # test_getSequence
#----------------------------------------------------------------------------------------------------
test_.matchForwardAndReverse <- function()
{
   printf("--- test_.matchForwardAndReverse")

   candidateFilterSpec <- create.vrk2.candidateFilterSpec(geneCentered=FALSE)
   hdf <- with(candidateFilterSpec,
               HumanDHSFilter(genomeName,
                              encodeTableName=encodeTableName,
                              #fimoDatabase.uri=fimoDB,
                              pwmMatchPercentageThreshold=85L,
                              geneInfoDatabase.uri=geneInfoDB,
                              geneCenteredSpec=geneCenteredSpec,
                              regionsSpec=regionsSpec))

   chrom <- "chr1"
   start <- 167829960
   end   <- 167830230

   motifName <- "MA0476.1"
   mtx <- query(MotifDb, motifName)[[1]];

   tbl.regions <- data.frame(chrom=chrom, chromStart=start, chromEnd=end, stringsAsFactors=FALSE)
   sequence <- getSequence(hdf, tbl.regions)

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

   tbl <- TReNA:::.matchForwardAndReverse(sequence, mtx, motifName, min.match.percentage=90, quiet=TRUE)
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
      # for easy comparison, round off the scores
   tbl$score         <- round(tbl$score, 2)
   tbl$maxScore      <- round(tbl$maxScore, 2)
   tbl$relativeScore <- round(tbl$relativeScore, 2)
   checkEquals(as.list(tbl), list(start=204, end=214, width=11, score=7.94, maxScore=8.35,
                                  relativeScore=0.95, motif="MA0478.1", match="CTGAGTCATGC",
                                  strand="-"))

        # now relax the score threshold
   tbl <- TReNA:::.matchForwardAndReverse(sequence, mtx, motifName, min.match.percentage=80, quiet=TRUE)
   tbl$score         <- round(tbl$score, 2)
   tbl$maxScore      <- round(tbl$maxScore, 2)
   tbl$relativeScore <- round(tbl$relativeScore, 2)
   checkEquals(tbl$start, c(56, 204))
   checkEquals(tbl$end,   c(66, 214))
   checkEquals(tbl$score, c(7.31, 7.94))
   checkEquals(tbl$maxScore, c(8.35, 8.35))
   checkEquals(tbl$relativeScore, c(0.88, 0.95))
   checkEquals(tbl$motif, c("MA0478.1", "MA0478.1"))
   checkEquals(tbl$match, c("AGCTGAGTCAT", "CTGAGTCATGC"))
   checkEquals(tbl$strand, c("+", "-"))


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
test_vrk2Promoter.incrementally <- function()
{
   printf("--- test_vrk2Promoter.incrementally")

    # chr5:88,813,245-88,832,344
   chrom <- "chr5"
   start <- 88824500
   end   <- 88832344

   df <- HumanDHSFilter("hg38");

   table <- "wgEncodeRegDnaseClustered"
   tbl.regions <- getRegulatoryRegions(df, table, chrom, start, end)
   checkEquals(dim(tbl.regions), c(18, 5))
   checkTrue(all(tbl.regions$chrom == chrom))
   checkTrue(all(tbl.regions$chromStart >= start))
   checkTrue(all(tbl.regions$chromEnd <= end))

   tbl.regions <- subset(tbl.regions, score >= 700)
   checkEquals(nrow(tbl.regions), 2)
   seqs <-getSequence(df, tbl.regions)
   checkEquals(with(tbl.regions, 1 + chromEnd - chromStart), nchar(seqs))  #  231 391

   motifs <- TReNA:::.getScoredMotifs(seqs, min.match.percentage=95)
   checkEquals(unlist(lapply(motifs, dim)), c(8,9,8,9))   # both

} # test_vrk2Promoter.incrementally
#----------------------------------------------------------------------------------------------------
test_getCandidates.vrk2.twoRegions <- function()
{
   printf("--- test_getCandidates.vrk2.twoRegions")

   cfSpec <- create.vrk2.candidateFilterSpec.twoRegions()
   hdf <- with(cfSpec, HumanDHSFilter(genomeName,
                                      encodeTableName=encodeTableName,
                                      pwmMatchPercentageThreshold=85L,
                                      geneInfoDatabase.uri=geneInfoDB,
                                      regionsSpec=regionsSpec,
                                      geneCenteredSpec=geneCenteredSpec,
                                      quiet=FALSE))

   x <- getCandidates(hdf)
   checkEquals(sort(names(x)), c("tbl", "tfs"))
   checkEquals(ncol(x$tbl), 11)
   checkTrue(nrow(x$tbl) > 150)    # 153 x 11
   checkTrue(length(x$tfs) > 400)  # 422

} # test_getCandidates.vrk2.twoRegions
#----------------------------------------------------------------------------------------------------
#  rs13384219  A->G
#  gtcagtagtggtggaaccagcatgc[A/G]aattagacaatgtgacttcatagcc
#  Chromosome: 2:57907323
#  vrk2 tss at chr2:57908651
#  from wgEncodeRegDnaseClustered:
#   chrom chromStart chromEnd count score
# 9  chr2   57907313 57907333     1    98

test_getCandidates.vrk2.rs13384219.neighborhood <- function()
{
   printf("--- test_getCandidates.vrk2.rs13384219.neighborhood")
   candidateFilterSpec <- create.vrk2.rs13384219.neighborhood.candidateFilterSpec()

   hdf <- with(candidateFilterSpec,
               HumanDHSFilter(genomeName,
                              encodeTableName=encodeTableName,
                              #fimoDatabase.uri=fimoDB,
                              pwmMatchPercentageThreshold=85L,
                              geneInfoDatabase.uri=geneInfoDB,
                              regionsSpec=regionsSpec,
                              geneCenteredSpec=geneCenteredSpec,
                              quiet=TRUE))
   x <- getCandidates(hdf)
   checkEquals(sort(names(x)), c("tbl", "tfs"))
   checkEquals(dim(x$tbl), c(66, 11))
   checkTrue(length(x$tfs) > 150)

} # test_getCandidates.vrk2.rs13384219.neighborhood
#------------------------------------------------------------------------------------------------------------------------
test_vrk2Promoter.normalUse <- function()
{
   printf("--- test_mef2cPromoter.normalUse")

#   hdcf <- HumanDHSFilter("hg38")
#    # chr5:88,813,245-88,832,344: has just a few high scoring clusters
#   chrom <- "chr5"
#   start <- 88824500
#   end   <- 88832344
#
#   filter.args <- list(chrom=chrom, start=start, end=end,
#                       region.score.threshold=700,
#                       motif.min.match.percentage=95,
#                       encode.table.name="wgEncodeRegDnaseClustered")
#
#   x <- getCandidates(hdcf, filter.args)
#   checkEquals(sort(names(x)), c("tbl", "tfs"))
#   checkTrue(all(c("chrom", "regionStart", "regionEnd", "regionScore", "sourceCount", "motif", "match",
#                   "motif.start", "motif.end", "motif.width", "motif.score", "strand", "tf") %in% colnames(x$tbl)))
#   checkTrue(nrow(x$tbl) > 15)     # 16 on (29 mar 2017)
#   checkTrue(length(x$tfs) > 200)  # 217 on (29 mar 2017)
#   checkEquals(length(which(duplicated(x$tfs))), 0)

} # test_vrk2Promoter.normalUse
#----------------------------------------------------------------------------------------------------
# a gene possibly involved in alzheimer's disease
notest_bin1 <- function()
{
   printf("--- notest_bin1")
   target.gene <- "BIN1"
   chrom <- "chr2"
   start <- 127107280
   end   <- 127107319
   start <- 127107538
   end   <- 127108337

   start <- 127107783
   end   <- 127107856


   dhs.args <- list(chrom=chrom,
                    start=start,
                    end=end,
                    region.score.threshold=0,
                    motif.min.match.percentage=85,
                    encode.table.name="wgEncodeRegDnaseClustered")

   dhsFilter <- HumanDHSFilter("hg38")
   dhs.out <- getCandidates(dhsFilter, dhs.args)

   genome.db.uri    <- "postgres://whovian/hg38"             # has gtf and motifsgenes tables
   footprint.db.uri <- "postgres://whovian/brain_hint"       # has hits and regions tables


   region <- sprintf("%s:%d-%d", chrom, start, end)
   recipe <- list(genomeDB=genome.db.uri,
                  footprintDB=footprint.db.uri,
                  geneCenteredSpec=list(),
                  regionsSpec=list(region))

   filter <- FootprintFilter(recipe$genomeDB,
                             recipe$footprintDB,
                             recipe$geneCenteredSpec,
                             recipe$regionsSpec,
                             quiet=TRUE)

   list.out <- getCandidates(filter)


   fpFilter <- FootprintFilter()
   fpf.out <- getCandidates(fpFilter, fpf.args)

      # the footprint finder, using fimo, finds more motifs than the dhsFilter, using Biostrings
      # explore.  first identify the motifs NOT found by the dhsFilter
   motifsNotFound <- setdiff(fpf.out$tbl$motifName, dhs.out$tbl$motif)
   tbl.notFound  <- subset(fpf.out$tbl, motifName %in% motifsNotFound)
   tbl.notFound$sequence <- getSequence(dhsFilter,

    with(tbl.notFound, data.frame(chrom=chrom, chromStart=start, chromEnd=end)))

    #library(FimoClient)
    #FIMO_HOST <- "whovian"
    #FIMO_PORT <- 5558
    #fc <<- FimoClient(FIMO_HOST, FIMO_PORT, quiet=TRUE)
    #requestMatch(fc, list(x="GGGGCGCGCGC"))

} # notest_bin1
#----------------------------------------------------------------------------------------------------
# one of the igap gwas alzheimer's snps, hand-verified to fall within a motif in a dhs region
test_rs34423320 <- function()
{
   printf("--- test_rs34423320")
     # chr1 167830170 167830170  rs34423320-C-T

   load(system.file(package="TReNA", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
   checkTrue(exists("mtx.sub"))
   hdcf <- HumanDHSFilter("hg38", mtx.sub)

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
   checkTrue(length(x$tfs) == 25)
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
