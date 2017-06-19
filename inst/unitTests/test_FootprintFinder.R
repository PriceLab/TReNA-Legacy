library(TReNA)
library(RUnit)
library(RPostgreSQL)
#----------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_.parseDatabaseUri()
   test_constructor()
   test_getGtfGeneBioTypes()
   test_getGtfMoleculeTypes()
   test_getChromLoc()
   test_getGenePromoterRegion()
   test_getFootprintsInRegion()
   test_getFootprintsInRegionWithVariants()
   test_getFootprintsForGene()
   test_mapMotifsToTFsMergeIntoTable()

} # runTests
#----------------------------------------------------------------------------------------------------
test_.parseDatabaseUri <- function()
{
   printf("--- test_.parseDatabaseUri")
   db.address <- system.file(package="TReNA", "extdata")
   genome.db.uri <- paste("sqlite:/",db.address,"genome.sub.db", sep = "/")
   project.db.uri <- paste("sqlite:/",db.address,"project.sub.db", sep = "/")

   x <- TReNA:::.parseDatabaseUri(genome.db.uri)
   checkEquals(x$brand, "sqlite")
   checkEquals(x$host,  db.address)
   checkEquals(x$name,  "genome.sub.db")

   x <- TReNA:::.parseDatabaseUri(project.db.uri)
   checkEquals(x$brand, "sqlite")
   checkEquals(x$host,  db.address)
   checkEquals(x$name,  "project.sub.db")

} # test_.parseDatabaseUri
#----------------------------------------------------------------------------------------------------
test_constructor <- function()
{
   printf("--- test_constructor")

   db.address <- system.file(package="TReNA", "extdata")
   genome.db.uri <- paste("sqlite:/",db.address,"genome.sub.db", sep = "/")
   project.db.uri <- paste("sqlite:/",db.address,"project.sub.db", sep = "/")

   fp <- FootprintFinder(genome.db.uri, project.db.uri, quiet=TRUE)
   closeDatabaseConnections(fp)

} # test_constructor
#----------------------------------------------------------------------------------------------------
test_getGtfGeneBioTypes <- function()
{
   printf("--- test_getGtfGeneBioTypes")

   db.address <- system.file(package="TReNA", "extdata")
   genome.db.uri <- paste("sqlite:/",db.address,"genome.sub.db", sep = "/")
   project.db.uri <- paste("sqlite:/",db.address,"project.sub.db", sep = "/")

   fp <- FootprintFinder(genome.db.uri, project.db.uri, quiet=TRUE)
   types <- getGtfGeneBioTypes(fp)
#   checkTrue(length(types) >= 40)
#   some.expected.types <- c("processed_transcript", "protein_coding", "pseudogene", "rRNA", "ribozyme", "sRNA")
   #   checkTrue(all(some.expected.types %in% types))
   checkEquals(types, "protein_coding")
   closeDatabaseConnections(fp)

} # test_getGtfGeneBioTypes
#----------------------------------------------------------------------------------------------------
test_getGtfMoleculeTypes <- function()
{
   printf("--- test_getGtfMoleculeTypes")

   db.address <- system.file(package="TReNA", "extdata")
   genome.db.uri <- paste("sqlite:/",db.address,"genome.sub.db", sep = "/")
   project.db.uri <- paste("sqlite:/",db.address,"project.sub.db", sep = "/")

   fp <- FootprintFinder(genome.db.uri, project.db.uri, quiet=TRUE)
   types <- getGtfMoleculeTypes(fp)

   checkEquals(types, "gene")
   #checkTrue(length(types) >= 9)
   #some.expected.types <- c("CDS", "exon", "five_prime_utr", "gene", "start_codon", "stop_codon",
   #                         "three_prime_utr", "transcript")
   #checkTrue(all(some.expected.types %in% types))
   closeDatabaseConnections(fp)

} # test_getGtfMoleculeTypes
#----------------------------------------------------------------------------------------------------
test_getChromLoc <- function()
{
   printf("--- test_getChromLoc")

   db.address <- system.file(package="TReNA", "extdata")
   genome.db.uri <- paste("sqlite:/",db.address,"genome.sub.db", sep = "/")
   project.db.uri <- paste("sqlite:/",db.address,"project.sub.db", sep = "/")

   fp <- FootprintFinder(genome.db.uri, project.db.uri, quiet=TRUE)
   tbl.loc <- getChromLoc(fp, "MEF2C", biotype="protein_coding", moleculetype="gene")
   checkEquals(dim(tbl.loc), c(1, 6))
   expected <- list(gene_id="ENSG00000081189",
                    gene_name="MEF2C",
                    chr="chr5",
                    start=88717117,
                    endpos=88904257,
                    strand="-")
   checkEquals(as.list(tbl.loc), expected)
      # now use default values for biotype and moleculetype
   tbl.loc <- getChromLoc(fp, "MEF2C")
   checkEquals(nrow(tbl.loc), 1)
   checkEquals(as.list(tbl.loc), expected)

      # repeat with ensembl gene id
   tbl.loc <- getChromLoc(fp, "ENSG00000081189")
   checkEquals(nrow(tbl.loc), 1)
   checkEquals(as.list(tbl.loc), expected)
   closeDatabaseConnections(fp)

} # test_getChromLoc
#------------------------------------------------------------------------------------------------------------------------
test_getGenePromoterRegion <- function()
{
   printf("--- test_getGenePromoterRegion")

   db.address <- system.file(package="TReNA", "extdata")
   genome.db.uri <- paste("sqlite:/",db.address,"genome.sub.db", sep = "/")
   project.db.uri <- paste("sqlite:/",db.address,"project.sub.db", sep = "/")

   fp <- FootprintFinder(genome.db.uri, project.db.uri, quiet=TRUE)

      # TREM2: prepare for test of both gene symbol and ensembl gene id
   tbl.loc <- getChromLoc(fp, "MEF2C")
   ensg <- tbl.loc$gene_id[1]

      # test this gene, two ids, zero length and 50 bp length regions

   region <- getGenePromoterRegion(fp, "MEF2C", 0, 0)
   checkEquals(region$chr, "chr5")
   checkEquals(region$start, 88904257)
   checkEquals(region$end,   88904257)

   region <- getGenePromoterRegion(fp, ensg, 0, 0)
   checkEquals(region$chr, "chr5")
   checkEquals(region$start, 88904257)
   checkEquals(region$end,   88904257)

   region <- getGenePromoterRegion(fp, "MEF2C", 20, 30)
   checkEquals(region$chr, "chr5")
   checkEquals(region$start, 88904227)  #
   checkEquals(region$end,   88904277)  # 20 bases upstream from the "end" TSS

   region <- getGenePromoterRegion(fp, ensg, 20, 30)
   checkEquals(region$chr, "chr5")
   checkEquals(region$start, 88904227)  #
   checkEquals(region$end,   88904277)  # 20 bases upstream from the "end" TSS

     # now a plus strand gene
#   region <- getGenePromoterRegion(fp, "SP1", 0, 0)
#   checkEquals(region$chr, "chr12")
#   checkEquals(region$start, 53380176)
#   checkEquals(region$end,   53380176)

#   region <- getGenePromoterRegion(fp, "SP1", 1000, 1)
#   checkEquals(region$chr, "chr12")
#   checkEquals(region$start, 53379176)
#   checkEquals(region$end,   53380177)

     # now try a lincRNA, with explicit biotype
#   region <- getGenePromoterRegion(fp, "LINC01254", 1000, 1000, biotype="lincRNA")
#   checkEquals(region$chr, "chr18")
#   checkEquals(region$start, 10413515)
#   checkEquals(region$end,   10415515)

     # now try a lincRNA, with implicit biotype, "protein_coding"
#   suppressWarnings(checkTrue(is.na(getGenePromoterRegion(fp, "LINC01254", 1000, 1000))))

   closeDatabaseConnections(fp)

} # test_getGenePromoterRegion
#----------------------------------------------------------------------------------------------------
test_getFootprintsForGene <- function()
{
   printf("--- test_getFootprintsForGene")

   db.address <- system.file(package="TReNA", "extdata")
   genome.db.uri <- paste("sqlite:/",db.address,"genome.sub.db", sep = "/")
   project.db.uri <- paste("sqlite:/",db.address,"project.sub.db", sep = "/")

   fp <- FootprintFinder(genome.db.uri, project.db.uri, quiet=TRUE)

      # get enembl gene id for MEF2C
   tbl.loc <- getChromLoc(fp, "MEF2C")
   mef2c.ensg <- tbl.loc$gene_id[1]

      # use MEF2C and the hg38 assembly
   tbl.fp <- getFootprintsForGene(fp, "MEF2C", size.upstream=0, size.downstream=0)
   checkEquals(dim(tbl.fp), c(0, 0))
   tbl.fp <- getFootprintsForGene(fp, mef2c.ensg, size.upstream=0, size.downstream=0)
   checkEquals(dim(tbl.fp), c(0, 0))

         # 3k up and downstream.  we expect more footprints upstream,  some downstream
   tbl <- getFootprintsForGene(fp, "MEF2C", size.upstream=3000,  size.downstream=1000)
   expected.colnames <- c("loc", "chrom", "start", "endpos", "type", "name", "length", "strand",
                          "sample_id", "method",  "provenance", "score1", "score2", "score3",
                          "score4",  "score5",  "score6")
   checkTrue(all(expected.colnames %in% colnames(tbl)))
   checkTrue(nrow(tbl) > 20)   # 1385

   tbl.ensg <- getFootprintsForGene(fp, mef2c.ensg, size.upstream=3000,  size.downstream=1000)
   checkEquals(dim(tbl), dim(tbl.ensg))

   closeDatabaseConnections(fp)

} # test_getFootprintsForGene
#----------------------------------------------------------------------------------------------------
test_getFootprintsInRegion <- function()
{
   printf("--- test_getFootprintsInRegion")

   db.address <- system.file(package="TReNA", "extdata")
   genome.db.uri <- paste("sqlite:/",db.address,"genome.sub.db", sep = "/")
   project.db.uri <- paste("sqlite:/",db.address,"project.sub.db", sep = "/")

   fp <- FootprintFinder(genome.db.uri, project.db.uri, quiet=TRUE)

      # use MEF2C and the hg38 assembly
   chromosome <- "chr5"
   tss <- 88904257

       # a region of size 1.  no footprints here
   tbl.fp <- getFootprintsInRegion(fp, chromosome, tss, tss)
   checkEquals(dim(tbl.fp), c(0, 0))

         # 3k up and downstream.  we expect more footprints upstream,  some downstream
   tbl <- getFootprintsInRegion(fp, chromosome, tss, tss + 3000)
   checkTrue(nrow(tbl) > 0)   # 59 before 10may16, 257 after
   closeDatabaseConnections(fp)

} # test_getFootprintsInRegion
#----------------------------------------------------------------------------------------------------
test_getFootprintsInRegionWithVariants <- function()
{
   printf("--- test_getFootprintsInRegionWithVariants")

   db.address <- system.file(package="TReNA", "extdata")
   genome.db.uri <- paste("sqlite:/",db.address,"genome.sub.db", sep = "/")
   project.db.uri <- paste("sqlite:/",db.address,"project.sub.db", sep = "/")

   fp <- FootprintFinder(genome.db.uri, project.db.uri, quiet=TRUE)

      # use MEF2C and the hg38 assembly
   chromosome <- "chr5"
   region.start <-  88894577 - 10
   region.end   <-  88894583 + 10

       # a region of size 1.  no footprints here
   tbl.fp <- getFootprintsInRegion(fp, chromosome, region.start, region.end)
   checkEquals(dim(tbl.fp), c(2, 17))
   checkTrue("MA0152.1" %in% tbl.fp$name)

   tbl.regions.noSeq <- data.frame(chrom="chr5", start=region.start, end=region.end, stringsAsFactors=FALSE)
   mm <- MotifMatcher("tester", "hg38")
   tbl.regions <- getSequence(mm, tbl.regions.noSeq)
   checkEquals(tbl.regions$seq,  "AGGATGAATTTTTTCCAAAAGTAAATC")

   mm.wt <- findMatchesByChromosomalRegion(mm, tbl.regions, pwmMatchMinimumAsPercentage=80)
      # perfect match:
   checkEquals(subset(mm.wt$tbl, motifName=="MA0152.1")$motifRelativeScore, 1)
      # what other motifs are reported by MotifMatcher (which uses the bioc matchPWM, different from FIMO)?
   motifs.wt.strong <- sort(subset(mm.wt$tbl, motifRelativeScore > .90)$motifName)
   checkEquals(motifs.wt.strong, c("MA0109.1", "MA0152.1", "MA0606.1", "MA0624.1", "MA0625.1"))

   mm.mut <- findMatchesByChromosomalRegion(mm, tbl.regions, pwmMatchMinimumAsPercentage=90, "chr5:88894580:T:G")
                                         # c("chr5:88894577:T:A", "chr5:88894580:T:G"))
   motifs.mut.strong <- sort(subset(mm.mut$tbl, motifRelativeScore > .90)$motifName)

   motifs.lost <- sort(setdiff(motifs.wt.strong, motifs.mut.strong))
   motifs.gained <- sort(setdiff(motifs.mut.strong, motifs.wt.strong))
   motifs.preserved <- sort(intersect(motifs.mut.strong, motifs.wt.strong))

   checkEquals(motifs.lost, c("MA0152.1", "MA0606.1", "MA0624.1", "MA0625.1"))
   checkEquals(motifs.gained, c("MA0161.1", "MA0670.1", "MA0671.1"))
   checkEquals(motifs.preserved, "MA0109.1")

} # test_getFootprintsInRegionWithVariants
#----------------------------------------------------------------------------------------------------
test_mapMotifsToTFsMergeIntoTable <- function()
{
    printf("--- test_mapMotifsToTFsMergeIntoTable")

   db.address <- system.file(package="TReNA", "extdata")
   genome.db.uri <- paste("sqlite:/",db.address,"genome.sub.db", sep = "/")
   project.db.uri <- paste("sqlite:/",db.address,"project.sub.db", sep = "/")

   fp <- FootprintFinder(genome.db.uri, project.db.uri, quiet=TRUE)

      # use MEF2C and the hg38 assembly
   chromosome <- "chr5"
   tss <- 88904257

         # 3k up and downstream.  we expect more footprints upstream,  some downstream
   tbl <- getFootprintsInRegion(fp, chromosome, tss-1000, tss + 1000)
   checkTrue(nrow(tbl) > 0)   # 11

   tbl.withTFs <- mapMotifsToTFsMergeIntoTable(fp, tbl)
   checkEquals(ncol(tbl.withTFs), 10)
   checkEquals(nrow(tbl.withTFs),nrow(tbl))
   closeDatabaseConnections(fp)

} # test_getFootprintsInRegion
#----------------------------------------------------------------------------------------------------
explore_variantsInFootprints <- function()
{

   # how does   rs13384219  A->G
   #   gtcagtagtggtggaaccagcatgc[A/G]aattagacaatgtgacttcatagcc, chr2:57907323
   # affect motifs in footprints?
   # vrk2 tss: 57907651 + strand
   # 1 + 57907651 - 57907323: 329 bp upstream of tss

   library(FimoClient)
   if(!exists("fimo.service"))
      fimo.service <-  FimoClient("whovian", 5558, quiet=TRUE)
   library(BSgenome.Hsapiens.UCSC.hg38)
   hg38 = BSgenome.Hsapiens.UCSC.hg38
   library(SNPlocs.Hsapiens.dbSNP144.GRCh38)
   snps <- SNPlocs.Hsapiens.dbSNP144.GRCh38
   as.data.frame(snpsById(snps, "rs13384219"))

     #   seqnames      pos strand  RefSNP_id alleles_as_ambig
     #      ch2 57907323      + rs13384219                R
   IUPAC_CODE_MAP[["R"]]   # AG

   ranges <- GRanges(seqnames="ch2", IRanges(tss-400, tss))
   gr.snpsOv <- snpsByOverlaps(snps, ranges, maxgap=0L, minoverlap=0L,
                               type="any",  #c("any", "start", "end", "within", "equal"),
                               drop.rs.prefix=FALSE)

   genome.db.uri <- "postgres://whovian/hg38"
   project.db.uri <-  "postgres://whovian/brain_hint"
   fp <- FootprintFinder(genome.db.uri, project.db.uri, quiet=TRUE)
   tss <- 57907651
   snp.loc <- 57907323
   tbl <- getFootprintsInRegion(fp, chromosome, tss-400, tss - 300)
   checkTrue(nrow(tbl) > 0)   # 11

   tbl.withTFs <- mapMotifsToTFsMergeIntoTable(fp, tbl)

} # explore_variantsInFootprints
#----------------------------------------------------------------------------------------------------
#    library(FimoClient)
#    if(!exists("fimo.service"))
#       fimo.service <-  FimoClient("whovian", 5558, quiet=TRUE)
#    library(BSgenome.Hsapiens.UCSC.hg38)
#    hg38 = BSgenome.Hsapiens.UCSC.hg38
#    library(SNPlocs.Hsapiens.dbSNP144.GRCh38)
#    snps <- SNPlocs.Hsapiens.dbSNP144.GRCh38
#    rsid <- "rs6718960"
#    chrom <- "chr2"
#    loc <- 165239218
#    ambiguity.code <- snpsById(snps, rsid)$alleles_as_ambig
#    elements.string <- IUPAC_CODE_MAP[[ambiguity.code]]
#    elements <- strsplit(elements.string,'')[[1]]
#    wt <- as.character(getSeq(hg38, chrom, loc, loc))
#    mut <- setdiff(elements, wt)
#    doComparativeFimo(chrom, loc, wt, mut, 10, quiet=FALSE) # gain
#      X.pattern.name sequence.name start stop strand   score  p.value q.value matched.sequence
#            MA0040.1           mut     2   12      + 12.4388 3.88e-05 0.00171      AAATGTTTAGA
#
#    rsid <- "rs7596642"
#    chrom <- "chr2"
#    loc <- 241150035
#    ambiguity.code <- snpsById(snps, rsid)$alleles_as_ambig
#    elements.string <- IUPAC_CODE_MAP[[ambiguity.code]]
#    elements <- strsplit(elements.string,'')[[1]]
#    wt <- as.character(getSeq(hg38, chrom, loc, loc))
#    mut <- setdiff(elements, wt)
#    doComparativeFimo(chrom, loc, wt, mut, 10, quiet=FALSE) # noMotif
#
#
#------------------------------------------------------------------------------------------------------------------------
if(!interactive()) runTests()
