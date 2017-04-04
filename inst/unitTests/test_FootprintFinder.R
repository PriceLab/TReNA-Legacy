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

   genome.db.uri <- "postgres://whovian/hg38"
   project.db.uri <-  "postgres://whovian/brain_wellington"
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
test_mapMotifsToTFsMergeIntoTable <- function()
{
    printf("--- test_mapMotifsToTFsMergeIntoTable")
    
   genome.db.uri <- "postgres://whovian/hg38"
   project.db.uri <-  "postgres://whovian/brain_wellington"

   fp <- FootprintFinder(genome.db.uri, project.db.uri, quiet=TRUE)

      # use MEF2C and the hg38 assembly
   chromosome <- "chr5"
   tss <- 88904257

         # 3k up and downstream.  we expect more footprints upstream,  some downstream
   tbl <- getFootprintsInRegion(fp, chromosome, tss-1000, tss + 1000)
   checkTrue(nrow(tbl) > 0)   # 11 

   tbl.withTFs <- mapMotifsToTFsMergeIntoTable(fp, tbl)
   checkEquals(ncol(tbl.withTFs), 1 + ncol(tbl))
   checkTrue(nrow(tbl.withTFs) > nrow(tbl))
   closeDatabaseConnections(fp)

} # test_getFootprintsInRegion
#----------------------------------------------------------------------------------------------------
if(!interactive()) runTests()
