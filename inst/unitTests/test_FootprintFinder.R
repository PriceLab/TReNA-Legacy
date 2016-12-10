library(TReNA)
library(RUnit)
library(RPostgreSQL)
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
   
   #test_getFootprintsForEnsemblGenes()

} # runTests
#----------------------------------------------------------------------------------------------------
test_.parseDatabaseUri <- function()
{
   printf("--- test_.parseDatabaseUri")
   genome.db.uri <- "postgres://whovian/hg38"
   #project.db.uri <-  "postgres://whovian/lymphoblast"
   project.db.uri <-  "postgres://whovian/brain_wellington"

   x <- TReNA:::.parseDatabaseUri(genome.db.uri)
   checkEquals(x$brand, "postgres")
   checkEquals(x$host,  "whovian")
   checkEquals(x$name,  "hg38")

   x <- TReNA:::.parseDatabaseUri(project.db.uri)
   checkEquals(x$brand, "postgres")
   checkEquals(x$host,  "whovian")
   checkEquals(x$name,  "brain_wellington")

} # test_.parseDatabaseUri
#----------------------------------------------------------------------------------------------------
test_constructor <- function()
{
   printf("--- test_constructor")

   genome.db.uri <- "postgres://whovian/hg38"
   project.db.uri <-  "postgres://whovian/brain_wellington"
   fp <- FootprintFinder(genome.db.uri, project.db.uri, quiet=TRUE)
   closeDatabaseConnections(fp)

} # test_constructor
#----------------------------------------------------------------------------------------------------
test_database.hg38.whovian <- function()
{
   printf("--- test_database.hg38.whovian")
   db <- dbConnect(PostgreSQL(), user= "trena", password="trena", dbname="hg38", host="whovian")
   checkTrue("DBIConnection" %in% is(db))
   checkTrue("gtf" %in% dbListTables(db))
   rowCount <-  dbGetQuery(db, "select count(*) from gtf")[1, 1]
   checkTrue(rowCount >  2.5 * 10^6)

} # test_database.hg38.whovian
#----------------------------------------------------------------------------------------------------
test_getGtfGeneBioTypes <- function()
{
   printf("--- test_getGtfGeneBioTypes")

   genome.db.uri <- "postgres://whovian/hg38"
   project.db.uri <-  "postgres://whovian/brain_wellington"
   fp <- FootprintFinder(genome.db.uri, project.db.uri, quiet=TRUE)
   types <- getGtfGeneBioTypes(fp)
   checkTrue(length(types) >= 40)
   some.expected.types <- c("processed_transcript", "protein_coding", "pseudogene", "rRNA", "ribozyme", "sRNA")
   checkTrue(all(some.expected.types %in% types))
   closeDatabaseConnections(fp)


} # test_getGtfGeneBioTypes
#----------------------------------------------------------------------------------------------------
test_getGtfMoleculeTypes <- function()
{
   printf("--- test_getGtfMoleculeTypes")

   genome.db.uri <- "postgres://whovian/hg38"
   project.db.uri <-  "postgres://whovian/brain_wellington"
   fp <- FootprintFinder(genome.db.uri, project.db.uri, quiet=TRUE)
   types <- getGtfMoleculeTypes(fp)

   checkTrue(length(types) >= 9)
   some.expected.types <- c("CDS", "exon", "five_prime_utr", "gene", "start_codon", "stop_codon",
                            "three_prime_utr", "transcript")
   checkTrue(all(some.expected.types %in% types))
   closeDatabaseConnections(fp)


} # test_getGtfMoleculeTypes
#----------------------------------------------------------------------------------------------------
test_getChromLoc <- function()
{
   printf("--- test_getChromLoc")

   genome.db.uri <- "postgres://whovian/hg38"
   project.db.uri <-  "postgres://whovian/brain_wellington"
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

   genome.db.uri <- "postgres://whovian/hg38"
   project.db.uri <-  "postgres://whovian/brain_wellington"
   fp <- FootprintFinder(genome.db.uri, project.db.uri, quiet=TRUE)

      # TREM2: prepare for test of both gene symbol and ensembl gene id
   tbl.loc <- getChromLoc(fp, "TREM2")
   ensg <- tbl.loc$gene_id[1]

      # test this gene, two ids, zero length and 50 bp length regions

   region <- getGenePromoterRegion(fp, "TREM2", 0, 0)
   checkEquals(region$chr, "chr6")
   checkEquals(region$start, 41163186)
   checkEquals(region$end,   41163186)

   region <- getGenePromoterRegion(fp, ensg, 0, 0)
   checkEquals(region$chr, "chr6")
   checkEquals(region$start, 41163186)
   checkEquals(region$end,   41163186)

   region <- getGenePromoterRegion(fp, "TREM2", 20, 30)
   checkEquals(region$chr, "chr6")
   checkEquals(region$start, 41163156)  #
   checkEquals(region$end,   41163206)  # 20 bases upstream from the "end" TSS

   region <- getGenePromoterRegion(fp, ensg, 20, 30)
   checkEquals(region$chr, "chr6")
   checkEquals(region$start, 41163156)  #
   checkEquals(region$end,   41163206)  # 20 bases upstream from the "end" TSS


     # now a plus strand gene
   region <- getGenePromoterRegion(fp, "SP1", 0, 0)
   checkEquals(region$chr, "chr12")
   checkEquals(region$start, 53380176)
   checkEquals(region$end,   53380176)

   region <- getGenePromoterRegion(fp, "SP1", 1000, 1)
   checkEquals(region$chr, "chr12")
   checkEquals(region$start, 53379176)
   checkEquals(region$end,   53380177)

     # now try a lincRNA, with explicit biotype
   region <- getGenePromoterRegion(fp, "LINC01254", 1000, 1000, biotype="lincRNA")
   checkEquals(region$chr, "chr18")
   checkEquals(region$start, 10413515)
   checkEquals(region$end,   10415515)

     # now try a lincRNA, with implicit biotype, "protein_coding"
   suppressWarnings(checkTrue(is.na(getGenePromoterRegion(fp, "LINC01254", 1000, 1000))))

   closeDatabaseConnections(fp)

} # test_getGenePromoterRegion
#----------------------------------------------------------------------------------------------------
test_getFootprintsForGene <- function()
{
   printf("--- test_getFootprintsForGene")

   genome.db.uri <- "postgres://whovian/hg38"
   project.db.uri <-  "postgres://whovian/brain_wellington"
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
   project.db.uri <-  "postgres://whovian/brain_hint"
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
# FootprintFinder originally accepted only HUGO gene symbols, which is still the expected case
# however, ensembl reports, and we currently have expression data for, a variety of DNA elements
# including miRNA, linkRNA, antisense genes, pseudogenes of various sorts, etc.
# these each have a unique ENSG id, which we test out here
# the constructor of FootprintFinder needs to recognise these identifiers, and make a corresponding
# call to biomart
test_getFootprintsForEnsemblGenes <- function()
{
   printf("--- test_getFootprintsForEnsemblGenes")

   genome.db.uri <- "postgres://whovian/hg38"
   project.db.uri <-  "postgres://whovian/brain_wellington"
   fp <- FootprintFinder(genome.db.uri, project.db.uri, quiet=TRUE)

   genes <- c("ENSG00000267051", "ENSG00000264503", "ENSG00000273141", "ENSG00000212712",
              "ENSG00000236396", "ENSG00000154889", "ENSG00000267794",  "ENSG00000264843",
              "ENSG00000260759", "ENSG00000154856")

   goi <- genes[3]
   loc <- getGenePromoterRegion(fp, goi, 250, 0)
   tbl <- getFootprints(fp, goi, 250, 0)
   checkTrue(all(tbl$mfpStart >= loc$start))
   checkTrue(all(tbl$mfpStart <= loc$end))
   checkTrue(all(tbl$mfpEnd >= loc$start))
   checkTrue(all(tbl$mfpEnd <= loc$end))

} # test_getFootprintsForEnsemblGenes
#----------------------------------------------------------------------------------------------------
if(!interactive()) runTests()
