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
   test_getFootprintsForGeneSymbol()

   #test_getFootprintsForEnsemblGenes()

} # runTests
#----------------------------------------------------------------------------------------------------
test_.parseDatabaseUri <- function()
{
   printf("--- test_.parseDatabaseUri")
   genome.db.uri <- "postgres://whovian/hg38"
   project.db.uri <-  "postgres://whovian/lymphoblast"

   x <- TReNA:::.parseDatabaseUri(genome.db.uri)
   checkEquals(x$brand, "postgres")
   checkEquals(x$host,  "whovian")
   checkEquals(x$name,  "hg38")

   x <- TReNA:::.parseDatabaseUri(project.db.uri)
   checkEquals(x$brand, "postgres")
   checkEquals(x$host,  "whovian")
   checkEquals(x$name,  "lymphoblast")

} # test_.parseDatabaseUri
#----------------------------------------------------------------------------------------------------
test_constructor <- function()
{
   printf("--- test_constructor")

   genome.db.uri <- "postgres://whovian/hg38"
   project.db.uri <-  "postgres://whovian/lymphoblast"
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
   project.db.uri <-  "postgres://whovian/lymphoblast"
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
   project.db.uri <-  "postgres://whovian/lymphoblast"
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
   project.db.uri <-  "postgres://whovian/lymphoblast"
   fp <- FootprintFinder(genome.db.uri, project.db.uri, quiet=TRUE)
   tbl.loc <- getChromLoc(fp, "MEF2C", biotype="protein_coding", moleculetype="gene")
   checkEquals(dim(tbl.loc), c(1, 5))
   expected <- list(gene_name="MEF2C",
                    chr="chr5",
                    start=88717117,
                    endpos=88904257,
                    strand="-")
   checkEquals(as.list(tbl.loc), expected)
   tbl.loc <- getChromLoc(fp, "MEF2C")
   checkEquals(nrow(tbl.loc), 1)
   checkEquals(as.list(tbl.loc), expected)
   closeDatabaseConnections(fp)

} # test_getChromLoc
#------------------------------------------------------------------------------------------------------------------------
test_getGenePromoterRegion <- function()
{
   printf("--- test_getGenePromoterRegion")

   genome.db.uri <- "postgres://whovian/hg38"
   project.db.uri <-  "postgres://whovian/lymphoblast"
   fp <- FootprintFinder(genome.db.uri, project.db.uri, quiet=TRUE)
   region <- getGenePromoterRegion(fp, "TREM2", 0, 0)

   checkEquals(region$chr, "chr6")
   checkEquals(region$start, 41163186)
   checkEquals(region$end,   41163186)

   region <- getGenePromoterRegion(fp, "TREM2", 20, 30)
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
   closeDatabaseConnections(fp)

} # test_getGenePromoterRegion
#----------------------------------------------------------------------------------------------------
test_getFootprintsForGeneSymbol <- function()
{
   printf("--- test_getFootprintsForGeneSymbol")

   genome.db.uri <- "postgres://whovian/hg38"
   project.db.uri <-  "postgres://whovian/wholeBrain"
   fp <- FootprintFinder(genome.db.uri, project.db.uri, quiet=TRUE)

      # use MEF2C and the hg38 assembly
   tbl.fp <- getFootprintsForGeneSymbol(fp, "MEF2C", size.upstream=0, size.downstream=0)
   checkEquals(dim(tbl.fp), c(0, 0))

         # 3k up and downstream.  we expect more footprints upstream,  some downstream
   tbl <- getFootprintsForGeneSymbol(fp, "MEF2C", size.upstream=3000,  size.downstream=1000)
   checkTrue(nrow(tbl) > 50)   #69
   closeDatabaseConnections(fp)

} # test_getFootprintsForGeneSymbol
#----------------------------------------------------------------------------------------------------
test_getFootprintsInRegion <- function()
{
   printf("--- test_getFootprintsInRegion")

   genome.db.uri <- "postgres://whovian/hg38"
   project.db.uri <-  "postgres://whovian/wholeBrain"
   fp <- FootprintFinder(genome.db.uri, project.db.uri, quiet=TRUE)

      # use MEF2C and the hg38 assembly
   chromosome <- "chr5"
   tss <- 88904257

       # a region of size 1.  no footprints here
   tbl.fp <- getFootprintsInRegion(fp, chromosome, tss, tss)
   checkEquals(dim(tbl.fp), c(0, 0))

         # 3k up and downstream.  we expect more footprints upstream,  some downstream
   tbl <- getFootprintsInRegion(fp, chromosome, tss, tss + 3000)
   checkTrue(nrow(tbl) > 50)   # 59 before 10may16, 257 after
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
   db.uri <- sprintf("sqlite:%s", system.file(package="TReNA.brain", "extdata", "fpTf.sqlite"))
   genes <- c("ENSG00000267051", "ENSG00000264503", "ENSG00000273141", "ENSG00000212712",
              "ENSG00000236396", "ENSG00000154889", "ENSG00000267794",  "ENSG00000264843",
              "ENSG00000260759", "ENSG00000154856")
   fp <- FootprintFinder(ensembl.hg38, db.uri, genes)
   checkEquals(getGenes(fp), genes)
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
