library(TReNA)
library(RUnit)
library(biomaRt)
library(RPostgreSQL)
#----------------------------------------------------------------------------------------------------
#http://uswest.ensembl.org/
#if(!exists("ensembl.hg38"))
#    ensembl.hg38 <- useMart("ensembl", dataset="hsapiens_gene_ensembl", host="www.ensembl.org")

#----------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_.parseDatabaseUri()
   test_constructor()
   test_getPromoterRegion()
   # test_privateWholeBrainFootprintDatabaseOnWhovian()
   # test_privateLymphoblastFootprintDatabaseOnWhovian()
   # test_privateWholeBrainFootprintDatabaseOnWhovian_manualJoin()
   # test_privateLymphoblastFootprintDatabaseOnWhovian_manualJoin()

   #test_.getGeneLocations()

   #test_getFootprints()
   #test_constructWithAllChromosome5genes()
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
# make sure that the gtf tables can be read and queried by the trena user
test_privateGtfDatabaseOnWhovian <- function()
{
   printf("--- test_privateGtfDatabaseOnWhovian")

   db <- dbConnect(PostgreSQL(), user= "trena", password="trena", dbname="gtf", host="whovian")

     # two tables at present: footprints and motifsGenes.  footprints first
   checkTrue("hg38human" %in% dbListTables(db))


} # test_privateGtfDatabaseOnWhovian
#----------------------------------------------------------------------------------------------------
# make sure that the wholeBrain tables can be read and queried by the trena user
test_privateWholeBrainFootprintDatabaseOnWhovian <- function()
{
   printf("--- test_privateWholeBrainFootprintDatabaseOnWhovian")

   db <- dbConnect(PostgreSQL(), user= "trena", password="trena", dbname="wholeBrain", host="whovian")

     # two tables at present: footprints and motifsGenes.  footprints first
   checkEquals(sort(dbListTables(db)), c("footprints", "motifsgenes"))

   checkTrue(dbExistsTable(db, "footprints"))
   query <- "select count(*) from footprints"
   tbl <- dbGetQuery(db, query)
   row.count <- tbl[1, "count"]
   checkTrue(row.count > 4e6)

   tbl <- dbGetQuery(db, "select * from footprints limit 2")
   checkEquals(nrow(tbl), 2)
   checkTrue(ncol(tbl) > 10)

     # now motifsGenes
   checkTrue(dbExistsTable(db, "motifsgenes"))
   query <- "select count(*) from motifsgenes"
   tbl <- dbGetQuery(db, query)
   row.count <- tbl[1, "count"]
   checkTrue(row.count > 9000)

   dbDisconnect(db)

} # test_privateWholeBrainFootprintDatabaseOnWhovian
#----------------------------------------------------------------------------------------------------
# make sure that the lymphoblast tables can be read and queried by the trena user
test_privateLymphoblastFootprintDatabaseOnWhovian <- function()
{
   printf("--- test_privateLymphoblastFootprintDatabaseOnWhovian")

   db <- dbConnect(PostgreSQL(), user= "trena", password="trena", dbname="lymphoblast", host="whovian")

     # two tables at present: footprints and motifsGenes.  footprints first
   checkEquals(sort(dbListTables(db)), c("footprints", "motifsgenes"))

   checkTrue(dbExistsTable(db, "footprints"))
   query <- "select count(*) from footprints"
   tbl <- dbGetQuery(db, query)
   row.count <- tbl[1, "count"]
   checkTrue(row.count > 4e6)

   tbl <- dbGetQuery(db, "select * from footprints limit 2")
   checkEquals(nrow(tbl), 2)
   checkTrue(ncol(tbl) > 10)

     # now motifsGenes
   checkTrue(dbExistsTable(db, "motifsgenes"))
   query <- "select count(*) from motifsgenes"
   tbl <- dbGetQuery(db, query)
   row.count <- tbl[1, "count"]
   checkTrue(row.count > 9000)

   dbDisconnect(db)

} # test_privateLymphoblastFootprintDatabaseOnWhovian
#----------------------------------------------------------------------------------------------------
# emulate the join which takes place behind the scenes in FootprintFinder
test_privateWholeBrainFootprintDatabaseOnWhovian_manualJoin <- function()
{
   printf("--- test_privateWholeBrainFootprintDatabaseOnWhovian_manualJoin")
   db <- dbConnect(PostgreSQL(), user= "trena", password="trena", dbname="wholeBrain", host="whovian")
   loc <- list(chr="chr5", start=88824783, end=88825650)

   query <- paste(c("select fp.chr, fp.mfpStart, fp.mfpEnd, fp.motifName, fp.pval, mg.motif, mg.tf",
                    "from footprints fp",
                    "inner join motifsgenes mg",
                     "on fp.motifName=mg.motif",
                     sprintf("where fp.chr = '%s' and  fp.mfpStart > %d and fp.mfpEnd < %d",
                             loc$chr, loc$start, loc$end)),
                    collapse=" ")

   tbl <- dbGetQuery(db, query)
   checkEquals(ncol(tbl), 7)
   checkTrue(nrow(tbl) > 100)

} # test_privateWholeBrainFootprintDatabaseOnWhovian_manualJoin
#----------------------------------------------------------------------------------------------------
# emulate the join which takes place behind the scenes in FootprintFinder
test_privateLymphoblastFootprintDatabaseOnWhovian_manualJoin <- function()
{
   printf("--- test_privateLymphoblastFootprintDatabaseOnWhovian_manualJoin")
   db <- dbConnect(PostgreSQL(), user= "trena", password="trena", dbname="lymphoblast", host="whovian")
   loc <- list(chr="chr5", start=88824783, end=88825650)

   query <- paste(c("select fp.chr, fp.mfpStart, fp.mfpEnd, fp.motifName, fp.pval, mg.motif, mg.tf",
                    "from footprints fp",
                    "inner join motifsgenes mg",
                     "on fp.motifName=mg.motif",
                     sprintf("where fp.chr = '%s' and  fp.mfpStart > %d and fp.mfpEnd < %d",
                             loc$chr, loc$start, loc$end)),
                    collapse=" ")

   tbl <- dbGetQuery(db, query)
   checkEquals(ncol(tbl), 7)
   checkTrue(nrow(tbl) > 30)

} # test_privateLymphoblastFootprintDatabaseOnWhovian_manualJoin
#----------------------------------------------------------------------------------------------------
test_getPromoterRegion <- function()
{
   printf("--- test_getPromoterRegion")

   genome.db.uri <- "postgres://whovian/hg38"
   project.db.uri <-  "postgres://whovian/lymphoblast"
   fp <- FootprintFinder(genome.db.uri, project.db.uri, quiet=TRUE)
   region <- getPromoterRegion(fp, "TREM2", 0, 0)

   checkEquals(region$chr, "chr6")
   checkEquals(region$start, 41163186)
   checkEquals(region$end,   41163186)

   region <- getPromoterRegion(fp, "TREM2", 20, 30)
   checkEquals(region$chr, "chr6")
   checkEquals(region$start, 41163156)  #
   checkEquals(region$end,   41163206)  # 20 bases upstream from the "end" TSS

     # now a plus strand gene
   region <- getPromoterRegion(fp, "SP1", 0, 0)
   checkEquals(region$chr, "chr12")
   checkEquals(region$start, 53380176)
   checkEquals(region$end,   53380176)

   region <- getPromoterRegion(fp, "SP1", 1000, 1)
   checkEquals(region$chr, "chr12")
   checkEquals(region$start, 53379176)
   checkEquals(region$end,   53380177)

} # test_getPromoterRegion
#----------------------------------------------------------------------------------------------------
test_getFootprints <- function()
{
   printf("--- test_getFootprints")

   genome.db.uri <- "postgres://whovian/hg38"
   project.db.uri <-  "postgres://whovian/lymphoblast"
   fp <- FootprintFinder(genome.db.uri, project.db.uri, quiet=TRUE)

   x <- getFootprints(fp, "MEF2C", 30000, 0); dim(x)
   dim(x) # 6 x 7   # should be MANY more!

   #tbl.fp <- getFootprints(fp, "MEF2C", 5000, 1000)
   #checkEquals(ncol(tbl.fp), 7)
   #checkTrue(nrow(tbl.fp) > 8)

  #    # an empty genomic span should return zero footprints
   #checkEquals(dim(getFootprints(fp, "MEF2C", 0, 0)), c(0, 7))

  #   # 3k up and downstream.  we expect more footpring upstream
   #                                     # some downstream
   #checkTrue(nrow(getFootprints(fp, "MEF2C", 3000, 0)) > 50)
   #checkTrue(nrow(getFootprints(fp, "MEF2C", 0, 3000)) < 20)

  #   # try the other gene we provided to constructor
   #checkTrue(nrow(getFootprints(fp, "ZNF454", 500, 500)) > 100)

  #   # now try a gene not provided to the constructor
   #checkException(getFootprints(fp, "IL1A", 500, 500), silent=TRUE)

} # test_getFootprints
#----------------------------------------------------------------------------------------------------
test_constructWithAllChromosome5genes <- function()
{
   printf("--- test_constructWithAllChromosome5genes")

   ensembl.hg38 <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
   tbl <- getBM(attributes=c("entrezgene", "chromosome_name", "hgnc_symbol"), filters="chromosome_name",
                values="5", mart=ensembl.hg38)

   genes.chr5 <- tbl$hgnc_symbol
   missing <- which(genes.chr5 == "")
   if(length(missing) > 0)
      genes.chr5 <- genes.chr5[-missing]  # 1515

   #db.uri <- "sqlite:~/github/snpFoot/inst/misc/sqlite.explorations/fpTf.sqlite"
   db.uri <- sprintf("sqlite:%s", system.file(package="TReNA.brain", "extdata", "fpTf.sqlite"))
   fp <- FootprintFinder(ensembl.hg38, db.uri, genes.chr5)
   goi <- genes.chr5[length(genes.chr5)]
   tbl.fp <- getFootprints(fp, goi, 1000, 1000)
   checkTrue(nrow(tbl.fp) > 200)

} # test_constructWithAllChromosome5genes
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
   loc <- getPromoterRegion(fp, goi, 250, 0)
   tbl <- getFootprints(fp, goi, 250, 0)
   checkTrue(all(tbl$mfpStart >= loc$start))
   checkTrue(all(tbl$mfpStart <= loc$end))
   checkTrue(all(tbl$mfpEnd >= loc$start))
   checkTrue(all(tbl$mfpEnd <= loc$end))

} # test_getFootprintsForEnsemblGenes
#----------------------------------------------------------------------------------------------------
if(!interactive()) runTests()
