library(TReNA)
library(RUnit)
library(biomaRt)
library(TReNA.brain)  # sqlite footprint database and expession data
library(TReNA.lymphoblast)
#----------------------------------------------------------------------------------------------------
if(!exists("ensembl.hg38"))
    ensembl.hg38 <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

#----------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_constructor()
   test_getPromoterRegion()
   test_getFootprints()
   test_constructWithAllChromosome5genes()
   test_getFootprintsForEnsemblGenes()

} # runTests
#----------------------------------------------------------------------------------------------------
test_constructor <- function()
{
   printf("--- test_constructor")
   db.uri <- sprintf("sqlite:%s", system.file(package="TReNA.brain", "extdata", "fpTf.sqlite"))
   genes <- c("MEF2C", "ZNF454")
   fp <- FootprintFinder(ensembl.hg38, db.uri, genes)
   checkEquals(getGenes(fp), genes)
   invisible(fp)

} # test_constructor
#----------------------------------------------------------------------------------------------------
test_getPromoterRegion <- function()
{
   printf("--- test_getPromoterRegion")
   #db.uri <- "sqlite:~/github/snpFoot/inst/misc/sqlite.explorations/fpTf.sqlite"
   db.uri <- sprintf("sqlite:%s", system.file(package="TReNA.brain", "extdata", "fpTf.sqlite"))
   genes <- c("SP1", "TREM2")
   fp <- FootprintFinder(ensembl.hg38, db.uri, genes)

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
   db.uri <- sprintf("sqlite:%s", system.file(package="TReNA.brain", "extdata", "fpTf.sqlite"))
   genes <- c("MEF2C", "ZNF454")
   fp <- FootprintFinder(ensembl.hg38, db.uri, genes)
   tbl.fp <- getFootprints(fp, "MEF2C", 1000, 1000)
   checkEquals(ncol(tbl.fp), 7)
   checkTrue(nrow(tbl.fp) > 8)

      # an empty genomic span should return zero footprints
   checkEquals(dim(getFootprints(fp, "MEF2C", 0, 0)), c(0, 7))

     # 3k up and downstream.  we expect more footpring upstream
                                        # some downstream
   checkTrue(nrow(getFootprints(fp, "MEF2C", 3000, 0)) > 50)
   checkTrue(nrow(getFootprints(fp, "MEF2C", 0, 3000)) < 20)

     # try the other gene we provided to constructor
   checkTrue(nrow(getFootprints(fp, "ZNF454", 500, 500)) > 100)

     # now try a gene not provided to the constructor
   checkException(getFootprints(fp, "IL1A", 500, 500), silent=TRUE)

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
