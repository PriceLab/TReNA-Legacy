library(TReNA)
library(RUnit)
library(biomaRt)
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

} # runTests
#----------------------------------------------------------------------------------------------------
test_constructor <- function()
{
   printf("--- test_constructor")
   db.uri <- "sqlite:~/github/snpFoot/inst/misc/sqlite.explorations/fpTf.sqlite"
   genes <- c("MEF2C", "ZNF454")
   fp <- FootprintFinder(ensembl.hg38, db.uri, genes)
   checkEquals(getGenes(fp), genes)
   invisible(fp)

} # test_constructor
#----------------------------------------------------------------------------------------------------
test_getPromoterRegion <- function()
{
   printf("--- test_getPromoterRegion")
   db.uri <- "sqlite:~/github/snpFoot/inst/misc/sqlite.explorations/fpTf.sqlite"
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
   db.uri <- "sqlite:~/github/snpFoot/inst/misc/sqlite.explorations/fpTf.sqlite"
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

   db.uri <- "sqlite:~/github/snpFoot/inst/misc/sqlite.explorations/fpTf.sqlite"
   fp <- FootprintFinder(ensembl.hg38, db.uri, genes.chr5)
   goi <- genes.chr5[length(genes.chr5)]
   tbl.fp <- getFootprints(fp, goi, 1000, 1000)
   checkTrue(nrow(tbl.fp) > 200)

} # test_constructWithAllChromosome5genes
#----------------------------------------------------------------------------------------------------
