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

} # runTests
#----------------------------------------------------------------------------------------------------
test_constructor <- function()
{
   printf("--- test_constructor")
   db.uri <- "sqlite:~/github/snpFoot/inst/misc/sqlite.explorations/fpTf.sqlite"
   genes <- c("MEF2C", "ZNF454")
   fp <- Footprints(ensembl.hg38, db.uri, genes)
   checkEquals(getGenes(fp), genes)
   invisible(fp)

} # test_constructor
#----------------------------------------------------------------------------------------------------
test_getPromoterRegion <- function()
{
   print("--- test_getPromoterRegion")
   db.uri <- "sqlite:~/github/snpFoot/inst/misc/sqlite.explorations/fpTf.sqlite"
   genes <- c("MEF2C", "ZNF454")
   fp <- Footprints(ensembl.hg38, db.uri, genes)

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
   fp <- Footprints(ensembl.hg38, db.uri, genes)
   tbl.fp <- getFootprints(fp, "MEF2C", 1000, 1000)
   checkEquals(ncol(tbl.fp), 7)
   checkTrue(nrow(tbl.fp) > 8)

} # test_getFootprints
#----------------------------------------------------------------------------------------------------
