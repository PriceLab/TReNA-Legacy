library(TReNA)
library(RUnit)
#----------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
if(!exists("mtx")){
   load(system.file(package="TReNA", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
   mtx <- asinh(mtx.sub)
   }

candidate.tfs <- c("ATF7", "NR3C2", "MAFB", "PRRX1", "E2F8", "XBP1"); # from gtex skin data(!)
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_RandomForestSolverFewCandidates()
   test_ampAD.mef2c.154tfs.278samples.randomForest()

} # runTests
#----------------------------------------------------------------------------------------------------
test_RandomForestSolverFewCandidates <- function()
{
   printf("--- test_RandomForestSolverFewCandidates")
   solver <- RandomForestSolver(mtx, targetGene="MEF2C", candidateRegulators=candidate.tfs)
   x <- run(solver)
   checkTrue(all(c("edges", "r2") %in% names(x)))
   tbl <- x$edges
   checkEquals(dim(tbl), c(3, 2))
   checkEquals(colnames(tbl), c("IncNodePurity", "gene.cor"))
   checkEquals(rownames(tbl), c("ATF7", "NR3C2", "PRRX1"))

} # test_RandomForestSolverFewCandidates
#----------------------------------------------------------------------------------------------------
test_ampAD.mef2c.154tfs.278samples.randomForest <- function()
{
   printf("--- test_ampAD.mef2c.154tfs.278samples.randomForest")

   targetGene <- "MEF2C"
   candidate.tfs <- setdiff(rownames(mtx), targetGene)
   solver <- RandomForestSolver(mtx, targetGene=targetGene, candidateRegulators=candidate.tfs)
   x <- run(solver)
   tbl <- x$edges
      # check just the highest scores
   tbl.10 <- subset(tbl, IncNodePurity > 10)
   checkEquals(rownames(tbl.10), c("HLF", "STAT4", "SATB2", "SATB1"))

} # test_ampAD.mef2c.154tfs.278samples.randomForest
#----------------------------------------------------------------------------------------------------
if(!interactive()) runTests()
