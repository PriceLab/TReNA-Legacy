library(TReNA)
library(RUnit)
#----------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_SpearmanSolverConstructor()
   test_ampAD.mef2c.154tfs.278samples.spearman()
   
} # runTests
#----------------------------------------------------------------------------------------------------
test_SpearmanSolverConstructor <- function()
{
    printf("--- test_SpearmanSolverConstructor")
    solver <- SpearmanSolver()
    checkEquals(getSolverName(solver), "SpearmanSolver")
    checkTrue(all(c("SpearmanSolver", "Solver") %in% is(solver)))

} # test_SpearmanSolverConstructor
#----------------------------------------------------------------------------------------------------
test_ampAD.mef2c.154tfs.278samples.spearman <- function()
{
   printf("--- test_ampAD.mef2c.154tfs.278samples.spearman")

   # Load matrix and transform via arcsinh
   load(system.file(package="TReNA", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
   target.gene <- "MEF2C"
   mtx.asinh <- asinh(mtx.sub)
   #print(fivenum(mtx.asinh)  # [1] 0.000000 1.327453 3.208193 4.460219 7.628290)
   
   trena <- TReNA(mtx.assay=mtx.asinh, solver="spearman", quiet=FALSE)
   tfs <- setdiff(rownames(mtx.asinh), "MEF2C")
   tbl <- solve(trena, target.gene, tfs)

   # Check for empirical values
   checkTrue(nrow(subset(tbl, abs(coefficient) > 0.8)) > 7)

} # test_ampAD.mef2c.154tfs.278samples.spearman
#----------------------------------------------------------------------------------------------------
if(!interactive()) runTests()
