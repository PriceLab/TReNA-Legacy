library(TReNA)
library(RUnit)
#----------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_EnsembleSolverConstructor()
   test_ampAD.mef2c.154tfs.278samples.ensemble()
   test_selectedSolversOnly()

} # runTests
#----------------------------------------------------------------------------------------------------
test_EnsembleSolverConstructor <- function()
{
    printf("--- test_EnsembleSolverConstructor")

    # Construct the EnsembleSolver and check that it's correct

    #solver <- EnsembleSolver()
    solver <- TReNA(matrix(), solver = "ensemble")
    checkEquals(getSolverName(solver), "EnsembleSolver")
    #checkTrue(all(c("EnsembleSolver", "Solver") %in% is(solver)))
}

# test_EnsembleSolverConstructor
#----------------------------------------------------------------------------------------------------
test_ampAD.mef2c.154tfs.278samples.ensemble <- function()
{
   printf("--- test_ampAD.mef2c.154tfs.278samples.ensemble")

   set.seed(122113)
   # Load matrix and transform via arcsinh
   load(system.file(package="TReNA", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
   target.gene <- "MEF2C"
   mtx.asinh <- asinh(mtx.sub)
   #print(fivenum(mtx.asinh)  # [1] 0.000000 1.327453 3.208193 4.460219 7.628290)

   trena <- TReNA(mtx.assay=mtx.asinh, solver="ensemble", quiet=FALSE)
   tfs <- setdiff(rownames(mtx.asinh), "MEF2C")
   tbl <- solve(trena, target.gene, tfs)

   # Check for empirical values
   checkTrue(min(tbl$pcaMax) > 1.2)
   checkTrue(max(tbl$pcaMax) < 5.2)
   checkTrue(min(tbl$concordance) > 0.4)
   checkTrue(max(tbl$concordance) < 0.7)
   checkTrue(c("HLF") %in% tbl$gene)

} # test_ampAD.mef2c.154tfs.278samples.ensemble
#----------------------------------------------------------------------------------------------------
test_selectedSolversOnly <- function()
{
   printf("--- test_selectedSolversOnly")

   set.seed(122113)
   # Load matrix and transform via arcsinh
   load(system.file(package="TReNA", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
   target.gene <- "MEF2C"
   mtx.asinh <- asinh(mtx.sub)
   #print(fivenum(mtx.asinh)  # [1] 0.000000 1.327453 3.208193 4.460219 7.628290)

      # solvers <- c("lasso", "ridge", "randomForest",  "lassopv", "pearson", "spearman") # "sqrtlasso",

   solver <- "ensemble:lasso;ridge;randomForest"
   trena <- TReNA(mtx.assay=mtx.asinh, solver=solver, quiet=FALSE)
   tfs <- setdiff(rownames(mtx.asinh), "MEF2C")
   tbl <- solve(trena, target.gene, tfs)

   # Check for empirical values
   checkTrue(min(tbl$extr) > 1.2)
   checkTrue(max(tbl$extr) < 5.0)
   checkTrue(c("HLF") %in% tbl$gene)

} # test_selectedSolversOnly
#----------------------------------------------------------------------------------------------------
if(!interactive()) runTests()
