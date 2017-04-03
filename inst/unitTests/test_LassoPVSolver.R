library(TReNA)
library(RUnit)
#----------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_LassoPVSolverConstructor()
   test_ampAD.mef2c.154tfs.278samples.lassopv()
   
} # runTests
#----------------------------------------------------------------------------------------------------    
test_LassoPVSolverConstructor <- function()
{
    printf("--- test_LassoPvSolverConstructor")

    # Construct the SqrtLassoSolver and check that it's correct

    #solver <- LassoPVSolver()
    solver <- TReNA(matrix(), solve = "lassopv")
    checkEquals(getSolverName(solver), "LassoPVSolver")
    #checkTrue(all(c("LassoPVSolver", "Solver") %in% is(solver)))
}

# test_LassoPVSolverConstructor
#----------------------------------------------------------------------------------------------------
test_ampAD.mef2c.154tfs.278samples.lassopv <- function()
{
   printf("--- test_ampAD.mef2c.154tfs.278samples.lassopv")

   # Load matrix and transform via arcsinh
   load(system.file(package="TReNA", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
   target.gene <- "MEF2C"
   mtx.asinh <- asinh(mtx.sub)
   #print(fivenum(mtx.asinh)  # [1] 0.000000 1.327453 3.208193 4.460219 7.628290)
   
   trena <- TReNA(mtx.assay=mtx.asinh, solver="lassopv", quiet=FALSE)
   tfs <- setdiff(rownames(mtx.asinh), "MEF2C")
   tbl <- solve(trena, target.gene, tfs)

   # Check for significant P-values; make sure they match the empirical value
   sig.genes <- tbl$p.values[tbl$p.values < 0.01]
   checkEquals(length(sig.genes),30)
   
} # test_ampAD.mef2c.154tfs.278samples.lassopv
#----------------------------------------------------------------------------------------------------
if(!interactive()) runTests()
