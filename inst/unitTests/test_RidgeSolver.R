library(TReNA)
library(RUnit)
#----------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_RidgeSolverConstructor()
   test_ampAD.mef2c.154tfs.278samples.ridge()

} # runTests
#----------------------------------------------------------------------------------------------------
test_RidgeSolverConstructor <- function()
{
    printf("--- test_RidgeSolverConstructor")

    # Construct the RidgeSolver and check that it's correct
    #solver <- RidgeSolver()
    solver <- TReNA(matrix(), solver = "ridge")
    checkEquals(getSolverName(solver), "RidgeSolver")
    #checkTrue(all(c("RidgeSolver", "Solver") %in% is(solver)))
}

# test_RidgeSolverConstructor
#----------------------------------------------------------------------------------------------------
test_ampAD.mef2c.154tfs.278samples.ridge <- function()
{
   printf("--- test_ampAD.mef2c.154tfs.278samples.ridge")

   # Load matrix and transform via arcsinh
   load(system.file(package="TReNA", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
   target.gene <- "MEF2C"
   mtx.asinh <- asinh(mtx.sub)
   #print(fivenum(mtx.asinh)  # [1] 0.000000 1.327453 3.208193 4.460219 7.628290)
   
   trena <- TReNA(mtx.assay=mtx.asinh, solver="ridge", quiet=FALSE)
   tfs <- setdiff(rownames(mtx.asinh), "MEF2C")
   tbl <- solve(trena, target.gene, tfs)

   # Check for empirical values
   checkTrue(min(tbl$beta) > -0.15)
   checkTrue(max(tbl$beta) < 0.15)
   checkTrue(c("FOXP1") %in% rownames(subset(tbl, abs(beta) > 0.08)))

} # test_ampAD.mef2c.154tfs.278samples.ridge
#----------------------------------------------------------------------------------------------------
if(!interactive()) runTests()
