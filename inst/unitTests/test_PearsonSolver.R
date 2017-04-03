library(TReNA)
library(RUnit)
#----------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_PearsonSolverConstructor()
   test_ampAD.mef2c.154tfs.278samples.pearson()
   
} # runTests
#----------------------------------------------------------------------------------------------------
test_PearsonSolverConstructor <- function()
{
    printf("--- test_PearsonSolverConstructor")
    #solver <- PearsonSolver()
    solver <- TReNA(matrix(), solver = "pearson")
    checkEquals(getSolverName(solver), "PearsonSolver")
    #checkTrue(all(c("PearsonSolver", "Solver") %in% is(solver)))

} # test_PearsonSolverConstructor
#----------------------------------------------------------------------------------------------------
test_ampAD.mef2c.154tfs.278samples.pearson <- function()
{
   printf("--- test_ampAD.mef2c.154tfs.278samples.pearson")

   # Load matrix and transform via arcsinh
   load(system.file(package="TReNA", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
   target.gene <- "MEF2C"
   mtx.asinh <- asinh(mtx.sub)
   #print(fivenum(mtx.asinh)  # [1] 0.000000 1.327453 3.208193 4.460219 7.628290)
   
   trena <- TReNA(mtx.assay=mtx.asinh, solver="pearson", quiet=FALSE)
   tfs <- setdiff(rownames(mtx.asinh), "MEF2C")
   tbl <- solve(trena, target.gene, tfs)

   # Check for empirical values
   checkTrue(nrow(subset(tbl, abs(coefficient) > 0.8)) > 7)

} # test_ampAD.mef2c.154tfs.278samples.pearson
#----------------------------------------------------------------------------------------------------
if(!interactive()) runTests()
