library(TReNA)
library(RUnit)
#----------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_SqrtLassoSolverConstructor()
   test_ampAD.mef2c.154tfs.278samples.sqrtlasso()
   
} # runTests
#----------------------------------------------------------------------------------------------------
test_SqrtLassoSolverConstructor <- function()
{
    printf("--- test_SqrtLassoSolverConstructor")

    # Construct the SqrtLassoSolver and check that it's correct
    #solver <- SqrtLassoSolver()
    solver <- TReNA(matrix(), solver = "sqrtlasso")
    checkEquals(getSolverName(solver), "SqrtLassoSolver")
    #checkTrue(all(c("SqrtLassoSolver", "Solver") %in% is(solver)))
}

# test_SqrtLassoSolverConstructor   
#----------------------------------------------------------------------------------------------------    
test_ampAD.mef2c.154tfs.278samples.sqrtlasso <- function()
{
   printf("--- test_ampAD.mef2c.154tfs.278samples.sqrtlasso")

   # Load matrix and transform via arcsinh
   load(system.file(package="TReNA", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
   target.gene <- "MEF2C"
   mtx.asinh <- asinh(mtx.sub)
   #print(fivenum(mtx.asinh)  # [1] 0.000000 1.327453 3.208193 4.460219 7.628290)
   
   trena <- TReNA(mtx.assay=mtx.asinh, solver="sqrtlasso", quiet=FALSE)
   tfs <- setdiff(rownames(mtx.asinh), "MEF2C")
   tbl <- solve(trena, target.gene, tfs, extraArgs = list("num.cores" = parallel::detectCores()))

   # Check for empirical values
   tbl <- tbl[order(abs(tbl$beta), decreasing=TRUE),, drop = FALSE]
   expected.genes <- sort(c("SATB2","STAT4","HLF","TSHZ3", "FOXP2"))
   actual.genes <- sort(rownames(subset(tbl, abs(beta) > 0.1)))
   printf("Top 5 genes: %s", paste(rownames(tbl[1:5,]),collapse=","))
   checkEquals(expected.genes,actual.genes)
   
} # test_ampAD.mef2c.154tfs.278samples.sqrtlasso
#----------------------------------------------------------------------------------------------------    
if(!interactive()) runTests()
