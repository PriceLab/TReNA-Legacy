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

    mtx <- matrix(1:9,nrow=3)   
    rownames(mtx) <- c("gene1","gene2","gene3")    
    solver <- PearsonSolver(mtx,targetGene = "gene1",                          
                          candidateRegulators = c("gene2","gene3"))    
   
    checkEquals(class(solver)[1], "PearsonSolver")    
    checkTrue(all(c("PearsonSolver", "Solver") %in% is(solver)))
    
} # test_PearsonSolverConstructor
#----------------------------------------------------------------------------------------------------
test_ampAD.mef2c.154tfs.278samples.pearson <- function()
{
   printf("--- test_ampAD.mef2c.154tfs.278samples.pearson")

   # Load matrix and transform via arcsinh
   load(system.file(package="TReNA", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
   
   mtx.asinh <- asinh(mtx.sub)
   target.gene <- "MEF2C"
   tfs <- setdiff(rownames(mtx.asinh), "MEF2C")
   #print(fivenum(mtx.asinh)  # [1] 0.000000 1.327453 3.208193 4.460219 7.628290)   

   pearson.solver <- PearsonSolver(mtx.asinh, target.gene, tfs)
   tbl <- run(pearson.solver)

   # Check for empirical values
   checkTrue(nrow(subset(tbl, abs(coefficient) > 0.8)) == 8)

} # test_ampAD.mef2c.154tfs.278samples.pearson
#----------------------------------------------------------------------------------------------------
if(!interactive()) runTests()
