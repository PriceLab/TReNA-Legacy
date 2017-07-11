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

    mtx <- matrix(1:9,nrow=3)   
    rownames(mtx) <- c("gene1","gene2","gene3")    
    solver <- SpearmanSolver(mtx,targetGene = "gene1",                          
                          candidateRegulators = c("gene2","gene3"))    
   
    checkEquals(class(solver)[1], "SpearmanSolver")    
    checkTrue(all(c("SpearmanSolver", "Solver") %in% is(solver)))

} # test_SpearmanSolverConstructor
#----------------------------------------------------------------------------------------------------
test_ampAD.mef2c.154tfs.278samples.spearman <- function()
{
   printf("--- test_ampAD.mef2c.154tfs.278samples.spearman")

   # Load matrix and transform via arcsinh
   load(system.file(package="TReNA", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
   
   mtx.asinh <- asinh(mtx.sub)
   target.gene <- "MEF2C"
   tfs <- setdiff(rownames(mtx.asinh), "MEF2C")
   #print(fivenum(mtx.asinh)  # [1] 0.000000 1.327453 3.208193 4.460219 7.628290)   
   
   spearman.solver <- SpearmanSolver(mtx.asinh, target.gene, tfs)
   tbl <- run(spearman.solver)

   # Check for empirical values
   checkTrue(nrow(subset(tbl, abs(coefficient) > 0.8)) == 9)

} # test_ampAD.mef2c.154tfs.278samples.spearman
#----------------------------------------------------------------------------------------------------
if(!interactive()) runTests()
