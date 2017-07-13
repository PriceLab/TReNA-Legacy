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
    mtx <- matrix(1:9,nrow=3)
    rownames(mtx) <- c("gene1","gene2","gene3")
    solver <- EnsembleSolver(mtx,targetGene = "gene1",
                             candidateRegulators = c("gene2","gene3"))
    checkEquals(class(solver)[1], "EnsembleSolver")
    checkTrue(all(c("EnsembleSolver", "Solver") %in% is(solver)))
}

# test_EnsembleSolverConstructor
#----------------------------------------------------------------------------------------------------
test_ampAD.mef2c.154tfs.278samples.ensemble <- function()
{
   printf("--- test_ampAD.mef2c.154tfs.278samples.ensemble")

   set.seed(122113)
   # Load matrix and transform via arcsinh
   load(system.file(package="TReNA", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
    mtx.asinh <- asinh(mtx.sub)
   #print(fivenum(mtx.asinh)  # [1] 0.000000 1.327453 3.208193 4.460219 7.628290)

   target.gene <- "MEF2C"
   tfs <- setdiff(rownames(mtx.asinh), "MEF2C")
   solver <- EnsembleSolver(mtx.asinh,target.gene,tfs)
   tbl <- run(solver)

   # Check for empirical values
   checkTrue(min(tbl$pcaMax) > 0.6)
   checkTrue(max(tbl$pcaMax) < 2.5)
   checkTrue(min(tbl$concordance) > 0.3)
   checkTrue(max(tbl$concordance) < 0.55)
   checkTrue(c("HLF") %in% tbl$gene)

} # test_ampAD.mef2c.154tfs.278samples.ensemble
#----------------------------------------------------------------------------------------------------
test_selectedSolversOnly <- function()
{
   printf("--- test_selectedSolversOnly")

   set.seed(122113)
   # Load matrix and transform via arcsinh
   load(system.file(package="TReNA", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))   
   mtx.asinh <- asinh(mtx.sub)
   #print(fivenum(mtx.asinh)  # [1] 0.000000 1.327453 3.208193 4.460219 7.628290)

   target.gene <- "MEF2C"
   tfs <- setdiff(rownames(mtx.asinh), "MEF2C")
   solvers <- c("lasso", "ridge", "lassopv", "pearson", "spearman") # "sqrtlasso",

   solver <- EnsembleSolver(mtx.asinh,target.gene,tfs,solverNames=solvers)
   tbl <- run(solver)

   # Check for empirical values
   checkTrue(min(tbl$pcaMax) > 0.8)
   checkTrue(max(tbl$pcaMax) < 1.9)
   checkTrue(min(tbl$concordance) > 0.35)
   checkTrue(max(tbl$concordance) < 0.55)
   checkTrue(c("HLF") %in% tbl$gene)

} # test_selectedSolversOnly
#----------------------------------------------------------------------------------------------------
if(!interactive()) runTests()
