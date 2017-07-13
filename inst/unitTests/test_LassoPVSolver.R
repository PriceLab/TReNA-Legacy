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

    mtx <- matrix(1:9,nrow=3)   
    rownames(mtx) <- c("gene1","gene2","gene3")    
    solver <- LassoPVSolver(mtx,targetGene = "gene1",                          
                          candidateRegulators = c("gene2","gene3"))    
   
    checkEquals(class(solver)[1], "LassoPVSolver")    
    checkTrue(all(c("LassoPVSolver", "Solver") %in% is(solver)))
}

# test_LassoPVSolverConstructor
#----------------------------------------------------------------------------------------------------
test_ampAD.mef2c.154tfs.278samples.lassopv <- function()
{
   printf("--- test_ampAD.mef2c.154tfs.278samples.lassopv")

   # Load matrix and transform via arcsinh
   load(system.file(package="TReNA", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
   mtx.asinh <- asinh(mtx.sub)
   target.gene <- "MEF2C"
   tfs <- setdiff(rownames(mtx.asinh), "MEF2C")   
   #print(fivenum(mtx.asinh)  # [1] 0.000000 1.327453 3.208193 4.460219 7.628290)

   lassopv.solver <- LassoPVSolver(mtx.asinh, target.gene, tfs)
   tbl <- run(lassopv.solver)

   # Check for significant P-values; make sure they match the empirical value
   sig.genes <- tbl$p.values[tbl$p.values < 0.01]
   checkEquals(length(sig.genes),30)
   checkTrue(max(-log10(sig.genes)) > 54)
   
} # test_ampAD.mef2c.154tfs.278samples.lassopv
#----------------------------------------------------------------------------------------------------
if(!interactive()) runTests()
