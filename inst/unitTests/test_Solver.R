library(TReNA)
library(RUnit)
#----------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
    test_getAssayData()
    test_getTarget()
    test_getRegulators()
    test_eliminateSelfTFs()    
    test_MatrixWarnings()    
   
} # runTests
#----------------------------------------------------------------------------------------------------
test_getAssayData <- function()
{
    printf("--- test_getAssayData")
    mtx <- matrix(1:9, nrow = 3)
    rownames(mtx) <- c("gene1","gene2","gene3")
    solver <- Solver(mtx, "gene1", c("gene2","gene3"))    
    checkEquals(class(getAssayData(solver)), "matrix")
    checkEquals(mtx, getAssayData(solver))
    
    } # test_getAssayData
#----------------------------------------------------------------------------------------------------
test_getTarget <- function()
{
    printf("--- test_getTarget")
    mtx <- matrix(1:9, nrow = 3)
    rownames(mtx) <- c("gene1","gene2","gene3")
    solver <- Solver(mtx, "gene1", c("gene2","gene3"))    
    checkEquals(class(getTarget(solver)), "character")
    checkEquals("gene1", getTarget(solver))
    
    } # test_getTarget
#----------------------------------------------------------------------------------------------------
test_getRegulators <- function()
{
    printf("--- test_getRegulators")
    mtx <- matrix(1:9, nrow = 3)
    rownames(mtx) <- c("gene1","gene2","gene3")
    solver <- Solver(mtx, "gene1", c("gene2","gene3"))    
    checkEquals(length(getRegulators(solver)), 2)
    checkEquals(c("gene2","gene3"), getRegulators(solver))
    
    } # test_getRegulators
#----------------------------------------------------------------------------------------------------
test_eliminateSelfTFs <- function()
{
   printf("--- test_eliminateSelfTFs")

   set.seed(10045)
   load(system.file(package="TReNA", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
   
   mtx.asinh <- asinh(mtx.sub)
   target.gene <- "MEF2C"
   tfs <- rownames(mtx.asinh)
   
   solver <- PearsonSolver(mtx.asinh, target.gene, tfs)
   
   checkTrue(target.gene %in% tfs)         # our test case
   tbl.betas <- run(solver)
   checkTrue(!target.gene %in% rownames(tbl.betas))

   solver2 <- LassoSolver(mtx.asinh, target.gene, tfs)
   tbl.betas2 <- run(solver2)
   checkTrue(!target.gene %in% rownames(tbl.betas2))

} # test_eliminateSelfTFs
#----------------------------------------------------------------------------------------------------
test_MatrixWarnings <- function()
{
    printf("--- test_MatrixWarnings")

    # Change warnings to errors
    options(warn = 2)

    # Check that a skewed matrix produces an error
    test.mtx <- matrix(1:10000, nrow = 100)
    test.mtx[1,1] <- 1e7
    rownames(test.mtx) <- paste0("gene",1:100)
    target.gene <- "gene1"
    tfs <- setdiff(rownames(test.mtx),target.gene)
    checkException(solver(test.mtx,target.gene,tfs), silent = TRUE)

    # Check that a matrix with a row of 0's produces an error for most solvers
    test.mtx[1,] <- 0
    checkException(BayesSpikeSolver(test.mtx, target.gene, tfs), silent = TRUE)
    checkException(LassoPVSolver(test.mtx, target.gene, tfs), silent = TRUE)
    checkException(SqrtLassoSolver(test.mtx, target.gene, tfs), silent = TRUE)
    checkException(RandomForestSolver(test.mtx, target.gene, tfs), silent = TRUE)
    checkException(PearsonSolver(test.mtx, target.gene, tfs), silent = TRUE)
    checkException(SpearmanSolver(test.mtx, target.gene, tfs), silent = TRUE)
    checkException(EnsembleSolver(test.mtx, target.gene, tfs), silent = TRUE)

    # Check that a target gene with low expression causes a warning for a solver
    test.mtx[1,] <- 0.1
    solver <- EnsembleSolver(test.mtx, target.gene, tfs)
    checkException(run(solver), silent = TRUE)    
    
    # Change warnings back to warnings
    options(warn = 1)

} #test_MatrixWarnings
#----------------------------------------------------------------------------------------------------
if(!interactive()) runTests()
