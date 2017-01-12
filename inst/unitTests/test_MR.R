library(TReNA)
library(RUnit)
                                        #---------------
printf <- function(...) print(noquote(sprintf(...)))
                                        #--------------
runTests <- function()
{
    test_PearsonSolverConstructor()
    test_fitDummyData.Pearson()
} # runTests
                                        #---------------
test_PearsonSolverConstructor <- function()
{
    printf("--- test_PearsonSolverConstructor")
    solver <- PearsonSolver()
    checkEquals(getSolverName(solver), "PearsonSolver")
    checkTrue(all(c("PearsonSolver", "Solver") %in% is(solver)))

} # test_PearsonSolverConstructor
                                        #----------------
test_fitDummyData.Pearson <- function()
{
    printf("--- test_fitDummyData.Pearson")
                                        # Grab the dummy data function and create the data
    # Hard code path for now; fix this later
    source("/local/mrichard/github/TReNA/inst/unitTests/test_TReNA.R")
    x <- test_developAndFitDummyTestData(quiet=TRUE)
    mtx <- x$assay
    tfs <- x$tf.genes
    target.gene <- x$correlated.target

    # Setup and solve using Pearson solver
    trena <- TReNA(mtx.assay=mtx,solver="pearson",quiet=FALSE)    
    coeffs <- solve(trena,target.gene,tfs)    

                                        # Select only genes and coefficients above 0.2
    sub.coeffs <- coeffs[abs(coeffs) > 0.2]

    checkTrue(length(sub.coeffs) > 5)
    
} #test_fitDummyData.Pearson            

                                        #--------------------

                                        # Run tests when not interactive
if(!interactive()) runTests()
                                      
