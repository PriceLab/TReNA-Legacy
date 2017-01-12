library(TReNA)
library(RUnit)
                                        #---------------
printf <- function(...) print(noquote(sprintf(...)))
                                        #--------------
runTests <- function()
{
    test_PearsonSolverConstructor()

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
                                        # Run tests when not interactive
if(!interactive()) runTests()
