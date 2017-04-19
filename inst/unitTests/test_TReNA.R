library(TReNA)
library(RUnit)
#----------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
    test_emptyConstructor()
    test_getSolver()
   
} # runTests
#----------------------------------------------------------------------------------------------------
test_emptyConstructor <- function()
{
   printf("--- test_emptyConstructor")
   trena <- TReNA()
   checkEquals(is(trena), "TReNA")

   #mtx <- getAssayData(trena)
   #checkTrue("matrix" %in% is(mtx))
   #checkEquals(dim(mtx), c(1,1))   # the default (and implicilty empty) matrix
   #checkTrue(is.na(mtx[1,1]))

   #mtx.priors <- getPriors(trena)
   #checkTrue("matrix" %in% is(mtx.priors))
   #checkEquals(dim(mtx.priors), c(1,1))

} # test_emptyConstructor
#----------------------------------------------------------------------------------------------------
test_getSolver <- function()
{
    printf("--- test_getSolver")
    trena <- TReNA(matrix(), solver = "lasso")
    checkTrue("LassoSolver" %in% class(getSolverObject(trena)))
    checkEquals("LassoSolver", getSolverName(trena))

} # test_getSolver
#----------------------------------------------------------------------------------------------------
if(!interactive()) runTests()
