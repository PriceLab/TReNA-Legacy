#----------------------------------------------------------------------------------------------------
# Unit Tests for Naive Solver
library(TReNA)
library(RUnit)
#----------------------------------------------------------------------------------------------------
# Run All Tests
runTests <- function() {
  test_NaiveSolverConstructor()
  test_ampAD.mef2c.154tfs.278samples.naive()
}
#----------------------------------------------------------------------------------------------------
# Constructor Test
test_NaiveSolverConstructor <- function() {
  printf("--- test_NaiveSolverConstructor")
  
  # Construct solver & get name
  solver <- TReNA(matrix(), solver = "naive")
  checkEquals(getSolverName(solver), "NaiveSolver")
}
#----------------------------------------------------------------------------------------------------
# MEF2C Data Test
test_ampAD.mef2c.154tfs.278samples.naive <- function() {
  printf("--- test_ampAD.mef2c.154tfs.278samples.naive")
  
  # Load matrix and transform via arcsinh
  load(system.file(package="TReNA", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
  target.gene <- "MEF2C"
  mtx.asinh <- asinh(mtx.sub)
  
  trena <- TReNA(mtx.assay=mtx.asinh, solver="naive", quiet=FALSE)
  tfs <- setdiff(rownames(mtx.asinh), "MEF2C")
  tbl <- solve(trena, target.gene, tfs)
  
  # Checks
  checkTrue(min(tbl$beta) > -0.3)
  checkTrue(max(tbl$beta) < 0.3)
  checkTrue(min(tbl$p.value) > 1e-14)
  checkTrue(max(tbl$p.value) < 1e-01)
}
