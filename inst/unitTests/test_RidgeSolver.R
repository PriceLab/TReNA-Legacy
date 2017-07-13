library(TReNA)
library(RUnit)
#----------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_RidgeSolverConstructor()
   test_ampAD.mef2c.154tfs.278samples.ridge()
   test_alpha.ridge()
   test_lambda.ridge()
   test_keep.metrics.ridge()

} # runTests
#----------------------------------------------------------------------------------------------------
test_RidgeSolverConstructor <- function()
{
    printf("--- test_RidgeSolverConstructor")

   mtx <- matrix(1:9,nrow=3)   
   rownames(mtx) <- c("gene1","gene2","gene3")   
   solver <- RidgeSolver(mtx,targetGene = "gene1",
                            candidateRegulators = c("gene2","gene3"))
   
   checkEquals(class(solver)[1], "RidgeSolver")   
   checkTrue(all(c("RidgeSolver", "Solver") %in% is(solver)))   
}

# test_RidgeSolverConstructor
#----------------------------------------------------------------------------------------------------
test_ampAD.mef2c.154tfs.278samples.ridge <- function()
{
   printf("--- test_ampAD.mef2c.154tfs.278samples.ridge")

   # Load matrix and transform via arcsinh
   load(system.file(package="TReNA", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
   mtx.asinh <- asinh(mtx.sub)
   target.gene <- "MEF2C"
   tfs <- setdiff(rownames(mtx.asinh), "MEF2C")
   #print(fivenum(mtx.asinh)  # [1] 0.000000 1.327453 3.208193 4.460219 7.628290)

   set.seed(88117)
   ridge.solver <- RidgeSolver(mtx.asinh, target.gene, tfs)
   tbl <- run(ridge.solver)
   
   # Check for empirical values
   checkTrue(min(tbl$beta) > -0.15)
   checkTrue(max(tbl$beta) < 0.15)
   checkTrue(c("FOXP1") %in% rownames(subset(tbl, abs(beta) > 0.08)))

} # test_ampAD.mef2c.154tfs.278samples.ridge
#----------------------------------------------------------------------------------------------------
test_alpha.ridge <- function()
{
   printf("--- test_alpha.ridge")

   # Load matrix and transform via arcsinh
   load(system.file(package="TReNA", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
   mtx.asinh <- asinh(mtx.sub)
   target.gene <- "MEF2C"
   tfs <- setdiff(rownames(mtx.asinh), "MEF2C")
   #print(fivenum(mtx.asinh)  # [1] 0.000000 1.327453 3.208193 4.460219 7.628290)
   
   # check for expected non-sensical values
   # I think this is now mostly unnecessary
   #checkTrue(min(tbl$beta) < -7)
   #checkTrue(max(tbl$beta) > 10)

   set.seed(901188)
   ridge.solver <- RidgeSolver(mtx.asinh, target.gene, tfs, alpha = 0.1)
   tbl <- run(ridge.solver)

   # Check for empirical values
   checkTrue(nrow(tbl) < 40)
   checkTrue(min(tbl$beta) > -0.06)
   checkTrue(max(tbl$beta) < 0.1)
   checkTrue(c("SATB2") %in% rownames(subset(tbl, abs(beta) > 0.08)))

} # test_alpha.ridge
#----------------------------------------------------------------------------------------------------
test_lambda.ridge <- function()
{
   printf("--- test_lambda.ridge")

   # Load matrix and transform via arcsinh
   load(system.file(package="TReNA", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
   mtx.asinh <- asinh(mtx.sub)
   target.gene <- "MEF2C"
   tfs <- setdiff(rownames(mtx.asinh), "MEF2C")
   #print(fivenum(mtx.asinh)  # [1] 0.000000 1.327453 3.208193 4.460219 7.628290)
   
   # check for expected non-sensical values
   # I think this is now mostly unnecessary
   #checkTrue(min(tbl$beta) < -7)
   #checkTrue(max(tbl$beta) > 10)

   ridge.solver <- RidgeSolver(mtx.asinh, target.gene, tfs, lambda = 0.1)
   tbl <- run(ridge.solver)

   # Check for empirical values
   checkTrue(min(tbl$beta) > -0.15)
   checkTrue(max(tbl$beta) < 0.15)
   checkTrue(c("SATB2") %in% rownames(subset(tbl, abs(beta) > 0.08)))

} # test_lambda.ridge
#----------------------------------------------------------------------------------------------------
test_keep.metrics.ridge <- function()
{
   printf("--- test_keep.metrics.ridge")

   # Load matrix and transform via arcsinh
   load(system.file(package="TReNA", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
   mtx.asinh <- asinh(mtx.sub)
   target.gene <- "MEF2C"
   tfs <- setdiff(rownames(mtx.asinh), "MEF2C")
   #print(fivenum(mtx.asinh)  # [1] 0.000000 1.327453 3.208193 4.460219 7.628290)
   
   # check for expected non-sensical values
   # I think this is now mostly unnecessary
   #checkTrue(min(tbl$beta) < -7)
   #checkTrue(max(tbl$beta) > 10)

   set.seed(344)
   ridge.solver <- RidgeSolver(mtx.asinh, target.gene, tfs, keep.metrics=TRUE)
   tbl <- run(ridge.solver)

   # Check for empirical values
   checkTrue(min(tbl$mtx.beta$beta) > -0.1)
   checkTrue(max(tbl$mtx.beta$beta) < 0.1)
   checkTrue(c("FOXP1") %in% rownames(subset(tbl$mtx.beta, abs(beta) > 0.08)))
   checkTrue(tbl$lambda < 0.9)
   checkTrue(tbl$r2 > 0.95)

} # test_keep.metrics.ridge
#----------------------------------------------------------------------------------------------------
if(!interactive()) runTests()
