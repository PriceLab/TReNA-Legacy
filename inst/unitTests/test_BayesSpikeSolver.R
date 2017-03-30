library(TReNA)
library(RUnit)
#----------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_BayesSpikeSolverConstructor()
   test_ampAD.mef2c.154tfs.278samples.bayesSpike()
   
} # runTests
#----------------------------------------------------------------------------------------------------
test_BayesSpikeSolverConstructor <- function()
{
   printf("--- test_BayesSpikeSolverConstructor")
   solver <- BayesSpikeSolver()
   checkEquals(getSolverName(solver), "BayesSpikeSolver")
   checkTrue(all(c("BayesSpikeSolver", "Solver") %in% is(solver)))

} # test_BayesSpikeSolverConstructor
#----------------------------------------------------------------------------------------------------
test_ampAD.mef2c.154tfs.278samples.bayesSpike <- function()
{
   printf("--- test_ampAD.mef2c.154tfs.278samples.bayesSpike")

   load(system.file(package="TReNA", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
   target.gene <- "MEF2C"
   mtx.asinh <- asinh(mtx.sub)
   #print(fivenum(mtx.asinh)  # [1] 0.000000 1.327453 3.208193 4.460219 7.628290)

   trena <- TReNA(mtx.assay=mtx.asinh, solver="bayesSpike", quiet=FALSE)
   tfs <- setdiff(rownames(mtx.asinh), "MEF2C")
   tbl <- solve(trena, target.gene, tfs)
   tbl.trimmed <- subset(tbl, abs(beta) > 0.1 & pval < 0.01)

   # I don't think we need these anymore
   #betas <- tbl.trimmed$beta
   #big.abs.betas <- betas[abs(betas) > 1]
   #checkTrue(length(big.abs.betas) > 20)

   # Check number of results and correlation of results
   checkTrue(nrow(tbl.trimmed) > 6)
   checkTrue(cor(tbl.trimmed$beta, tbl.trimmed$gene.cor) > 0.6)

} # test_ampAD.mef2c.154tfs.278samples.bayesSpike
#----------------------------------------------------------------------------------------------------
if(!interactive()) runTests()
