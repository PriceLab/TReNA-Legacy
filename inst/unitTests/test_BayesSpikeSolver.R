library(TReNA)
library(RUnit)
#----------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_BayesSpikeSolverConstructor()
   test_ampAD.mef2c.154tfs.278samples.bayesSpike()
   test_nOrderings()
   
} # runTests
#----------------------------------------------------------------------------------------------------
test_BayesSpikeSolverConstructor <- function()
{
   printf("--- test_BayesSpikeSolverConstructor")

   mtx <- matrix(1:9,nrow=3)   
   rownames(mtx) <- c("gene1","gene2","gene3")   
   solver <- BayesSpikeSolver(mtx,targetGene = "gene1",
                            candidateRegulators = c("gene2","gene3"))
   
   checkEquals(class(solver)[1], "BayesSpikeSolver")   
   checkTrue(all(c("BayesSpikeSolver", "Solver") %in% is(solver)))   

} # test_BayesSpikeSolverConstructor
#----------------------------------------------------------------------------------------------------
test_ampAD.mef2c.154tfs.278samples.bayesSpike <- function()
{
   printf("--- test_ampAD.mef2c.154tfs.278samples.bayesSpike")

   set.seed(12415)
   load(system.file(package="TReNA", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
   
   mtx.asinh <- asinh(mtx.sub)
   target.gene <- "MEF2C"
   tfs <- setdiff(rownames(mtx.asinh), "MEF2C")
   #print(fivenum(mtx.asinh)  # [1] 0.000000 1.327453 3.208193 4.460219 7.628290)

   bayes.solver <- BayesSpikeSolver(mtx.asinh,target.gene,tfs)
   tbl <- run(bayes.solver)
   tbl.trimmed <- subset(tbl, abs(beta) > 0.1 & pval < 0.01)

   #betas <- tbl.trimmed$beta
   #big.abs.betas <- betas[abs(betas) > 1]
   #checkTrue(length(big.abs.betas) > 20)

   # Check number of results and correlation of results
   checkTrue(nrow(tbl.trimmed) == 12)
   checkTrue(cor(tbl.trimmed$beta, tbl.trimmed$gene.cor) > 0.8)

} # test_ampAD.mef2c.154tfs.278samples.bayesSpike
#----------------------------------------------------------------------------------------------------
test_nOrderings <- function()
{
   printf("--- test_nOrderings")

   set.seed(12415)
   load(system.file(package="TReNA", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
   
   mtx.asinh <- asinh(mtx.sub)
   target.gene <- "MEF2C"
   tfs <- setdiff(rownames(mtx.asinh), "MEF2C")
   #print(fivenum(mtx.asinh)  # [1] 0.000000 1.327453 3.208193 4.460219 7.628290)

   # Use 100 orderings instead of 10
   bayes.solver <- BayesSpikeSolver(mtx.asinh,target.gene,tfs,nOrderings = 100)
   tbl <- run(bayes.solver)
   tbl.trimmed <- subset(tbl, abs(beta) > 0.1 & pval < 0.01)

   #betas <- tbl.trimmed$beta
   #big.abs.betas <- betas[abs(betas) > 1]
   #checkTrue(length(big.abs.betas) > 20)

   # Check number of results and correlation of results
   checkTrue(nrow(tbl.trimmed) == 10)
   checkTrue(cor(tbl.trimmed$beta, tbl.trimmed$gene.cor) > 0.75)

} # test_nOrderings
#----------------------------------------------------------------------------------------------------
if(!interactive()) runTests()
