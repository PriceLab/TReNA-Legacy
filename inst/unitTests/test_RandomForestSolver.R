library(TReNA)
library(RUnit)
#----------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
if(!exists("mtx")){
   load(system.file(package="TReNA", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
   mtx <- asinh(mtx.sub)
   }

candidate.tfs <- c("ATF7", "NR3C2", "MAFB", "PRRX1", "E2F8", "XBP1"); # from gtex skin data(!)
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_RandomForestSolverConstructor()
   test_ampAD.mef2c.154tfs.278samples.randomForest()

} # runTests
#----------------------------------------------------------------------------------------------------
test_RandomForestSolverConstructor <- function()
{
   printf("--- test_RandomForestSolverConstructor")
   solver <- RandomForestSolver(mtx, targetGene="MEF2C", candidateRegulators=candidate.tfs)

   solver <- TReNA(matrix(), solver = "randomForest")
   checkEquals(getSolverName(solver), "RandomForestSolver")
   #checkTrue(all(c("RandomForestSolver", "Solver") %in% is(solver)))

} # test_RandomForestSolverConstructor
#----------------------------------------------------------------------------------------------------
test_ampAD.mef2c.154tfs.278samples.randomForest <- function()
{
   printf("--- test_ampAD.mef2c.154tfs.278samples.randomForest")

   target.gene <- "MEF2C"
   mtx.asinh <- asinh(mtx.sub)
   #print(fivenum(mtx.asinh)  # [1] 0.000000 1.327453 3.208193 4.460219 7.628290)

   # Set the seed and solve
   set.seed(10101)
   trena <- TReNA(mtx.assay=mtx.asinh, solver="randomForest", quiet=FALSE)
   tfs <- setdiff(rownames(mtx.asinh), "MEF2C")
   rf.result <- solve(trena, target.gene, tfs)
   tbl.scores <- rf.result$edges

   tbl.scores <- tbl.scores[order(tbl.scores$IncNodePurity, decreasing=TRUE),, drop=FALSE]

   # a loose test, ignoring rank of these 7 genes for now
   actual.genes.reported <- sort(rownames(subset(tbl.scores, IncNodePurity > 2.8)))
   # Updated to reflect asinh results
   expected.genes <- sort(c("HLF", "STAT4", "SATB1", "SATB2", "FOXP2", "DRGX","TSHZ3"))
   printf("expected: %s", paste(expected.genes, collapse=","))
   printf("actual: %s", paste(actual.genes.reported, collapse=","))
   checkEquals(actual.genes.reported, expected.genes)


} # test_ampAD.mef2c.154tfs.278samples.randomForest
#----------------------------------------------------------------------------------------------------
if(!interactive()) runTests()
