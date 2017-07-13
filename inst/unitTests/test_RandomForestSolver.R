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
   test_RandomForestSolverFewCandidates()
   test_ampAD.mef2c.154tfs.278samples.randomForest()

} # runTests
#----------------------------------------------------------------------------------------------------
test_RandomForestSolverConstructor <- function()
{
    printf("--- test_RandomForestSolverConstructor")

    mtx <- matrix(1:9,nrow=3)   
    rownames(mtx) <- c("gene1","gene2","gene3")    
    solver <- RandomForestSolver(mtx,targetGene = "gene1",                          
                          candidateRegulators = c("gene2","gene3"))    
   
    checkEquals(class(solver)[1], "RandomForestSolver")    
    checkTrue(all(c("RandomForestSolver", "Solver") %in% is(solver)))
    
} # test_RandomForestSolverConstructor
#----------------------------------------------------------------------------------------------------
test_RandomForestSolverFewCandidates <- function()
{
    printf("--- test_RandomForestSolverFewCandidates")
    set.seed(17)
   solver <- RandomForestSolver(mtx, targetGene="MEF2C", candidateRegulators=candidate.tfs)
   x <- run(solver)
   checkTrue(all(c("edges", "r2") %in% names(x)))
   tbl <- x$edges
   checkEquals(dim(tbl), c(3, 2))
   checkEquals(colnames(tbl), c("IncNodePurity", "gene.cor"))
   checkEquals(rownames(tbl), c("ATF7", "NR3C2", "PRRX1"))

} # test_RandomForestSolverFewCandidates
#----------------------------------------------------------------------------------------------------
test_ampAD.mef2c.154tfs.278samples.randomForest <- function()
{
   printf("--- test_ampAD.mef2c.154tfs.278samples.randomForest")

   set.seed(334)
   targetGene <- "MEF2C"
   candidate.tfs <- setdiff(rownames(mtx), targetGene)
   solver <- RandomForestSolver(mtx, targetGene=targetGene, candidateRegulators=candidate.tfs)
   x <- run(solver)
   tbl <- x$edges
      # check just the highest scores
   tbl.10 <- subset(tbl, IncNodePurity > 10)
   checkEquals(rownames(tbl.10), c("HLF", "STAT4", "SATB2", "SATB1"))

} # test_ampAD.mef2c.154tfs.278samples.randomForest
#----------------------------------------------------------------------------------------------------
if(!interactive()) runTests()
