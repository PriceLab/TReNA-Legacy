library(TReNA)
library(RUnit)
#----------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_LassoSolverConstructor()   
   test_developAndFitDummyTestData()
   test_fitDummyData()   
   test_ampAD.mef2c.154tfs.278samples.lasso()
   test_scalePredictorPenalties.lasso()
   
} # runTests
#----------------------------------------------------------------------------------------------------
test_LassoSolverConstructor <- function()
{
   printf("--- test_LassoSolverConstructor")
   #solver <- LassoSolver()
   solver <- TReNA(matrix(), solver = "lasso")
   checkEquals(getSolverName(solver), "LassoSolver")
   #checkTrue(all(c("LassoSolver", "TReNA") %in% is(solver)))

} # test_LassoSolverConstructor
#----------------------------------------------------------------------------------------------------
test_developAndFitDummyTestData <- function(quiet=FALSE)
{
   if(!quiet)
      printf("--- test_developAndFitDummyTestData")

   set.seed(37)

   gene.count <- 50
   sample.count <- 20

   mtx <- matrix(100 * abs(rnorm(sample.count * gene.count)), sample.count, gene.count)
   colnames(mtx) <- sprintf("gene.%02d", 1:ncol(mtx))
   rownames(mtx) <- sprintf("samp.%02d", 1:nrow(mtx))

      # arbitrarily designate 10 of the genes as transcription factors
   TF.genes <- sort(sprintf("gene.%02d", sample(1:gene.count, 10)))
   target.genes <- setdiff(colnames(mtx), TF.genes)
   TF.1 <- TF.genes[1]
   TF.2 <- TF.genes[2]
   target.gene <- target.genes[1]

   mtx[, target.gene] <- jitter((mtx[, TF.1]+mtx[, TF.2]), amount=10)
   mtx[, TF.2] <- (mtx[, TF.1] - mtx[, target.gene])

   # make sure that the target is the sum of the two TFs
   checkTrue(all( mtx[, target.gene] == mtx[, TF.1] - mtx[, TF.2]))

   # make sure other correlations are low
   exclude.these.columns <- unlist(lapply(c(TF.1, TF.2, target.gene), function(g) grep(g, colnames(mtx))))
   mtx.sub <- mtx[, -exclude.these.columns]
   other.correlations <- apply(mtx.sub, 2, function(col) cor(col, mtx[, TF.1]))
      # random chance could produce another gene well-correlated to TF.1, but with our seed has not
   checkTrue(max(other.correlations) < 0.5)

   target.col <- grep(target.gene, colnames(mtx))
   target   <- mtx[,  target.col]
   features <- mtx[ , -target.col]

       # learn lambda.min
   cv.out <- cv.glmnet(features, target, grouped=FALSE)
   #suppressWarnings(cv.out <- cv.glmnet(features, target, grouped=FALSE))
   lambda.min <- cv.out$lambda.min
   weights <- rep(1, nrow(features))
   fit = glmnet(features, target, weights=weights, lambda=lambda.min)
       # extract the exponents of the fit
   betas <- as.matrix(t(coef(cv.out, s="lambda.min")))

       # only TF.1 should contribute to a model of the target gene
   checkTrue(betas[1, "(Intercept)"] > 1)
   checkTrue(betas[1, TF.1] > 0.9)
   checkTrue(betas[1, TF.2] < -0.9)

      # return this for other tests to use.
      # learned belatedly:  genes as rownames, samples as colnames is the standard
      # so transpose this matrix before returning it
   invisible(list(assay=t(mtx), tf.genes=TF.genes, target.genes=target.genes,
                  correlated.tfs=c(TF.1, TF.2), correlated.target=target.gene))

} # test_developAndFitDummyTestData
#----------------------------------------------------------------------------------------------------
test_fitDummyData <- function()
{
   printf("--- test_fitDummyData")

   x <- test_developAndFitDummyTestData(quiet=TRUE)
   mtx <- x$assay

   #asinh-transform the data   
   mtx <- asinh(mtx)  
   target.gene <- x$correlated.target
   tfs <- x$tf.genes

   trena <- TReNA(mtx.assay=mtx, solver="lasso", quiet=FALSE)

   tf1 <- x$correlated.tfs[1]
   tf2 <- x$correlated.tfs[2]
   target.gene <- x$correlated.target

   target.values <- as.numeric(x$assay[target.gene,])
   tf1.values    <- as.numeric(x$assay[tf1,])
   tf2.values    <- as.numeric(x$assay[tf2,])

     # we expect an intercept and a coef for tfs gene.02 and gene.03
     # which predict the value of the target.gene

   tbl.betas <- solve(trena, target.gene, tfs, extraArgs =list(alpha=1.0, lambda=NULL))
   checkTrue(all(c(tf1,tf2) %in% rownames(tbl.betas)))
   checkEquals(colnames(tbl.betas), c("beta", "intercept", "gene.cor"))
   intercept <- tbl.betas[1, "intercept"]
   coef.tf1  <- tbl.betas[tf1, "beta"]
   coef.tf2  <- tbl.betas[tf2, "beta"]
   predicted <- intercept + (coef.tf1 * mtx[tf1,]) + (coef.tf2 * mtx[tf2,])
   actual    <- mtx[target.gene, ]

      # on average, the prediction should be reasonable
   checkEqualsNumeric(sum(actual - predicted), 0, tol=1e-8)

} # test_fitDummyData
#----------------------------------------------------------------------------------------------------
test_ampAD.mef2c.154tfs.278samples.lasso <- function()
{
   printf("--- test_ampAD.mef2c.154tfs.278samples.lasso")

   # Load matrix and transform via arcsinh
   load(system.file(package="TReNA", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
   target.gene <- "MEF2C"
   mtx.asinh <- asinh(mtx.sub)
   #print(fivenum(mtx.asinh)  # [1] 0.000000 1.327453 3.208193 4.460219 7.628290)
   
   # check for expected non-sensical values
   # I think this is now mostly unnecessary
   #checkTrue(min(tbl$beta) < -7)
   #checkTrue(max(tbl$beta) > 10)

   trena <- TReNA(mtx.assay=mtx.asinh, solver="lasso", quiet=FALSE)
   tfs <- setdiff(rownames(mtx.asinh), "MEF2C")
   tbl <- solve(trena, target.gene, tfs)

   # Check for empirical values
   checkTrue(min(tbl$beta) > -0.1)
   checkTrue(max(tbl$beta) < 0.3)
   checkTrue(c("SATB2") %in% rownames(subset(tbl, abs(beta) > 0.15)))

} # test_ampAD.mef2c.154tfs.278samples.lasso
#----------------------------------------------------------------------------------------------------
# one possible source of down-weighting data from TFs is the frequency of their putative
# binding sites across the genome.  the SP1-n family has a motif-in-footprint about every
# 5k, for example.
# db <- dbConnect(dbDriver("SQLite"), "~/github/snpFoot/inst/misc/sqlite.explorations/fpTf.sqlite")
# tbl.fpTfFreqs <- dbGetQuery(db, "select * from fpTfFreqs")
# as.integer(fivenum(values))  # [1]    241   4739   9854  22215 658334
test_scalePredictorPenalties.lasso <- function()
{
   printf("--- test_scalePredictorPenalties.lasso")
   raw.values <-  c(241, 4739, 9854, 22215, 658334)
   ls <- LassoSolver(matrix())
   min.observed <- 1        # just one footprint in the genome for some possible gene
   max.observed <- 658334   # max observed putative binding sites for SPx family of tfs

   cooked.values <- rescalePredictorWeights(ls, rawValue.min=1, rawValue.max=1000000, raw.values)
   checkEqualsNumeric(cooked.values[1], 0.99976, tol=1e-3)
   checkEqualsNumeric(cooked.values[2], 0.99526, tol=1e-3)
   checkEqualsNumeric(cooked.values[3], 0.99015, tol=1e-3)
   checkEqualsNumeric(cooked.values[4], 0.97778, tol=1e-3)
   checkEqualsNumeric(cooked.values[5], 0.34166, tol=1e-3)

} # test_scalePredictorPenalties.lasso
#----------------------------------------------------------------------------------------------------
if(!interactive()) runTests()
