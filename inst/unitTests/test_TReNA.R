library(TReNA)
library(RUnit)
#----------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_emptyConstructor()
   test_developAndFitDummyTestData()
   test_LassoSolverConstructor()
   test_fitDummyData()
   test_fitDREAM5_yeast.lasso()
   test_fitDREAM5_yeast.randomForest()

} # runTests
#----------------------------------------------------------------------------------------------------
test_emptyConstructor <- function()
{
   printf("--- test_emptyConstructor")
   trena <- TReNA()

   mtx <- getAssayData(trena)
   checkTrue("matrix" %in% is(mtx))
   checkEquals(dim(mtx), c(1,1))   # the default (and implicilty empty) matrix
   checkTrue(is.na(mtx[1,1]))

   #mtx.priors <- getPriors(trena)
   #checkTrue("matrix" %in% is(mtx.priors))
   #checkEquals(dim(mtx.priors), c(1,1))

} # test_emptyConstructor
#----------------------------------------------------------------------------------------------------
test_LassoSolverConstructor <- function()
{
   printf("--- test_LassoSolverConstructor")
   solver <- LassoSolver()
   checkEquals(getSolverName(solver), "LassoSolver")
   checkTrue(all(c("LassoSolver", "Solver") %in% is(solver)))

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

   mtx[, target.gene] <- jitter(mtx[, TF.1], amount=10)
   mtx[, TF.2] <- mtx[, TF.1] - mtx[, target.gene]

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
   fit = glmnet(features, target, lambda=lambda.min)
       # extract the exponents of the fit
   betas <- as.matrix(t(coef(cv.out, s="lambda.min")))

       # only TF.1 should contribute to a model of the target gene
   checkTrue(betas[1, "(Intercept)"] > 1)
   checkTrue(betas[1, TF.1] > 0.9)
   checkTrue(betas[1, TF.2] < -0.5)

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
   assay <- x$assay
   target.gene <- x$correlated.target
   tfs <- x$tf.genes

   trena <- TReNA(mtx.assay=assay, solver="lasso", quiet=FALSE)

   tf1 <- x$correlated.tfs[1]
   tf2 <- x$correlated.tfs[2]
   target.gene <- x$correlated.target

   target.values <- as.numeric(x$assay[target.gene,])
   tf1.values    <- as.numeric(x$assay[tf1,])
   tf2.values    <- as.numeric(x$assay[tf2,])

     # we expect an intercept and a coef for tfs gene.02 and gene.03
     # which predict the value of the target.gene

   result <- solve(trena, target.gene, tfs)
   checkEquals(sort(names(result)), c("rSquared", "scores"))
   betas <- result$scores
   rSquared <- result$rSquared
   checkTrue(all(c("(Intercept)", tf1, tf2) %in% rownames(betas)))
   intercept <- betas["(Intercept)", 1]
   coef.tf1  <- betas[tf1, 1]
   coef.tf2  <- betas[tf2, 1]
   predicted <- intercept + (coef.tf1 * assay[tf1,]) + (coef.tf2 * assay[tf2,])
   actual    <- assay[target.gene, ]

      # on average, the prediction should be reasonable
   checkEqualsNumeric(sum(actual - predicted), 0, tol=1e-8)

      # but there is lots of slop in the prediction: no overfitting here!
   checkTrue(sum(abs(actual-predicted)) > 30)

   checkTrue(mean(rSquared[, "s0"]) > 0.98)

} # test_fitDummyData
#----------------------------------------------------------------------------------------------------
test_fitDREAM5_yeast.lasso <- function()
{
   printf("--- test_fitDREAM5_yeast.lasso")
   load(system.file(package="TReNA", "extdata", "dream5.net4.yeast.RData"))
   checkTrue(exists("mtx"))
   checkTrue(exists("tbl.gold"))
   checkTrue(exists("tbl.ids"))

   checkEquals(dim(mtx), c(5950, 536))

   trena <- TReNA(mtx.assay=mtx, solver="lasso", quiet=FALSE)
     # subset(tbl.gold, target=="MET2")
     #     TF target score        cor
     #   CBF1   MET2     1 -0.4746397
     #  MET32   MET2     1  0.8902950
     #  MET31   MET2     1  0.1245628
     #   MET4   MET2     1  0.5301484

   target.gene <- "MET2"
   tbl.gold.met2 <- subset(tbl.gold, target=="MET2")
   tfs <- tbl.gold.met2$TF
   result <- solve(trena, target.gene, tfs)
   tbl.betas <- as.data.frame(result$scores)
   tbl.betas$cor <- c(NA, tbl.gold.met2$cor)
   colnames (tbl.betas) <- c("beta", "cor")

     # is there some rough correlation between the calculated betas and the
     # measured correlation?
   checkTrue(with(tbl.betas[2:5, ], cor(beta,cor)) > 0.9)

} # test_fitDREAM5_yeast.lasso
#----------------------------------------------------------------------------------------------------
test_fitDREAM5_yeast.randomForest <- function()
{
   printf("--- test_fitDREAM5_yeast.randomForest")
   load(system.file(package="TReNA", "extdata", "dream5.net4.yeast.RData"))
   checkTrue(exists("mtx"))
   checkTrue(exists("tbl.gold"))
   checkTrue(exists("tbl.ids"))

   checkEquals(dim(mtx), c(5950, 536))

   trena <- TReNA(mtx.assay=mtx, solver="randomForest", quiet=FALSE)
     # subset(tbl.gold, target=="MET2")
     #     TF target score        cor
     #   CBF1   MET2     1 -0.4746397
     #  MET32   MET2     1  0.8902950
     #  MET31   MET2     1  0.1245628
     #   MET4   MET2     1  0.5301484

   target.gene <- "MET2"
   tbl.gold.met2 <- subset(tbl.gold, target=="MET2")
   tfs <- tbl.gold.met2$TF
   result <- solve(trena, target.gene, tfs)


   tbl.purity <- as.data.frame(result$scores)

     # is there some rough correlation between the importance
     # values returned by randomforest, and the directly
     # measured corrleation of the tfs to the tarte

   checkTrue(cor(tbl.purity$IncNodePurity, tbl.gold.met2$cor) > 0.7)

} # test_fitDREAM5_yeast.randomForest
#----------------------------------------------------------------------------------------------------
