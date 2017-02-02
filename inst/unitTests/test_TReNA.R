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
   test_fitDREAM5_yeast.lasso_weighted.tfs()
   test_fitDREAM5_yeast.randomForest()
   test_fitDREAM5_yeast.bayesSpike()

   test_trainAndPredict_DREAM5_yeast.lasso()
   test_scalePredictorPenalties.lasso()
   test_eliminateSelfTFs()

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
   weights <- rep(1, nrow(features))
   fit = glmnet(features, target, weights=weights, lambda=lambda.min)
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
   mtx <- x$assay

     # log transform the data

   mtx <- mtx - min(mtx) + 0.001
   mtx <- log2(mtx)
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
   checkTrue(all(c(tf1, tf2) %in% rownames(tbl.betas)))
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

   tbl.betas <- as.data.frame(solve(trena, target.gene, tfs))

   #tbl.betas <- as.data.frame(result)
     # 1st row of tbl.betas is the intercept, give it an NA correlation
   #tbl.betas$cor <- c(NA, unlist(lapply(rownames(tbl.betas)[-1], function(tf) subset(tbl.gold.met2, TF==tf)$cor)))
   #colnames (tbl.betas) <- c("beta", "cor")

     # is there some rough correlation between the calculated betas and the
     # measured correlation?
   checkTrue(cor(tbl.betas$beta, tbl.betas$gene.cor) > 0.7)

} # test_fitDREAM5_yeast.lasso
#----------------------------------------------------------------------------------------------------
test_fitDREAM5_yeast.lasso_weighted.tfs <- function()
{
   printf("--- test_fitDREAM5_yeast.lasso_weighted.tfs")

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

   tf.weights <- c(100000, 1.0, 100000, 1.0)
   mtx.betas <- solve(trena, target.gene, tfs, tf.weights)

   tbl.betas <- as.data.frame(mtx.betas)
   significant.tfs <- rownames(tbl.betas)
   tbl.gold.met2.trimmed <- subset(tbl.gold.met2, TF %in% significant.tfs)

     # is there some rough correlation between the calculated betas and the
     # measured correlation?
   checkTrue(any(tbl.betas$gene.cor > 0.8))

} # test_fitDREAM5_yeast.lasso_weighted.tfs
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
      # RandomForest returns its own structured data object
      # we respect that here rather than squeeze it into a lasso-like table of beta coefficients
   rf.result <- solve(trena, target.gene, tfs)
   tbl.importance  <- rf.result$edges

     # is there some rough correlation between the importance
     # values returned by randomforest, and the directly
     # measured corrleation of the tfs to the target?

      # random forest results matrix is sorted by IncNodePurity
      # extract those values in the same order as they appear in tbl.gold
   rf.score <- tbl.importance[tbl.gold.met2$TF, "IncNodePurity"]
   checkTrue(cor(rf.score, tbl.gold.met2$cor) > 0.7)

} # test_fitDREAM5_yeast.randomForest
#----------------------------------------------------------------------------------------------------
test_trainAndPredict_DREAM5_yeast.lasso <- function()
{
   printf("--- test_trainAndPredict_DREAM5_yeast.lasso")
   load(system.file(package="TReNA", "extdata", "dream5.net4.yeast.RData"))
   checkTrue(exists("mtx"))
   checkTrue(exists("tbl.gold"))
   checkTrue(exists("tbl.ids"))

   checkEquals(dim(mtx), c(5950, 536))

   trena <- TReNA(mtx.assay=mtx, solver="lasso", quiet=FALSE)
   target.gene <- "MET2"
   tbl.gold.met2 <- subset(tbl.gold, target=="MET2")
   tfs <- tbl.gold.met2$TF

   count <- as.integer(0.80 * ncol(mtx))
   set.seed(31)
   training.samples <- sample(colnames(mtx), count)
   test.samples <- setdiff(colnames(mtx), training.samples)

   weights <- rep(1, length(tfs))
   model <- trainModel(trena, target.gene, tfs, tf.weightstraining.samples)
   prediction <- predictFromModel(trena, model, tfs, test.samples)

   agreement <- cor(mtx["MET2", test.samples], as.numeric(prediction))
   checkTrue(agreement > 0.8)

} # test_trainAndPredict_DREAM5_yeast.lasso
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
# locate all tfs mapped within 3k bases of MEF2C's transcription start site
# get their expression data
prepare_predict.mef2c.regulators <- function()
{
   print("--- prepare_predict.mef2c.regulators")
   tbl.candidates <- getFootprints("MEF2C", 10000, 10000)[, c("chr", "mfpStart", "mfpEnd", "motif", "tfs")] # 59 x 5
   table(tbl.candidates$motif)
       #  MA0103.2  MA0715.1 ZBTB16.p2
       #        14        43         2
   candidates <- sort(c(unique(tbl.candidates$tfs), "MEF2C"))

       # need ENSG ids for these gene symbols, in order to extract a small
       # expression matrix for testing

   if(!exists("ensembl.hg38"))
      ensembl.hg38 <<- useMart("ensembl", dataset="hsapiens_gene_ensembl")

    tbl.ids <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"),
                   filters="hgnc_symbol", values=candidates, mart=ensembl.hg38)
    deleters <- which(duplicated(tbl.ids$hgnc_symbol))
    if(length(deleters) > 0)
       tbl.ids <- tbl.ids[-deleters,]

    print(load("~/s/data/priceLab/AD/ampADMayo.64253genes.278samples.RData"))  # "mtx"
    gene.medians <- apply(mtx, 1, median)
    gene.sds     <- apply(mtx, 1, sd)
    median.keepers <- which(gene.medians > 0.4)
    sd.keepers     <- which(gene.sds > 0.1)
    keepers <- sort(c(median.keepers, sd.keepers))
    mtx.keep <- mtx[keepers,]
    ens.goi <- intersect(tbl.ids$ensembl_gene_id, rownames(mtx.keep)) # 154
    mtx.sub <- mtx.keep[ens.goi,]   # 154 genes, 278 samples
    tbl.ids <- subset(tbl.ids, ensembl_gene_id %in% ens.goi)
    rownames(mtx.sub) <- tbl.ids$hgnc_symbol
    filename <- sprintf("../extdata/ampAD.%dgenes.mef2cTFs.%dsamples.RData", nrow(mtx.sub),
                        ncol(mtx.sub))
    printf("saving mtx.sub to %s", filename)
    save(mtx.sub, file=filename)

} # prepare_predict.mef2c.regulators
#----------------------------------------------------------------------------------------------------
test_predict.mef2c.regulators <- function()
{
   print("--- test_predict.mef2c.regulators")
   if(!exists("mtx.sub"))
       load(system.file(package="TReNA", "extdata/ampAD.58genes.mef2cTFs.278samples.RData"))
   if(!exists("tbl.mef2c.candidates"))
      tbl.mef2c.candidates <<- getFootprints("MEF2C", 3000, 0)[, c("chr", "mfpStart", "mfpEnd", "motif", "tfs")] # 59 x 5
   #candidates.sub <- subset(tbl.mef2c.candidates, motif != "MA0715.1")$tfs
   #mtx.sub <- mtx.sub[c(candidates.sub, "MEF2C"),]

   # mtx.sub <- log10(mtx.sub + 0.001)

   target.gene <- "MEF2C"

   trena <- TReNA(mtx.assay=mtx.sub, solver="lasso", quiet=FALSE)
   tfs <- setdiff(rownames(mtx.sub), "MEF2C")
   result <- solve(trena, target.gene, tfs)

} # test_predict.mef2c.regulators
#----------------------------------------------------------------------------------------------------
test_fitDREAM5_yeast.bayesSpike <- function()
{
   printf("--- test_fitDREAM5_yeast.bayesSpike")
   load(system.file(package="TReNA", "extdata", "dream5.net4.yeast.RData"))
   checkTrue(exists("mtx"))
   checkTrue(exists("tbl.gold"))
   checkTrue(exists("tbl.ids"))

   checkEquals(dim(mtx), c(5950, 536))

   trena <- TReNA(mtx.assay=mtx, solver="bayesSpike", quiet=FALSE)

     # subset(tbl.gold, target=="MET2")
     #     TF target score        cor
     #   CBF1   MET2     1 -0.4746397
     #  MET32   MET2     1  0.8902950
     #  MET31   MET2     1  0.1245628
     #   MET4   MET2     1  0.5301484

   target.gene <- "MET2"
   tbl.gold.met2 <- subset(tbl.gold, target=="MET2")
   tfs <- tbl.gold.met2$TF

   tbl <- solve(trena, target.gene, tfs)
   #tbl.betas <- data.frame(beta=result$beta, pval=result$pval, z=result$z, post=result$post)
   #rownames(tbl.betas) <- tfs
   #tbl.betas$score <- -log10(tbl.betas$pval)
   #tbl.betas$cor <- tbl.gold.met2$cor

     # is there some rough correlation between the calculated betas and the
     # measured correlation?
   checkTrue(with(tbl, cor(beta,gene.cor)) > 0.9)

} # test_fitDREAM5_yeast.bayesSpike
#----------------------------------------------------------------------------------------------------

test_eliminateSelfTFs <- function()
{
   printf("--- test_eliminateSelfTFs")

   load(system.file(package="TReNA", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
   target.gene <- "MEF2C"

   mtx.tmp <- mtx.sub - min(mtx.sub) + 0.001
   mtx.log2 <- log2(mtx.tmp)

   trena <- TReNA(mtx.assay=mtx.log2, solver="lasso", quiet=FALSE)
   tfs <- rownames(mtx.log2)
   checkTrue(target.gene %in% tfs)         # our test case
   tbl.betas <- solve(trena, target.gene, tfs)
   checkTrue(!target.gene %in% rownames(tbl.betas))
   checkTrue(cor(tbl.betas$beta[1:10], tbl.betas$gene.cor[1:10]) > 0.6)

   trena2 <- TReNA(mtx.assay=mtx.log2, solver="bayesSpike", quiet=FALSE)
   tbl.betas2 <- solve(trena2, target.gene, tfs)
   checkTrue(!target.gene %in% rownames(tbl.betas2))
   checkTrue(cor(tbl.betas2$beta[1:10], tbl.betas2$gene.cor[1:10]) > 0.7)

} # test_eliminateSelfTFs
#----------------------------------------------------------------------------------------------------

test_PearsonSolverConstructor <- function()
{
    printf("--- test_PearsonSolverConstructor")
    solver <- PearsonSolver()
    checkEquals(getSolverName(solver), "PearsonSolver")
    checkTrue(all(c("PearsonSolver", "Solver") %in% is(solver)))

} # test_PearsonSolverConstructor
#----------------------------------------------------------------------------------------------------
test_SpearmanSolverConstructor <- function()
{
    printf("--- test_SpearmanSolverConstructor")
    solver <- SpearmanSolver()
    checkEquals(getSolverName(solver), "SpearmanSolver")
    checkTrue(all(c("SpearmanSolver", "Solver") %in% is(solver)))

} # test_SpearmanSolverConstructor
#----------------------------------------------------------------------------------------------------
if(!interactive()) runTests()
