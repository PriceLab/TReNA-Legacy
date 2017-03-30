library(TReNA)
library(RUnit)
#----------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_emptyConstructor()   
   test_developAndFitDummyTestData()
   test_fitDummyData()   
   test_eliminateSelfTFs()
   test_MatrixWarnings()
   
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
test_eliminateSelfTFs <- function()
{
   printf("--- test_eliminateSelfTFs")

   set.seed(10045)
   load(system.file(package="TReNA", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
   target.gene <- "MEF2C"

   mtx.asinh <- asinh(mtx.sub)

   trena <- TReNA(mtx.assay=mtx.asinh, solver="lasso", quiet=FALSE)
   tfs <- rownames(mtx.asinh)
   checkTrue(target.gene %in% tfs)         # our test case
   tbl.betas <- solve(trena, target.gene, tfs)
   checkTrue(!target.gene %in% rownames(tbl.betas))
   checkTrue(cor(tbl.betas$beta, tbl.betas$gene.cor) > 0.6)

   trena2 <- TReNA(mtx.assay=mtx.asinh, solver="bayesSpike", quiet=FALSE)
   tbl.betas2 <- solve(trena2, target.gene, tfs)
   checkTrue(!target.gene %in% rownames(tbl.betas2))
   checkTrue(cor(tbl.betas2$beta[1:10], tbl.betas2$gene.cor[1:10]) > 0.6)

} # test_eliminateSelfTFs
#----------------------------------------------------------------------------------------------------
test_MatrixWarnings <- function()
{
    printf("--- test_MatrixWarnings")

    # Change warnings to errors
    options(warn = 2)

    # Check that a skewed matrix produces an error
    test.mtx <- matrix(1:10000, nrow = 100)
    test.mtx[1,1] <- 1e7
    checkException(TReNA(test.mtx), silent = TRUE)

    # Check that a matrix with a row of 0's produces an error for most solvers
    test.mtx[1,] <- 0
    checkException(TReNA(test.mtx, solver = "bayesSpike"), silent = TRUE)
    checkException(TReNA(test.mtx, solver = "lassopv"), silent = TRUE)
    checkException(TReNA(test.mtx, solver = "sqrtlasso"), silent = TRUE)
    checkException(TReNA(test.mtx, solver = "randomForest"), silent = TRUE)
    checkException(TReNA(test.mtx, solver = "pearson"), silent = TRUE)
    checkException(TReNA(test.mtx, solver = "spearman"), silent = TRUE)

    # Check that a target gene with low expression causes a warning for a solver
    test.mtx[1,] <- 0.1
    rownames(test.mtx) <- 1:100
    target.gene <- 1
    tfs <- 2:100
    trena <- TReNA(test.mtx, solver = "ensemble")
    checkException(solve(trena, target.gene, tfs), silent = TRUE)    
    
    # Change warnings back to warnings
    options(warn = 1)

} #test_MatrixWarnings
#----------------------------------------------------------------------------------------------------
if(!interactive()) runTests()
