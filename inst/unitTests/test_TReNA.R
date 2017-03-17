library(TReNA)
library(RUnit)
#----------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_emptyConstructor()
   test_LassoSolverConstructor()
   test_RidgeSolverConstructor()
   test_RandomForestSolverConstructor()
   test_BayesSpikeSolverConstructor()
   test_SqrtLassoSolverConstructor()
   test_LassoPVSolverConstructor()
   test_PearsonSolverConstructor()
   test_SpearmanSolverConstructor()
   test_EnsembleSolverConstructor()
   
   test_developAndFitDummyTestData()
   test_fitDummyData()   

   test_ampAD.mef2c.154tfs.278samples.lasso()
   test_ampAD.mef2c.154tfs.278samples.bayesSpike()
#   test_ampAD.mef2c.154tfs.278samples.bayesSpike.nonCodingGenes()
   test_ampAD.mef2c.154tfs.278samples.randomForest()
   test_ampAD.mef2c.154tfs.278samples.ridge()
#   test_ampAD.mef2c.154tfs.278samples.sqrtlasso()
   test_ampAD.mef2c.154tfs.278samples.lassopv()
   test_ampAD.mef2c.154tfs.278samples.pearson()
   test_ampAD.mef2c.154tfs.278samples.spearman()
   test_ampAD.mef2c.154tfs.278samples.ensemble()

   test_scalePredictorPenalties.lasso()
   test_eliminateSelfTFs()

   test_NullFilter()
   test_VarianceFilter()
   test_FootprintFilter()

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
test_LassoSolverConstructor <- function()
{
   printf("--- test_LassoSolverConstructor")
   solver <- LassoSolver()
   checkEquals(getSolverName(solver), "LassoSolver")
   checkTrue(all(c("LassoSolver", "Solver") %in% is(solver)))

} # test_LassoSolverConstructor
#----------------------------------------------------------------------------------------------------
test_RandomForestSolverConstructor <- function()
{
   printf("--- test_RandomForestSolverConstructor")
   solver <- RandomForestSolver()
   checkEquals(getSolverName(solver), "RandomForestSolver")
   checkTrue(all(c("RandomForestSolver", "Solver") %in% is(solver)))

} # test_RandomForestSolverConstructor
#----------------------------------------------------------------------------------------------------
test_BayesSpikeSolverConstructor <- function()
{
   printf("--- test_BayesSpikeSolverConstructor")
   solver <- BayesSpikeSolver()
   checkEquals(getSolverName(solver), "BayesSpikeSolver")
   checkTrue(all(c("BayesSpikeSolver", "Solver") %in% is(solver)))

} # test_BayesSpikeSolverConstructor
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
test_ampAD.mef2c.154tfs.278samples.randomForest <- function()
{
   printf("--- test_ampAD.mef2c.154tfs.278samples.randomForest")

   load(system.file(package="TReNA", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
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
test_ampAD.mef2c.154tfs.278samples.bayesSpike.nonCodingGenes <- function()
{
    printf("--- test_ampAD.mef2c.154tfs.278samples.bayesSpike.nonCodingGenes")
    
    print(load(system.file(package="TReNA", "extdata/mtx.AD.noncodingNearPiez02.RData")))    
    target.genes <- genes.noncoding.near.piez02.active    
    mtx <- log2(mtx.nonCoding + 0.0001)    
    tfs <- setdiff(rownames(mtx), target.genes)    
    trena <- TReNA(mtx.assay=mtx, solver="bayesSpike", quiet=FALSE)
    
    findings <- list()
    
    for(target.gene in target.genes){        
        tbl <- solve(trena, target.gene, tfs)        
        tbl <- subset(tbl, pval < 0.01)       
        findings[[target.gene]] <- tbl
    }
    
    tbl.trimmed <- subset(tbl, abs(beta) > 0.1 & pval < 0.01)
    betas <- tbl.trimmed$beta
    big.abs.betas <- betas[abs(betas) > 1]
    checkTrue(length(big.abs.betas) > 20)
    checkTrue(nrow(tbl) > 10)
    checkTrue(cor(tbl.trimmed$beta, tbl.trimmed$gene.cor) < 0.2)
    
    mtx.tmp <- mtx.sub - min(mtx.sub) + 0.001
    mtx.log2 <- log2(mtx.tmp)
    fivenum(mtx.log2)  # [1] -9.9657843  0.8107618  3.6262014  5.4345771 10.0052973
    
    trena <- TReNA(mtx.assay=mtx.log2, solver="bayesSpike", quiet=FALSE)
    tfs <- setdiff(rownames(mtx.log2), "MEF2C")
    tbl2 <- solve(trena, target.gene, tfs)
    tbl2.trimmed <- subset(tbl2, abs(beta) > 0.1 & pval < 0.01)
    betas2 <- tbl2.trimmed$beta
    big.abs.betas2 <- betas2[abs(betas2) > 1]
    checkEquals(length(big.abs.betas2), 0)
    checkTrue(cor(tbl2.trimmed$beta, tbl2.trimmed$gene.cor) > 0.6)

} # test_ampAD.mef2c.154tfs.278samples.bayesSpike.nonCodingGenes
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
test_NullFilter <- function()
{
    printf("--- test_NullFilter")

    # Check that dummy data returns all gene names
    x <- test_developAndFitDummyTestData(quiet=TRUE)
    null.filter <- NullFilter(mtx.assay = x$assay)
    checkEquals(getCandidates(null.filter),rownames(x$assay))

} # test_NullFilter
#----------------------------------------------------------------------------------------------------
test_VarianceFilter <- function()
{
    printf("--- test_VarianceFilter")

    # Create dummy data and filter it based on variance
    x <- test_developAndFitDummyTestData(quiet=TRUE)
    var.filter <- VarianceFilter(mtx.assay = x$assay)
    checkTrue(length(getCandidates(var.filter, x$correlated.target)) > 0)

} # test_VarianceFilter
#----------------------------------------------------------------------------------------------------
test_SqrtLassoSolverConstructor <- function()
{
    printf("--- test_SqrtLassoSolverConstructor")

    # Construct the SqrtLassoSolver and check that it's correct
    solver <- SqrtLassoSolver()
    checkEquals(getSolverName(solver), "SqrtLassoSolver")
    checkTrue(all(c("SqrtLassoSolver", "Solver") %in% is(solver)))
}

# test_SqrtLassoSolverConstructor   
#----------------------------------------------------------------------------------------------------    
test_ampAD.mef2c.154tfs.278samples.sqrtlasso <- function()
{
   printf("--- test_ampAD.mef2c.154tfs.278samples.sqrtlasso")

   # Load matrix and transform via arcsinh
   load(system.file(package="TReNA", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
   target.gene <- "MEF2C"
   mtx.asinh <- asinh(mtx.sub)
   #print(fivenum(mtx.asinh)  # [1] 0.000000 1.327453 3.208193 4.460219 7.628290)
   
   trena <- TReNA(mtx.assay=mtx.asinh, solver="sqrtlasso", quiet=FALSE)
   tfs <- setdiff(rownames(mtx.asinh), "MEF2C")
   tbl <- solve(trena, target.gene, tfs)

   # Check for empirical values
   tbl <- tbl[order(abs(tbl$beta), decreasing=TRUE),, drop = FALSE]
   expected.genes <- sort(c("SATB2","STAT4","HLF","TSHZ3", "FOXP2"))
   actual.genes <- sort(rownames(subset(tbl, abs(beta) > 0.1)))
   printf("Top 5 genes: %s", paste(rownames(tbl[1:5,]),collapse=","))
   checkEquals(expected.genes,actual.genes)
   
} # test_ampAD.mef2c.154tfs.278samples.sqrtlasso
#----------------------------------------------------------------------------------------------------    
test_LassoPVSolverConstructor <- function()
{
    printf("--- test_LassoPvSolverConstructor")

    # Construct the SqrtLassoSolver and check that it's correct

    solver <- LassoPVSolver()
    checkEquals(getSolverName(solver), "LassoPVSolver")
    checkTrue(all(c("LassoPVSolver", "Solver") %in% is(solver)))
}

# test_LassoPVSolverConstructor
#----------------------------------------------------------------------------------------------------
test_ampAD.mef2c.154tfs.278samples.lassopv <- function()
{
   printf("--- test_ampAD.mef2c.154tfs.278samples.lassopv")

   # Load matrix and transform via arcsinh
   load(system.file(package="TReNA", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
   target.gene <- "MEF2C"
   mtx.asinh <- asinh(mtx.sub)
   #print(fivenum(mtx.asinh)  # [1] 0.000000 1.327453 3.208193 4.460219 7.628290)
   
   trena <- TReNA(mtx.assay=mtx.asinh, solver="lassopv", quiet=FALSE)
   tfs <- setdiff(rownames(mtx.asinh), "MEF2C")
   tbl <- solve(trena, target.gene, tfs)

   # Check for significant P-values; make sure they match the empirical value
   sig.genes <- tbl$p.values[tbl$p.values < 0.01]
   checkEquals(length(sig.genes),30)
   
} # test_ampAD.mef2c.154tfs.278samples.lassopv
#----------------------------------------------------------------------------------------------------
test_RidgeSolverConstructor <- function()
{
    printf("--- test_RidgeSolverConstructor")

    # Construct the RidgeSolver and check that it's correct

    solver <- RidgeSolver()
    checkEquals(getSolverName(solver), "RidgeSolver")
    checkTrue(all(c("RidgeSolver", "Solver") %in% is(solver)))
}

# test_RidgeSolverConstructor
#----------------------------------------------------------------------------------------------------
test_ampAD.mef2c.154tfs.278samples.ridge <- function()
{
   printf("--- test_ampAD.mef2c.154tfs.278samples.ridge")

   # Load matrix and transform via arcsinh
   load(system.file(package="TReNA", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
   target.gene <- "MEF2C"
   mtx.asinh <- asinh(mtx.sub)
   #print(fivenum(mtx.asinh)  # [1] 0.000000 1.327453 3.208193 4.460219 7.628290)
   
   trena <- TReNA(mtx.assay=mtx.asinh, solver="ridge", quiet=FALSE)
   tfs <- setdiff(rownames(mtx.asinh), "MEF2C")
   tbl <- solve(trena, target.gene, tfs)

   # Check for empirical values
   checkTrue(min(tbl$beta) > -0.15)
   checkTrue(max(tbl$beta) < 0.15)
   checkTrue(c("FOXP1") %in% rownames(subset(tbl, abs(beta) > 0.08)))

} # test_ampAD.mef2c.154tfs.278samples.ridge
#----------------------------------------------------------------------------------------------------
test_ampAD.mef2c.154tfs.278samples.pearson <- function()
{
   printf("--- test_ampAD.mef2c.154tfs.278samples.pearson")

   # Load matrix and transform via arcsinh
   load(system.file(package="TReNA", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
   target.gene <- "MEF2C"
   mtx.asinh <- asinh(mtx.sub)
   #print(fivenum(mtx.asinh)  # [1] 0.000000 1.327453 3.208193 4.460219 7.628290)
   
   trena <- TReNA(mtx.assay=mtx.asinh, solver="pearson", quiet=FALSE)
   tfs <- setdiff(rownames(mtx.asinh), "MEF2C")
   tbl <- solve(trena, target.gene, tfs)

   # Check for empirical values
   checkTrue(nrow(subset(tbl, abs(coefficient) > 0.8)) > 7)

} # test_ampAD.mef2c.154tfs.278samples.pearson
#----------------------------------------------------------------------------------------------------
test_ampAD.mef2c.154tfs.278samples.spearman <- function()
{
   printf("--- test_ampAD.mef2c.154tfs.278samples.spearman")

   # Load matrix and transform via arcsinh
   load(system.file(package="TReNA", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
   target.gene <- "MEF2C"
   mtx.asinh <- asinh(mtx.sub)
   #print(fivenum(mtx.asinh)  # [1] 0.000000 1.327453 3.208193 4.460219 7.628290)
   
   trena <- TReNA(mtx.assay=mtx.asinh, solver="spearman", quiet=FALSE)
   tfs <- setdiff(rownames(mtx.asinh), "MEF2C")
   tbl <- solve(trena, target.gene, tfs)

   # Check for empirical values
   checkTrue(nrow(subset(tbl, abs(coefficient) > 0.8)) > 7)

} # test_ampAD.mef2c.154tfs.278samples.spearman
#----------------------------------------------------------------------------------------------------
test_FootprintFilter <- function()
{
    printf("--- test_FootprintFilter")

    # Load ampAD data and filter it based on footprints
    load(system.file(package="TReNA", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
    footprint.filter <- FootprintFilter(mtx.assay = mtx.sub)

    db.address <- system.file(package="TReNA", "extdata")
    genome.db.uri <- paste("sqlite:/",db.address,"genome.sub.db", sep = "/")
    project.db.uri <- paste("sqlite:/",db.address,"project.sub.db", sep = "/")
    target.gene <- "MEF2C"
    tfs <- getCandidates(footprint.filter, target.gene,
                         genome.db.uri, project.db.uri,
                         size.upstream = 1000,
                         size.downstream = 1000)

    # Make sure it grabs the right number of genes
    checkEquals(length(tfs), 64)

} # test_FootprintFilter
#----------------------------------------------------------------------------------------------------
test_EnsembleSolverConstructor <- function()
{
    printf("--- test_EnsembleSolverConstructor")

    # Construct the EnsembleSolver and check that it's correct

    solver <- EnsembleSolver()
    checkEquals(getSolverName(solver), "EnsembleSolver")
    checkTrue(all(c("EnsembleSolver", "Solver") %in% is(solver)))
}

# test_EnsembleSolverConstructor
#----------------------------------------------------------------------------------------------------
test_ampAD.mef2c.154tfs.278samples.ensemble <- function()
{
   printf("--- test_ampAD.mef2c.154tfs.278samples.ensemble")

   set.seed(122113)
   # Load matrix and transform via arcsinh
   load(system.file(package="TReNA", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
   target.gene <- "MEF2C"
   mtx.asinh <- asinh(mtx.sub)
   #print(fivenum(mtx.asinh)  # [1] 0.000000 1.327453 3.208193 4.460219 7.628290)
   
   trena <- TReNA(mtx.assay=mtx.asinh, solver="ensemble", quiet=FALSE)
   tfs <- setdiff(rownames(mtx.asinh), "MEF2C")
   tbl <- solve(trena, target.gene, tfs)

   # Check for empirical values
   checkTrue(min(tbl$extr) > 1.2)
   checkTrue(max(tbl$extr) < 5.0)
   checkTrue(c("HLF") %in% tbl$gene)

} # test_ampAD.mef2c.154tfs.278samples.ensemble
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
