library(TReNA)
library(RUnit)
#----------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_ampAD.mef2c.154tfs.278samples.lasso()
   test_ampAD.mef2c.154tfs.278samples.bayesSpike()
#   test_ampAD.mef2c.154tfs.278samples.bayesSpike.nonCodingGenes()
   test_ampAD.mef2c.154tfs.278samples.randomForest()
   test_LCLs.build_genomewide_model.lasso()
   test_LCLs.build_genomewide_model.randomForest()

} # runTests
#----------------------------------------------------------------------------------------------------
test_ampAD.mef2c.154tfs.278samples.lasso <- function()
{
   printf("--- test_ampAD.mef2c.154tfs.278samples.lasso")

   load(system.file(package="TReNA", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
   target.gene <- "MEF2C"
   # print(fivenum(mtx.sub))   # 0.000000    1.753137   12.346965   43.247467 1027.765854

    print("note!  without log transform of the data")
    print("bayesSpike model is quite useless, even after")
    print("filtering for abs(beta) and pval")

   trena <- TReNA(mtx.assay=mtx.sub, solver="lasso", quiet=FALSE)
   tfs <- setdiff(rownames(mtx.sub), "MEF2C")
   tbl <- solve(trena, target.gene, tfs)
     # check for expected non-sensical values
   checkTrue(min(tbl$beta) < -7)
   checkTrue(max(tbl$beta) > 10)

      # with log transform, justified how?
      # good results are returned, as loosely checked
      # by correlating betas against  expression

   mtx.tmp <- mtx.sub - min(mtx.sub) + 0.001
   mtx.log2 <- log2(mtx.tmp)
   fivenum(mtx.log2)  # [1] -9.9657843  0.8107618  3.6262014  5.4345771 10.0052973

   trena <- TReNA(mtx.assay=mtx.log2, solver="lasso", quiet=FALSE)
   tfs <- setdiff(rownames(mtx.log2), "MEF2C")
   tbl2 <- solve(trena, target.gene, tfs)
   checkTrue(min(tbl2$beta) > -0.2)
   checkTrue(max(tbl2$beta) < 1)
   checkTrue(c("SATB2") %in% rownames(subset(tbl2, abs(beta) > 0.15)))

      # now transform mtx.sub with asinh

   mtx.asinh <- asinh(mtx.sub)
   fivenum(mtx.asinh)  # [1] 0.000000 1.327453 3.208193 4.460219 7.628290

   trena <- TReNA(mtx.assay=mtx.asinh, solver="lasso", quiet=FALSE)
   tfs <- setdiff(rownames(mtx.asinh), "MEF2C")
   tbl3 <- solve(trena, target.gene, tfs)
   checkTrue(min(tbl3$beta) > -0.2)
   checkTrue(max(tbl3$beta) < 1)
   checkTrue(c("SATB2") %in% rownames(subset(tbl2, abs(beta) > 0.15)))

   
} # test_ampAD.mef2c.154tfs.278samples.lasso
#----------------------------------------------------------------------------------------------------
test_ampAD.mef2c.154tfs.278samples.bayesSpike <- function()
{
   printf("--- test_ampAD.mef2c.154tfs.278samples.bayesSpike")

   load(system.file(package="TReNA", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
   target.gene <- "MEF2C"
   # print(fivenum(mtx.sub))   # 0.000000    1.753137   12.346965   43.247467 1027.765854

    print("note!  without log transform of the data")
    print("bayesSpike model is quite useless, even after")
    print("filtering for abs(beta) and pval")

   trena <- TReNA(mtx.assay=mtx.sub, solver="bayesSpike", quiet=FALSE)
   tfs <- setdiff(rownames(mtx.sub), "MEF2C")
   tbl <- solve(trena, target.gene, tfs)
   tbl.trimmed <- subset(tbl, abs(beta) > 0.1 & pval < 0.01)
   betas <- tbl.trimmed$beta
   big.abs.betas <- betas[abs(betas) > 1]
   checkTrue(length(big.abs.betas) > 20)

   checkTrue(nrow(tbl) > 10)
   checkTrue(cor(tbl.trimmed$beta, tbl.trimmed$gene.cor) < 0.2)

      # with log transform, justified how?
      # good results are returned, as loosely checked
      # by correlating betas against  expression

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
   
} # test_ampAD.mef2c.154tfs.278samples.bayesSpike
#----------------------------------------------------------------------------------------------------
test_ampAD.mef2c.154tfs.278samples.randomForest <- function()
{
   printf("--- test_ampAD.mef2c.154tfs.278samples.randomForest")

   load(system.file(package="TReNA", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
   target.gene <- "MEF2C"
   # print(fivenum(mtx.sub))   # 0.000000    1.753137   12.346965   43.247467 1027.765854

    print("note!  without log transform of the data")
    print("bayesSpike model is quite useless, even after")
    print("filtering for abs(beta) and pval")

   trena <- TReNA(mtx.assay=mtx.sub, solver="randomForest", quiet=FALSE)
   tfs <- setdiff(rownames(mtx.sub), "MEF2C")
   rf.result <- solve(trena, target.gene, tfs)
   tbl.scores <- rf.result$edges

   tbl.scores <- tbl.scores[order(tbl.scores$IncNodePurity, decreasing=TRUE),, drop=FALSE]

     # a loose test, ignoring rank of these 7 genes for now
   actual.genes.reported <- sort(rownames(subset(tbl.scores, IncNodePurity > 100000)))
   expected.genes <- sort(c("HLF", "STAT4", "SATB1", "SATB2", "FOXP2", "FOXO4","ATF2"))
   printf("1: expected: %s", paste(expected.genes, collapse=","))
   printf("1: actual: %s", paste(actual.genes.reported, collapse=","))
#   checkEquals(actual.genes.reported, expected.genes)

      # with log transform, justified how?
      # good results are returned, as loosely checked
      # by correlating betas against  expression

   mtx.tmp <- mtx.sub - min(mtx.sub) + 0.001
   mtx.log2 <- log2(mtx.tmp)
   fivenum(mtx.log2)  # [1] -9.9657843  0.8107618  3.6262014  5.4345771 10.0052973

   trena <- TReNA(mtx.assay=mtx.log2, solver="randomForest", quiet=FALSE)
   tfs <- setdiff(rownames(mtx.log2), "MEF2C")
   rf.result.2 <- solve(trena, target.gene, tfs)
   tbl.scores.2 <- rf.result.2$edges
   tbl.scores.2 <- tbl.scores.2[order(tbl.scores.2$IncNodePurity, decreasing=TRUE),, drop=FALSE]

     # a loose test, ignoring rank of these 7 genes for now
   actual.genes.reported <- sort(rownames(subset(tbl.scores.2, IncNodePurity > 100000)))
   expected.genes <- sort(c("HLF", "STAT4", "SATB1", "SATB2", "FOXP2", "DRGX","ATF2"))
   printf("2: expected: %s", paste(expected.genes, collapse=","))
   printf("2: actual: %s", paste(actual.genes.reported, collapse=","))
#   checkTrue( length( intersect(actual.genes.reported, expected.genes)) > 3 )

       # lasso reports, with log2 transformed data,
       # rownames(subset(tbl2, abs(beta) > 0.15)) "CUX1"   "FOXK2"  "SATB2"  "HLF"    "STAT5B" "ATF2"
   
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

      # with log transform, justified how?
      # good results are returned, as loosely checked
      # by correlating betas against  expression

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
test_LCLs.build_genomewide_model.lasso <- function()
{
   printf("--- test_LCLs.build_genomewide_model.lasso")


   load(system.file(package="TReNA", "extdata/lcl.13847genes.448samples.775TFsTFbs.76TFsChip.59TFsShRNA.RData"))

   expr = as.matrix(expr)

   gene.mean = rowMeans(expr)
   gene.sd = apply( expr , 1 , sd )
   enorm = ( expr - gene.mean ) / gene.sd
   checkTrue( median(apply( enorm , 1 , sd )) == 1 )

   tfbs = tfbs.sub[ , intersect( colnames(tfbs.sub) , rownames(expr) ) ]
   checkTrue( all( colnames(tfbs) %in% rownames(expr) ) )

   require( doParallel )
   ncores = detectCores()
   registerDoParallel( cores = floor(ncores/3) )

   checkTrue( nrow(expr) == 13847 )
   checkTrue( all(rownames(tfbs) == rownames(enorm) ) )

   # define candidate TFs from the distribution of TFBSs
   TfbsCountsQuantile = apply( tfbs , 2 , quantile , probs = 0.9 )
   candidate_regulators = t( t(tfbs) > 0.1*TfbsCountsQuantile )
   #candidate_regulators = tfbs > 0
   #checkTrue(
   #   all( tfbs[candidate_regulators[,1],1] > TfbsCountsQuantile[1] ))
   checkTrue( median(rowSums(candidate_regulators)) > 30 &
	median(rowSums(candidate_regulators)) < 200 )

   # select an appropriate lambda by evaluating a subset of genes
   trena = TReNA(mtx.assay=enorm,solver="lasso")

   fit.cv =
   foreach( target.gene=sample(rownames(enorm),100) ) %dopar% {
      tfs = names(which(candidate_regulators[target.gene,]==T))
      fit = solve(trena,target.gene,tfs,extraArgs=list(alpha=1,keep.metrics=T))
   }

   lambda = do.call( c ,
      lapply(1:length(fit.cv), function(i) fit.cv[[i]]$lambda))
   lambda.median = median(lambda,na.rm=T)
   checkTrue( lambda.median > 0 & lambda.median < 1 )

   # fit the model for all genes using the median lambda from fit.cv

   fit2 =
   foreach( target.gene=rownames(enorm)[1:100] ) %dopar% {
      # tfs = names(which(candidate_regulators[target.gene,]==T))
      tfs = names(which(candidate_regulators[target.gene,]==T))
      fit = solve(trena,target.gene,tfs,
        extraArgs=list(alpha=1,lambda=lambda.median,keep.metrics=T))
      if( length(fit) > 0 ) {
         if( nrow(fit$mtx.beta) > 0 ) {
            fit$mtx.beta$target = target.gene
            fit$mtx.beta$tf = rownames(fit$mtx.beta)
         }
      }
      return( fit )
   }

   r2 = do.call( c ,
      lapply(1:length(fit2), function(i) fit2[[i]]$r2))
   n.nonzero = do.call( c ,
      lapply(1:length(fit2), function(i) nrow(fit2[[i]]$mtx.beta)))
   trn = do.call( rbind ,
      lapply(1:length(fit2), function(i) fit2[[i]]$mtx.beta))

   checkTrue( any( r2 > 0.25 ) )
   checkTrue( median(n.nonzero) > 3 & median(n.nonzero) < 100 )
   checkTrue( ncol(trn) == 5 )
   checkTrue( nrow(trn) == sum(n.nonzero) )

} # test_LCLs.build_genomewide_model.lasso
#----------------------------------------------------------------------------------------------------
test_LCLs.build_genomewide_model.randomForest <- function()
{
   printf("--- test_LCLs.build_genomewide_model.randomForest")


   print(load(system.file(package="TReNA", "extdata/lcl.13847genes.448samples.775TFsTFbs.76TFsChip.59TFsShRNA.RData")))

   expr = as.matrix(expr)

   gene.mean = rowMeans(expr)
   gene.sd = apply( expr , 1 , sd )
   enorm = ( expr - gene.mean ) / gene.sd
   checkTrue( median(apply( enorm , 1 , sd )) == 1 )

   tfbs = tfbs.sub[ , intersect( colnames(tfbs.sub) , rownames(expr) ) ]
   checkTrue( all( colnames(tfbs) %in% rownames(expr) ) )

   require( doParallel )
   ncores = detectCores()
   registerDoParallel( cores = floor(ncores/2) )

   checkTrue( nrow(expr) == 13847 )
   checkTrue( all(rownames(tfbs) == rownames(enorm) ) )

   # define candidate TFs from the distribution of TFBSs
   TfbsCountsQuantile = apply( tfbs , 2 , quantile , probs = 0.9 )
   #candidate_regulators = t( t(tfbs) > 0.1*TfbsCountsQuantile )
   candidate_regulators = tfbs > 0

   trena <- TReNA(mtx.assay=enorm, solver="randomForest", quiet=FALSE)
   trn0 =
   foreach( target.gene=rownames(enorm)[1:100] ) %dopar% {
      tfs <- names(which( candidate_regulators[ target.gene , ] == T ))
      result <- solve(trena, target.gene, tfs)
      if( is.null(result) == F ) {
        result$edges$target = target.gene
        result$edges$tf = rownames(result$edges)
      }
      return(result)
   }

   r2 = do.call( c ,
      lapply(1:length(trn0), function(i) trn0[[i]]$r2))
   trn = do.call( rbind ,
      lapply(1:length(trn0), function(i) trn0[[i]]$edges))

   checkTrue( any( r2 > 0.1 ) )
   checkTrue( ncol(trn) == 3 )
   checkTrue( nrow(trn) > length(r2) )


}  # test_LCLs.build_genomewide_model.randomForest
#----------------------------------------------------------------------------------------------------
test_rosmapAssayMatrixUnderDifferentTransformations <- function()
{
    printf("--- test_rosmapAssayMatrixUnderDifferentTransformations");
    mtx.file <- "/local/mrichard/github/BDDS/trenadb/src/coryADpresentationNov2016/rosmapHuge.RData" # Changed hardcoded location
    checkTrue(file.exists(mtx.file))
    load(mtx.file)
    mtx <- tbl[, -1]
    mtx <- as.matrix(tbl[, -1])
    rownames(mtx) <- tbl$ensembl_id
    rows.to.delete <- setdiff(1:nrow(mtx), grep("ENSG", rownames(mtx)))
    mtx <- mtx[-rows.to.delete,]
    fivenum(mtx) # [1] 0.000000e+00 0.000000e+00 0.000000e+00 3.486310e-01 5.148396e+05

    trena <- TReNA(mtx.assay=mtx, solver="lasso", quiet=FALSE)

    gene.variances <- sort(apply(mtx, 1, var), decreasing=TRUE)
    target.gene <- names(gene.variances[100])  # a high, but otherwise arbitrary choice)
    set.seed(17)
    candidates <- rownames(mtx)[sample(1:nrow(mtx), 100)]
      # make sure the target is not among the candidates
    candidates <- setdiff(candidates, target.gene)
    tbl.betas <- solve(trena, target.gene, candidates, extraArgs=list(alpha=1.0, lambda=NULL))

    checkTrue(max(tbl.betas$beta) > 500)

    mtx2 <- asinh(mtx)
    trena <- TReNA(mtx.assay=mtx2, solver="lasso", quiet=FALSE)
    tbl.betas <- solve(trena, target.gene, candidates, extraArgs=list(alpha=1.0, lambda=NULL))
    checkTrue(max(tbl.betas$beta) < 10)
    
    mtx3 <- asinh(mtx2)
    trena <- TReNA(mtx.assay=mtx3, solver="lasso", quiet=FALSE)
    tbl.betas <- solve(trena, target.gene, candidates, extraArgs=list(alpha=1.0, lambda=NULL))
    checkTrue(max(tbl.betas$beta) < 1)

    # Repeat with Random Forest as the solver
    trena <- TReNA(mtx.assay = mtx, solver = "randomForest", quiet = "FALSE")
    rf.result <- solve(trena, target.gene, candidates)
    printf(max(rf.result))

    trena <- TReNA(mtx.assay = mtx2, solver = "randomForest", quiet = "FALSE")
    rf.result <- solve(trena, target.gene, candidates)
    printf(max(rf.result))

    trena <- TReNA(mtx.assay = mtx3, solver = "randomForest", quiet = "FALSE")
    rf.result <- solve(trena, target.gene, candidates)
    printf(max(rf.result))

    # Repeat with Bayes Spike as the solver
    trena <- TReNA(mtx.assay = mtx, solver = "bayesSpike", quiet = "FALSE")    
    tbl <- solve(trena, target.gene, tfs)
    tbl.trimmed <- subset(tbl, abs(beta) > 0.1 & pval < 0.01)
    betas <- tbl.trimmed$beta
    big.abs.betas <- betas[abs(betas) > 1]
    printf(length(big.abs.betas))

    trena <- TReNA(mtx.assay = mtx2, solver = "bayesSpike", quiet = "FALSE")    
    tbl <- solve(trena, target.gene, tfs)
    tbl.trimmed <- subset(tbl, abs(beta) > 0.1 & pval < 0.01)
    betas <- tbl.trimmed$beta
    big.abs.betas <- betas[abs(betas) > 1]
    printf(length(big.abs.betas))

    trena <- TReNA(mtx.assay = mtx3, solver = "bayesSpike", quiet = "FALSE")    
    tbl <- solve(trena, target.gene, tfs)
    tbl.trimmed <- subset(tbl, abs(beta) > 0.1 & pval < 0.01)
    betas <- tbl.trimmed$beta
    big.abs.betas <- betas[abs(betas) > 1]
    printf(length(big.abs.betas))

    
} # test_rosmapAssayMatrixUnderDifferentTransformations
#----------------------------------------------------------------------------------------------------
if(!interactive()) runTests()