library(TReNA)
library(RUnit)
#----------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
    test_fitDREAM5_yeast.lasso()
    test_fitDREAM5_yeast.lasso_weighted.tfs()
    test_fitDREAM5_yeast.randomForest()
    test_fitDREAM5_yeast.bayesSpike()
    test_trainAndPredict_DREAM5_yeast.lasso()
    
    test_ampAD.mef2c.154tfs.278samples.bayesSpike.nonCodingGenes()
    test_LCLs.build_genomewide_model.lasso()
    test_LCLs.build_genomewide_model.randomForest()
    

} # runTests
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
# FootprintFinder originally accepted only HUGO gene symbols, which is still the expected case
# however, ensembl reports, and we currently have expression data for, a variety of DNA elements
# including miRNA, linkRNA, antisense genes, pseudogenes of various sorts, etc.
# these each have a unique ENSG id, which we test out here
# the constructor of FootprintFinder needs to recognise these identifiers, and make a corresponding
# call to biomart
test_getFootprintsForEnsemblGenes <- function()
{
   printf("--- test_getFootprintsForEnsemblGenes")

   genome.db.uri <- "postgres://whovian/hg38"
   project.db.uri <-  "postgres://whovian/brain_wellington"
   fp <- FootprintFinder(genome.db.uri, project.db.uri, quiet=TRUE)

   genes <- c("ENSG00000267051", "ENSG00000264503", "ENSG00000273141", "ENSG00000212712",
              "ENSG00000236396", "ENSG00000154889", "ENSG00000267794",  "ENSG00000264843",
              "ENSG00000260759", "ENSG00000154856")

   goi <- genes[3]
   loc <- getGenePromoterRegion(fp, goi, 250, 0)
   tbl <- getFootprintsForGene(fp, goi, 250, 0)
   checkTrue(all(tbl$mfpStart >= loc$start))
   checkTrue(all(tbl$mfpStart <= loc$end))
   checkTrue(all(tbl$mfpEnd >= loc$start))
   checkTrue(all(tbl$mfpEnd <= loc$end))

} # test_getFootprintsForEnsemblGenes
#----------------------------------------------------------------------------------------------------
if(!interactive()) runTests()
