#------------------------------------------------------------------------------------------------------------------------
.BayesSpikeSolver <- setClass ("BayesSpikeSolver", contains="Solver")
#------------------------------------------------------------------------------------------------------------------------
BayesSpikeSolver <- function(mtx.assay=matrix(), quiet=TRUE)
{
    obj <- .BayesSpikeSolver(Solver(mtx.assay=mtx.assay, quiet=quiet))

    obj

} # TReNA, the constructor
#------------------------------------------------------------------------------------------------------------------------
setMethod("getSolverName", "BayesSpikeSolver",

  function (obj){
     return("BayesSpikeSolver")
     })

#----------------------------------------------------------------------------------------------------
setMethod("run", "BayesSpikeSolver",

  function (obj, target.gene, tfs, tf.weights=rep(1,length(tfs))){

        # we don't try to handle tf self-regulation
    deleters <- grep(target.gene, tfs)
    if(length(deleters) > 0){
       tfs <- tfs[-deleters]
       tf.weights <- tf.weights[-deleters]
       message(sprintf("BayesSpikeSolver removing target.gene from candidate regulators: %s", target.gene))
       }

    mtx <- obj@mtx.assay
    stopifnot(target.gene %in% rownames(mtx))
    stopifnot(all(tfs %in% rownames(mtx)))
    features <- t(mtx[tfs, ])
    target <- as.numeric(mtx[target.gene,])
    result <- vbsr(target, features, family='normal')
    tbl.out <- data.frame(beta=result$beta, pval=result$pval, z=result$z, post=result$post)
    rownames(tbl.out) <- tfs
    tbl.out$score <- -log10(tbl.out$pval)
    tbl.out <- tbl.out[order(tbl.out$score, decreasing=TRUE),]
    #browser()
    gene.cor <- sapply(rownames(tbl.out), function(tf) cor(mtx[tf,], mtx[target.gene,]))
    tbl.out$gene.cor <- as.numeric(gene.cor)
    tbl.out
    })


#----------------------------------------------------------------------------------------------------
setMethod("trainModel", "BayesSpikeSolver",

   function (obj, target.gene, tfs, training.samples, tf.weights=rep(1, length(tfs))){

       # transform to the 1..infinity range penalty.factor expects,
       # where infinity indicates "exclude this tf"

     tf.weights <- 1/tf.weights

     mtx <- obj@mtx.assay
     stopifnot(target.gene %in% rownames(mtx))
     stopifnot(all(tfs %in% rownames(mtx)))
     features <- t(mtx[tfs, ])
     target <- t(mtx[target.gene, , drop=FALSE])

     cv.out <- cv.glmnet(features, target, penalty.factor=tf.weights, grouped=FALSE)
     lambda.min <- cv.out$lambda.min
     glmnet(features, target, penalty.factor=tf.weights, lambda=lambda.min)
     })

#----------------------------------------------------------------------------------------------------
setMethod("predictFromModel", "BayesSpikeSolver",

   function (obj, model, tfs, test.samples){
      mtx <- obj@mtx.assay
      predict.glmnet(model, t(mtx[tfs, test.samples]))
      })

#----------------------------------------------------------------------------------------------------
# lasso penalizes predictors on a scale of 1 (full weight) to infinity (zero weight)
# hwere we wish to support incoming rawValues scaled between a possibly theoretical
# rawValue.min and rawValue.max
# we have empirical evidence that <large but non-infinite number> functions as a full penalty
# without distorting the scale so much that even good rawValues get reduced to nothing
# which is what .Machine$double.xmax would do
setMethod("rescalePredictorWeights", "BayesSpikeSolver",

   function (obj, rawValue.min, rawValue.max, rawValues){
      rawValues
      })

#----------------------------------------------------------------------------------------------------

