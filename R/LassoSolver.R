#------------------------------------------------------------------------------------------------------------------------
.LassoSolver <- setClass ("LassoSolver", contains="Solver")
#------------------------------------------------------------------------------------------------------------------------
LassoSolver <- function(mtx.assay=matrix(), quiet=TRUE)
{
    obj <- .LassoSolver(Solver(mtx.assay=mtx.assay, quiet=quiet))

    obj

} # TReNA, the constructor
#------------------------------------------------------------------------------------------------------------------------
setMethod("getSolverName", "LassoSolver",

  function (obj){
     return("LassoSolver")
     })

#----------------------------------------------------------------------------------------------------
setMethod("run", "LassoSolver",

  function (obj, target.gene, tfs, tf.weights=rep(1,length(tfs))){

     tf.weights <- 1/tf.weights

     mtx <- obj@mtx.assay
     stopifnot(target.gene %in% rownames(mtx))
     stopifnot(all(tfs %in% rownames(mtx)))
     features <- t(mtx[tfs, ])
     target <- as.numeric(mtx[target.gene,])

     if(!obj@quiet)
         printf("begining cross-validation for glmnet, using %d tfs, target %s", length(tfs), target.gene)

     cv.out <- cv.glmnet(features, target, penalty.factor=tf.weights, grouped=FALSE)
     lambda.min <- cv.out$lambda.min
     lambda.1se <-cv.out$lambda.1se

     fit.nolambda = glmnet(features, target, penalty.factor=tf.weights)
     fit = glmnet(features, target, penalty.factor=tf.weights, lambda=lambda.1se)

       # extract the exponents of the fit
     #tbl.out <- as.matrix(fit$beta)
     #deleters <- as.integer(which(tbl.out[,1] == 0))
     #if(length(deleters) > 0)
     #   tbl.out <- tbl.out[-deleters, , drop=FALSE]
     #colnames(tbl.out) <- "beta"

     mtx.beta <- as.matrix(coef(fit, s="lambda.min"))
     colnames(mtx.beta) <- "beta"
     deleters <- as.integer(which(mtx.beta[,1] == 0))
     if(length(deleters) > 0)
        mtx.beta <- mtx.beta[-deleters, , drop=FALSE]

        # put the intercept, admittedly with much redundancy, into its owh column
     intercept <- mtx.beta[1,1]
     mtx.beta <- mtx.beta[-1, , drop=FALSE]
     mtx.beta <- cbind(mtx.beta, intercept=rep(intercept, nrow(mtx.beta)))
     correlations.of.betas.to.targetGene <- unlist(lapply(rownames(mtx.beta), function(x) cor(mtx[x,], mtx[target.gene,])))
     #browser()
     mtx.beta <- cbind(mtx.beta, gene.cor=correlations.of.betas.to.targetGene)
     if(!obj@quiet)
        plot(fit.nolambda, xvar='lambda', label=TRUE)

     return(mtx.beta)
     })


#----------------------------------------------------------------------------------------------------
setMethod("trainModel", "LassoSolver",

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
setMethod("predictFromModel", "LassoSolver",

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
setMethod("rescalePredictorWeights", "LassoSolver",

   function (obj, rawValue.min, rawValue.max, rawValues){
      1 - ((rawValues-rawValue.min)/(rawValue.max-rawValue.min))
      })

#----------------------------------------------------------------------------------------------------

