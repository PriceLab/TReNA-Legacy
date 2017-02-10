# Extra functions trainModel and predictFromModel

# Generics set in TReNA.R

----------------------------------------------------------------------------------------------------
setGeneric("trainModel",               signature="obj", function(obj, target.gene, tfs, training.samples,
                                                                 tf.weights=rep(1, length(tfs)))
                                                           standardGeneric ("trainModel"))
setGeneric("predictFromModel",         signature="obj", function(obj, model, tfs, test.samples)
                                                           standardGeneric ("predictFromModel"))
#----------------------------------------------------------------------------------------------------
setMethod("trainModel", "TReNA",

   function (obj, target.gene, tfs, training.samples, tf.weights=rep(1, length(tfs))){
      fit <- trainModel(obj@solver, target.gene, tfs, training.samples, tf.weights)
      return(fit)
      })

#----------------------------------------------------------------------------------------------------
setMethod("predictFromModel", "TReNA",

   function (obj, model, tfs, test.samples){
      prediction <- predictFromModel(obj@solver, model, tfs, test.samples)
      return(prediction)
      })

#------------------------------------------------------------------------------------------------------------------------
# Inherited methods in BayesSpikeSolver.R

#' Train a model using Bayes Spike
#'
#' @rdname trainModel.BayesSpikeSolver
#' @aliases trainModel.BayesSpikeSolver
#' @param target.gene A designated target gene that should be part of the mtx.assay data
#' @param tfs The designated set of transcription factors that could be associated with the target gene.
#' @param training.samples ***These aren't currently used!
#' @param tf.weights A set of weights on the transcription factors (default = rep(1, length(tfs)))
#'
#' @return A model relating the target gene (response) to the transcription factors (predictors)

setGeneric("trainModel", function(obj, target.gene, tfs, training.samples, tf.weights=rep(1, length(tfs))) standardGeneric("trainModel"))
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
#' Predict expression using a Bayes Spike Model
#'
#' @aliases predictFromModel.BayesSpikeSolver
#' @param model A Bayes Spike model generated using \code{\\link{trainModel.BayesSpikeSolver}} that relates the transcription
#' factors to the target gene of interest
#' @param tfs A list of transcription factor predictors matching those used in the Bayes Spike model
#' @param test.samples A matrix of gene expression levels for the supplied transcription factors
#'
#' @return Predictions of the target gene's expression level for each sample in the test sample set
#' @seealso \code{\\link{trainModel.BayesSpikeSolver}}

setGeneric("predictFromModel", function(obj, model, tfs, test.samples) standardGeneric("predictFromModel"))
setMethod("predictFromModel", "BayesSpikeSolver",

   function (obj, model, tfs, test.samples){
      mtx <- obj@mtx.assay
      predict.glmnet(model, t(mtx[tfs, test.samples]))
      })

#----------------------------------------------------------------------------------------------------
# Inherited methods in LassoSolver.R
#' Train a model using LASSO
#'
#' @aliases trainModel.LassoSolver
#' @param target.gene A designated target gene that should be part of the mtx.assay data
#' @param tfs The designated set of transcription factors that could be associated with the target gene.
#' @param training.samples ***These aren't currently used!
#' @param tf.weights A set of weights on the transcription factors (default = rep(1, length(tfs)))
#'
#' @return A model relating the target gene (response) to the transcription factors (predictors)

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
#' Predict expression using a LASSO Model
#'
#' @aliases predictFromModel.LassoSolver
#' @param model A LASSO model generated using \code{\\link{trainModel}} that relates the transcription
#' factors to the target gene of interest
#' @param tfs A list of transcription factor predictors matching those used in the LASSO model
#' @param test.samples A matrix of gene expression levels for the supplied transcription factors
#'
#' @return Predictions of the target gene's expression level for each sample in the test sample set
#' @seealso \code{\\link{trainModel}}

setMethod("predictFromModel", "LassoSolver",

   function (obj, model, tfs, test.samples){
      mtx <- obj@mtx.assay
      predict.glmnet(model, t(mtx[tfs, test.samples]))
      })

#----------------------------------------------------------------------------------------------------

