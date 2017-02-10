#------------------------------------------------------------------------------------------------------------------------
#' An S4 class to represent a Bayes Spike solver
#'
#' @include Solver.R
#' @name BayesSpikeSolver-class

.BayesSpikeSolver <- setClass ("BayesSpikeSolver", contains="Solver")
#------------------------------------------------------------------------------------------------------------------------
#' Designate Bayes Spike as the TReNA Solver and Solve
#'
#' @param mtx.assay An assay matrix of gene expression data
#'
#' @return A Solver class object with Bayes Spike as the solver

BayesSpikeSolver <- function(mtx.assay=matrix(), quiet=TRUE)
{
    obj <- .BayesSpikeSolver(Solver(mtx.assay=mtx.assay, quiet=quiet))

    obj

} # TReNA, the constructor
#------------------------------------------------------------------------------------------------------------------------
#' Get Bayes Spike Solver name
#'
#' @param obj An object of the class BayesSpikeSolver
#' 
#' @return "BayesSpikeSolver"
#' 

setMethod("getSolverName", "BayesSpikeSolver",

  function (obj){
     return("BayesSpikeSolver")
     })

#----------------------------------------------------------------------------------------------------
#' Run the Bayes Spike Solver
#'
#' @aliases run.BayesSpikeSolver
#' @description Given a TReNA object with Bayes Spike as the solver, use the \code{\link{vbsr}} function to estimate coefficients
#' for each transcription factor as a predictor of the target gene's expression level.
#'
#' @param obj An object of the class BayesSpikeSolver
#' @param target.gene A designated target gene that should be part of the mtx.assay data
#' @param tfs The designated set of transcription factors that could be associated with the target gene.
#' @param tf.weights A set of weights on the transcription factors (default = rep(1, length(tfs)))
#' @param extraArgs Modifiers to the Bayes Spike solver
#'
#' @return A data frame containing the coefficients relating the target gene to each transcription factor, plus other fit parameters.
#'
#' @seealso \code{\link{vbsr}}

setMethod("run", "BayesSpikeSolver",

  function (obj, target.gene, tfs, tf.weights=rep(1,length(tfs), extraArgs=list())){

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
#' Rescale Bayes Spike Predictor Weights
#'
#' @aliases rescalePredictorWeights.BayesSpikeSolver
#'
#' @param obj An object of class BayesSpikeSolver
#' @param rawValue.min The minimum value of the raw expression values
#' @param rawValue.max The maximum value of the raw expression values
#' @param rawValues A matrix of raw expression values
#'
#' @return A matrix of the raw values re-scaled using the minimum and maximum values

# lasso penalizes predictors on a scale of 1 (full weight) to infinity (zero weight)
# hwere we wish to support incoming rawValues scaled between a possibly theoretical
# rawValue.min and rawValue.max
# we have empirical evidence that <large but non-infinite number> functions as a full penalty
# without distorting the scale so much that even good rawValues get reduced to nothing
# which is what .Machine$double.xmax would do

setGeneric("rescalePredictorWeights", function(obj, rawValue.min, rawValue.max, rawValues) standardGeneric("rescalePredictorWeights"))
setMethod("rescalePredictorWeights", "BayesSpikeSolver",

   function (obj, rawValue.min, rawValue.max, rawValues){
      rawValues
      })

#----------------------------------------------------------------------------------------------------

