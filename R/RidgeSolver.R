#------------------------------------------------------------------------------------------------------------------------
#' An S4 class to represent a Ridge Regression solver
#'
#' @include Solver.R
#' @name RidgeSolver-class

.RidgeSolver <- setClass ("RidgeSolver", contains="Solver")
#------------------------------------------------------------------------------------------------------------------------
#' Create a Solver class object using the Ridge Regression solver
#'
#' @param mtx.assay An assay matrix of gene expression data
#'
#' @return A Solver class object with Ridge Regression as the solver
#'
#' @examples
#' solver <- RidgeSolver()

RidgeSolver <- function(mtx.assay=matrix(), quiet=TRUE)
{
    obj <- .RidgeSolver(Solver(mtx.assay=mtx.assay, quiet=quiet))

    obj

} # RidgeSolver, the constructor
#------------------------------------------------------------------------------------------------------------------------
#' Get Ridge Solver name
#'
#' @param obj An object of class RidgeSolver
#' 
#' @return "RidgeSolver"
#'
#' @examples
#' solver <- RidgeSolver()
#' getSolverName(solver)

setMethod("getSolverName", "RidgeSolver",

  function (obj){
     return("RidgeSolver")
     })

#----------------------------------------------------------------------------------------------------
#' Run the Ridge Regression Solver
#' @aliases run.RidgeSolver
#' @description Given a TReNA object with Ridge Regression as the solver, use the \code{\link{glmnet}} function to estimate coefficients
#' for each transcription factor as a predictor of the target gene's expression level. 
#'
#' @param obj An object of class RidgeSolver
#' @param target.gene A designated target gene that should be part of the mtx.assay data
#' @param tfs The designated set of transcription factors that could be associated with the target gene.
#' @param tf.weights A set of weights on the transcription factors (default = rep(1, length(tfs)))
#' @param extraArgs Modifiers to the Ridge Regression solver
#'
#' @return A data frame containing the coefficients relating the target gene to each transcription factor, plus other fit parameters.
#'
#' @seealso \code{\link{glmnet}}
#'
#' @examples
#' 

setMethod("run", "RidgeSolver",

  function (obj, target.gene, tfs, tf.weights=rep(1,length(tfs)), extraArgs=list()){

      # Run LassoSolver, but use alpha = 0

      trena <- TReNA(mtx.assay = obj@solver@mtx.assay, solver = "lasso")
      tbl <- solve(trena, target.gene, tfs, extraArgs = list("alpha" = 0))

      return(tbl)
})
#----------------------------------------------------------------------------------------------------
#' Rescale Ridge Regression Predictor Weights
#'
#' @aliases rescalePredictorWeights.RidgeSolver
#'
#' @param obj An object of class RidgeSolver
#' @param rawValue.min The minimum value of the raw expression values
#' @param rawValue.max The maximum value of the raw expression values
#' @param rawValues A matrix of raw expression values
#'
#' @return A matrix of the raw values re-scaled using the minimum and maximum values

# lasso penalizes predictors on a scale of 1 (full weight) to infinity (zero weight)
# here we wish to support incoming rawValues scaled between a possibly theoretical
# rawValue.min and rawValue.max
# we have empirical evidence that <large but non-infinite number> functions as a full penalty
# without distorting the scale so much that even good rawValues get reduced to nothing
# which is what .Machine$double.xmax would do

setMethod("rescalePredictorWeights", "RidgeSolver",

   function (obj, rawValue.min, rawValue.max, rawValues){
      1 - ((rawValues-rawValue.min)/(rawValue.max-rawValue.min))
      })

#----------------------------------------------------------------------------------------------------

