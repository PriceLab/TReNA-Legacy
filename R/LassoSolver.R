#----------------------------------------------------------------------------------------------------
#' An S4 class to represent a LASSO solver
#'
#' @import glmnet
#' @include Solver.R
#' @import methods
#' 
#' @name LassoSolver-class

.LassoSolver <- setClass ("LassoSolver", contains="Solver")
#----------------------------------------------------------------------------------------------------
#' Create a Solver class object using the LASSO solver
#'
#' @param mtx.assay An assay matrix of gene expression data
#' @param quiet A logical denoting whether or not the solver should print output
#' @return A Solver class object with LASSO as the solver
#'
#' @export
#' 
#' @examples
#' solver <- LassoSolver()

LassoSolver <- function(mtx.assay=matrix(), quiet=TRUE)
{
    obj <- .LassoSolver(Solver(mtx.assay=mtx.assay, quiet=quiet))

    obj

} # LassoSolver, the constructor
#----------------------------------------------------------------------------------------------------
#' Get Lasso Solver name
#'
#' @param obj An object of class LassoSolver
#' 
#' @return "LassoSolver"
#'
#' @export
#' 
#' @examples
#' solver <- LassoSolver()
#' getSolverName(solver)

setMethod("getSolverName", "LassoSolver",

  function (obj){
     return("LassoSolver")
     })

#----------------------------------------------------------------------------------------------------
#' Run the LASSO Solver
#'
#' @rdname solve.Lasso
#' @aliases run.LassoSolver solve.Lasso
#' 
#' @description Given a TReNA object with LASSO as the solver, use the \code{\link{glmnet}} function to estimate coefficients
#' for each transcription factor as a predictor of the target gene's expression level. 
#'
#' @usage
#' tbl.out <- solve(obj, target.gene, tfs, tf.weights, extraArgs)
#' 
#' @param obj An object of class TReNA with "lasso" as the solver string
#' @param target.gene A designated target gene that should be part of the mtx.assay data
#' @param tfs The designated set of transcription factors that could be associated with the target gene.
#' @param tf.weights A set of weights on the transcription factors (default = rep(1, length(tfs)))
#' @param extraArgs Modifiers to the LASSO solver
#'
#' @return A data frame containing the coefficients relating the target gene to each transcription factor, plus other fit parameters.
#'
#' @seealso \code{\link{glmnet}}
#'
#' @family solver methods
#'
#' @examples
#' # Load included Alzheimer's data, create a TReNA object with LASSO as solver, and solve
#' load(system.file(package="TReNA", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
#' trena <- TReNA(mtx.assay = mtx.sub, solver = "lasso")
#' target.gene <- "MEF2C"
#' tfs <- setdiff(rownames(mtx.sub), target.gene)
#' tbl <- solve(trena, target.gene, tfs)

setMethod("run", "LassoSolver",

  function (obj, target.gene, tfs, tf.weights=rep(1,length(tfs)), extraArgs=list()){

      # Check if target.gene is in the bottom 10% in mean expression; if so, send a warning
      if(rowMeans(obj@mtx.assay)[target.gene] < quantile(rowMeans(obj@mtx.assay), probs = 0.1)){          
          warning("Target gene mean expression is in the bottom 10% of all genes in the assay matrix")          
      }      

      # Run the ElasticNetSolver function with alpha = 0.9 as default
      alpha = 0.9
      lambda <- NULL
      keep.metrics = FALSE

      if("alpha" %in% names(extraArgs))
          alpha <- extraArgs[["alpha"]]

      if("lambda" %in% names(extraArgs))
          lambda <- extraArgs[["lambda"]]

      if("keep.metrics" %in% names(extraArgs))
          keep.metrics <- extraArgs[["keep.metrics"]]
      
      mtx.beta <- .elasticNetSolver(obj, target.gene, tfs, tf.weights, alpha, lambda, keep.metrics)

      return(mtx.beta)
     })


#----------------------------------------------------------------------------------------------------
#' Rescale LASSO Predictor Weights
#'
#' @aliases rescalePredictorWeights.LassoSolver
#'
#' @param obj An object of class LassoSolver
#' @param rawValue.min The minimum value of the raw expression values
#' @param rawValue.max The maximum value of the raw expression values
#' @param rawValues A matrix of raw expression values
#'
#' @export
#' 
#' @return A matrix of the raw values re-scaled using the minimum and maximum values

# lasso penalizes predictors on a scale of 1 (full weight) to infinity (zero weight)
# here we wish to support incoming rawValues scaled between a possibly theoretical
# rawValue.min and rawValue.max
# we have empirical evidence that <large but non-infinite number> functions as a full penalty
# without distorting the scale so much that even good rawValues get reduced to nothing
# which is what .Machine$double.xmax would do

setMethod("rescalePredictorWeights", "LassoSolver",

   function (obj, rawValue.min, rawValue.max, rawValues){
      1 - ((rawValues-rawValue.min)/(rawValue.max-rawValue.min))
      })

#----------------------------------------------------------------------------------------------------

