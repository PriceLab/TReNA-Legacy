#----------------------------------------------------------------------------------------------------
#' Class RidgeSolver
#'
#' @include Solver.R
#' @import glmnet
#' @import methods
#' 
#' @name RidgeSolver-class
#' @rdname RidgeSolver-class

.RidgeSolver <- setClass ("RidgeSolver", contains="Solver")
#----------------------------------------------------------------------------------------------------
#' Create a Solver class object using the Ridge solver
#' 
#' @param mtx.assay An assay matrix of gene expression data
#' @param quiet A logical denoting whether or not the solver should print output
#'
#' @return A Solver class object with Ridge as the solver
#'
#' @seealso  \code{\link{solve.Ridge}}, \code{\link{getAssayData}}
#'
#' @family Solver class objects
#' 
#' @export
#' 
#' @examples
#' solver <- RidgeSolver()

RidgeSolver <- function(mtx.assay=matrix(), quiet=TRUE)
{
    obj <- .RidgeSolver(Solver(mtx.assay=mtx.assay, quiet=quiet))

    obj

} # RidgeSolver, the constructor
#----------------------------------------------------------------------------------------------------
#' Run the Ridge Regression Solver
#'
#' @rdname solve.Ridge
#' @aliases run.RidgeSolver solve.Ridge
#' 
#' @description Given a TReNA object with Ridge Regression as the solver,
#' use the \code{\link{glmnet}} function to estimate coefficients
#' for each transcription factor as a predictor of the target gene's expression level.
#' This method should be called using the \code{\link{solve}} method on an appropriate TReNA object.
#'
#' @param obj An object of class Solver with "ridge" as the solver string
#' @param target.gene A designated target gene that should be part of the mtx.assay data
#' @param tfs The designated set of transcription factors that could be associated with the target gene.
#' @param tf.weights A set of weights on the transcription factors (default = rep(1, length(tfs)))
#' @param extraArgs Modifiers to the Ridge Regression solver
#'
#' @return A data frame containing the coefficients relating the target gene to each
#' transcription factor, plus other fit parameters.
#'
#' @seealso \code{\link{glmnet}}, , \code{\link{RidgeSolver}}
#'
#' @family solver methods
#' 
#' @examples
#' # Load included Alzheimer's data, create a TReNA object with Bayes Spike as solver, and solve
#' load(system.file(package="TReNA", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
#' trena <- TReNA(mtx.assay = mtx.sub, solver = "ridge")
#' target.gene <- "MEF2C"
#' tfs <- setdiff(rownames(mtx.sub), target.gene)
#' tbl <- solve(trena, target.gene, tfs)

setMethod("run", "RidgeSolver",

  function (obj, target.gene, tfs, tf.weights=rep(1,length(tfs)), extraArgs=list()){

      # Check if target.gene is in the bottom 10% in mean expression; if so, send a warning      
      if(rowMeans(getAssayData(obj))[target.gene] < stats::quantile(rowMeans(getAssayData(obj)), probs = 0.1)){          
          warning("Target gene mean expression is in the bottom 10% of all genes in the assay matrix")          
      }      

      # Run Elastic Net, but use alpha = 0
      alpha = 0
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

