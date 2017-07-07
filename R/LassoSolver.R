#' Class LassoSolver
#'
#' @import glmnet
#' @include Solver.R
#' @import methods
#' 
#' @name LassoSolver-class

.LassoSolver <- setClass ("LassoSolver",
                          contains="Solver",
                          slots = c(regulatorWeights="numeric",
                                    alpha = "numeric",
                                    lambda = "numeric",
                                    keep.metrics = "logical")
                          )
#----------------------------------------------------------------------------------------------------
#' Create a Solver class object using the LASSO solver
#'
#' @param mtx.assay An assay matrix of gene expression data
#' @param quiet A logical denoting whether or not the solver should print output
#'
#' @return A Solver class object with LASSO as the solver
#'
#' @seealso  \code{\link{solve.Lasso}}, \code{\link{getAssayData}}
#'
#' @family Solver class objects
#' 
#' @export
#' 
#' @examples
#' solver <- LassoSolver()

LassoSolver <- function(mtx.assay=matrix(), targetGene, candidateRegulators,
                        regulatorWeights=rep(1, length(candidateRegulators)),
                        alpha = 0.9, lambda = numeric(0),
                        keep.metrics = FALSE, quiet=TRUE)
{
   if(any(grepl(targetGene, candidateRegulators)))
      candidateRegulators <- candidateRegulators[-grep(targetGene, candidateRegulators)]

   candidateRegulators <- intersect(candidateRegulators, rownames(mtx.assay))
   stopifnot(length(candidateRegulators) > 0)

   obj <- .LassoSolver(mtx.assay=mtx.assay,
                              quiet=quiet,
                              targetGene=targetGene,
                              candidateRegulators=candidateRegulators,
                              regulatorWeights=regulatorWeights,
                              alpha = alpha,
                              lambda = lambda,
                              keep.metrics = keep.metrics
                       )
    obj

} # LassoSolver, the constructor
#----------------------------------------------------------------------------------------------------
#' Run the LASSO Solver
#'
#' @rdname solve.Lasso
#' @aliases run.LassoSolver solve.Lasso
#' 
#' @description Given a TReNA object with LASSO as the solver, use the \code{\link{glmnet}} function
#' to estimate coefficients for each transcription factor as a predictor of the target gene's
#' expression level. This method should be called using the \code{\link{solve}} method on an
#' appropriate TReNA object.
#'
#' @param obj An object of class Solver with "lasso" as the solver string

#' @param tfs The designated set of transcription factors that could be associated with the target gene.
#' @param tf.weights A set of weights on the transcription factors (default = rep(1, length(tfs)))
#' @param extraArgs Modifiers to the LASSO solver
#'
#' @return A data frame containing the coefficients relating the target gene to each transcription factor, plus other fit parameters.
#'
#' @seealso \code{\link{glmnet}},, \code{\link{LassoSolver}}
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

  function (obj){

      mtx <- getAssayData(obj)
      target.gene <- getTarget(obj)
      tfs <- getRegulators(obj)
      
      # Check if target.gene is in the bottom 10% in mean expression; if so, send a warning
      if(rowMeans(mtx)[target.gene] < stats::quantile(rowMeans(mtx), probs = 0.1)){          
          warning("Target gene mean expression is in the bottom 10% of all genes in the assay matrix")          
      }      
      
      mtx.beta <- .elasticNetSolver(obj, target.gene, tfs,
                                    obj@regulatorWeights,
                                    obj@alpha,
                                    obj@lambda,
                                    obj@keep.metrics)

      return(mtx.beta)
     })
#----------------------------------------------------------------------------------------------------

