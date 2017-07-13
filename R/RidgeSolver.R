#----------------------------------------------------------------------------------------------------
#' Class RidgeSolver
#'
#' @include Solver.R
#' @import glmnet
#' @import methods
#' 
#' @name RidgeSolver-class
#' @rdname RidgeSolver-class

.RidgeSolver <- setClass ("RidgeSolver",
                          contains="Solver",
                          slots = c(regulatorWeights = "numeric",
                                    alpha = "numeric",
                                    lambda = "numeric",
                                    keep.metrics = "logical")
                          )
#----------------------------------------------------------------------------------------------------
#' Create a Solver class object using the Ridge solver
#' 
#' @param mtx.assay An assay matrix of gene expression data
#' @param target.gene A designated target gene that should be part of the mtx.assay data
#' @param candidateRegulators The designated set of transcription factors that could be associated
#' with the target gene
#' @param regulatorWeights A set of weights on the transcription factors
#' (default = rep(1, length(tfs)))
#' @param alpha A parameter from 0-1 that determines the proportion of LASSO to ridge used in the
#' elastic net solver, with 0 being fully ridge and 1 being fully LASSO (default = 0.9)
#' @param lambda A tuning parameter that determines the severity of the penalty function imposed
#' on the elastic net regression. If unspecified, lambda will be determined via
#' permutation testing (default = numeric(0)).
#' @param keep.metrics A logical denoting whether or not to keep the initial supplied metrics
#' versus determining new ones
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

RidgeSolver <- function(mtx.assay=matrix(), targetGene, candidateRegulators,
                        regulatorWeights = rep(1, length(candidateRegulators)),
                        alpha = 0, lambda = numeric(0),
                        keep.metrics = FALSE, quiet=TRUE)
{
   if(any(grepl(targetGene, candidateRegulators)))
      candidateRegulators <- candidateRegulators[-grep(targetGene, candidateRegulators)]

   candidateRegulators <- intersect(candidateRegulators, rownames(mtx.assay))
   stopifnot(length(candidateRegulators) > 0)

   obj <- .RidgeSolver(mtx.assay=mtx.assay,
                              quiet=quiet,
                              targetGene=targetGene,
                              candidateRegulators=candidateRegulators,
                              regulatorWeights=regulatorWeights,
                              alpha = alpha,
                              lambda = lambda,
                              keep.metrics = keep.metrics
                       )
    obj

} # RidgeSolver, the constructor
#----------------------------------------------------------------------------------------------------
setMethod('show', 'RidgeSolver',

    function(obj) {
       regulator.count <- length(getRegulators(obj))
       if(regulator.count > 10){
          regulatorString <- paste(getRegulators(obj)[1:10], collapse=",")
          regulatorString <- sprintf("%s...", regulatorString);
          }
       else
          regulatorString <- paste(getRegulators(obj), collapse=",")

       msg = sprintf("RidgeSolver with mtx.assay (%d, %d), targetGene %s, %d candidate regulators %s, alpha = %f",
                     nrow(getAssayData(obj)), ncol(getAssayData(obj)),
                     getTarget(obj), regulator.count, regulatorString, obj@alpha)
       cat (msg, '\n', sep='')
    })
#----------------------------------------------------------------------------------------------------
#' Run the Ridge Regression Solver
#'
#' @rdname solve.Ridge
#' @aliases run.RidgeSolver solve.Ridge
#' 
#' @description Given a TReNA object with Ridge Regression as the solver,
#' use the \code{\link{glmnet}} function to estimate coefficients
#' for each transcription factor as a predictor of the target gene's expression level.
#'
#' @param obj An object of class RidgeSolver
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
#' target.gene <- "MEF2C"
#' tfs <- setdiff(rownames(mtx.sub), target.gene)
#' ridge.solver <- LassoSolver(mtx.sub, target.gene, tfs)
#' tbl <- run(ridge.solver)

setMethod("run", "RidgeSolver",

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

