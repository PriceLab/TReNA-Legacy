#----------------------------------------------------------------------------------------------------
#' An S4 class to represent a Spearman solver
#'
#' @include Solver.R
#' @name SpearmanSolver-class
#' 

.SpearmanSolver <- setClass ("SpearmanSolver", contains = "Solver")
#----------------------------------------------------------------------------------------------------
#' Create a Solver class object using Spearman correlation coefficients as the solver
#'
#' @param mtx.assay An assay matrix of gene expression data
#'
#' @return A Solver class object with Spearman correlation coefficients as the solver
#'
#' @examples
#' solver <- SpearmanSolver()

SpearmanSolver <- function(mtx.assay = matrix(), quiet=TRUE)
{
    obj <- .SpearmanSolver(Solver(mtx.assay=mtx.assay, quiet=quiet))

    obj

} #SpearmanSolver, the constructor
#----------------------------------------------------------------------------------------------------
#' Get Spearman Solver Name
#'
#' @param obj An object of class SpearmanSolver
#' 
#' @return "SpearmanSolver"
#'
#' @examples
#' solver <- SpearmanSolver()
#' getSolverName(solver)

setMethod("getSolverName", "SpearmanSolver",

          function (obj){
              return("SpearmanSolver")
          })

#----------------------------------------------------------------------------------------------------
#' Run the Spearman Solver
#'
#' @aliases run.SpearmanSolver
#' @description Given a TReNA object with Spearman as the solver, use the \code{\link{cor}} function with
#' \code{method = "spearman"} to esimate coefficients for each transcription factor as a predictor of the target
#' gene's expression level
#'
#' @param obj An object of class SpearmanSolver
#' @param target.gene A designated target gene that should be part of the mtx.assay data
#' @param tfs The designated set of transcription factors that could be associated with the target gene.
#'
#' @return fit The set of Spearman Correlation Coefficients between each transcription factor and the target gene.
#'
#' @seealso \code{\link{cor}}

setMethod("run", "SpearmanSolver",

          function (obj, target.gene, tfs, tf.weights = rep(1,length(tfs), extraArgs=list())){

              mtx <- obj@mtx.assay
              # Check that target gene and tfs are all part of the matrix
              stopifnot(target.gene %in% rownames(mtx))
              stopifnot(all(tfs %in% rownames(mtx)))
              # If given no tfs, return nothing
              if (length(tfs)==0) return(NULL)

              # Don't handle tf self-regulation, so take target gene out of tfs
              deleters <- grep(target.gene, tfs)
              if(length(deleters) > 0){
                  tfs <- tfs[-deleters]
                  tf.weights = tf.weights[-deleters]
              }
              # If target gene was the only tf, then return nothing
              if(length(tfs)==0) return(NULL)

              x = t(mtx[tfs,,drop=F])
              y = as.vector(t(mtx[target.gene,])) # Make target gene levels into a vector

              # Calculate Spearman correlation coefficients
              fit <- cor( x = x, y = y, method = "spearman")

                                        # For now, just return the fit
              return(fit)
})
