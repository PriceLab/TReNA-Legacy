#----------------------------------------------------------------------------------------------------
#' An S4 class to represent a Pearson solver
#'
#' @include Solver.R
#' @name PearsonSolver-class
#' 

.PearsonSolver <- setClass ("PearsonSolver", contains = "Solver")
#----------------------------------------------------------------------------------------------------
#' Create a Solver class object using  Pearson correlation coefficients as the solver
#'
#' @param mtx.assay An assay matrix of gene expression data
#'
#' @return A Solver class object with Pearson correlation coefficients as the solver
#'
#' @examples
#' solver <- PearsonSolver()

PearsonSolver <- function(mtx.assay = matrix(), quiet=TRUE)
{
    obj <- .PearsonSolver(Solver(mtx.assay=mtx.assay, quiet=quiet))

    obj

} #PearsonSolver, the constructor
#----------------------------------------------------------------------------------------------------
#' Get Pearson Solver name
#'
#' @param obj An object of class PearsonSolver
#'
#' @return "PearsonSolver
#'
#' @examples
#' solver <- PearsonSolver()
#' getSolverName(solver)

setMethod("getSolverName", "PearsonSolver",

          function (obj){
              return("PearsonSolver")
          })
#----------------------------------------------------------------------------------------------------
#' Run the Pearson Solver
#' 
#' @aliases run.PearsonSolver
#' @description Given a TReNA object with Pearson as the solver, use the \code{\link{cor}} function to
#' estimate coefficients for each transcription factor as a perdictor of the target gene's expression level
#'
#' @param obj An object of class PearsonSolver
#' @param target.gene A designated target gene that should be part of the mtx.assay data
#' @param tfs The designated set of transcription factors that could be associated with the target gene.
#'
#' @return fit The set of Pearson Correlation Coefficients between each transcription factor and the target gene.
#'
#' @seealso \code{\link{cor}}

setMethod("run", "PearsonSolver",

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

              # Calculate Pearson correlation coefficients
              fit <- cor( x = x, y = y)

              # Return the coefficients as a data frame 
              tbl <- data.frame(row.names = rownames(fit)[order(abs(fit), decreasing = TRUE)],
                                coefficient = fit[order(abs(fit), decreasing = TRUE)])

              return(tbl)
          })
#----------------------------------------------------------------------------------------------------
