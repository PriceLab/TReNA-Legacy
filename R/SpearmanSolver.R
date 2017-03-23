#----------------------------------------------------------------------------------------------------
#' An S4 class to represent a Spearman solver
#'
#' @include Solver.R
#'
#' @import methods
#' 
#' @name SpearmanSolver-class
#' 

.SpearmanSolver <- setClass ("SpearmanSolver", contains = "Solver")
#----------------------------------------------------------------------------------------------------
#' Create a Solver class object using Spearman correlation coefficients as the solver
#'
#' @param mtx.assay An assay matrix of gene expression data
#' @param quiet A logical denoting whether or not the solver should print output
#' 
#' @return A Solver class object with Spearman correlation coefficients as the solver
#'
#' @export
#' 
#' @examples
#' solver <- SpearmanSolver()

SpearmanSolver <- function(mtx.assay = matrix(), quiet=TRUE)
{
    obj <- .SpearmanSolver(Solver(mtx.assay=mtx.assay, quiet=quiet))

    # Send a warning if there's a row of zeros
    if(!is.na(max(mtx.assay)) & any(rowSums(mtx.assay) == 0))
       warning("One or more gene has zero expression; this may yield 'NA' results and warnings when using Spearman correlations")

    obj

} #SpearmanSolver, the constructor
#----------------------------------------------------------------------------------------------------
#' Get Spearman Solver Name
#'
#' @param obj An object of class SpearmanSolver
#' 
#' @return "SpearmanSolver"
#'
#' @export
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
#' @rdname solve.Spearman
#' @aliases run.SpearmanSolver solve.Spearman
#' 
#' @description Given a TReNA object with Spearman as the solver, use the \code{\link{cor}} function with
#' \code{method = "spearman"} to esimate coefficients for each transcription factor as a predictor of the target
#' gene's expression level
#'
#' @usage
#' tbl.out <- solve(obj, target.gene, tfs)
#' 
#' @param obj An object of class TReNA with "spearman" as the solver string
#' @param target.gene A designated target gene that should be part of the mtx.assay data
#' @param tfs The designated set of transcription factors that could be associated with the target gene.
#'
#' @return fit The set of Spearman Correlation Coefficients between each transcription factor and the target gene.
#'
#' @seealso \code{\link{cor}}
#'
#' @family solver methods
#' 
#' @examples
#' # Load included Alzheimer's data, create a TReNA object with Bayes Spike as solver, and solve
#' load(system.file(package="TReNA", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
#' trena <- TReNA(mtx.assay = mtx.sub, solver = "pearson")
#' target.gene <- "MEF2C"
#' tfs <- setdiff(rownames(mtx.sub), target.gene)
#' tbl <- solve(trena, target.gene, tfs)

setMethod("run", "SpearmanSolver",

          function (obj, target.gene, tfs){

              # Check if target.gene is in the bottom 10% in mean expression; if so, send a warning              
              if(rowMeans(obj@mtx.assay)[target.gene] < quantile(rowMeans(obj@mtx.assay), probs = 0.1)){                  
                  warning("Target gene mean expression is in the bottom 10% of all genes in the assay matrix")                  
              }              
              
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
              }
              # If target gene was the only tf, then return nothing
              if(length(tfs)==0) return(NULL)

              x = t(mtx[tfs,,drop=F])
              y = as.vector(t(mtx[target.gene,])) # Make target gene levels into a vector

              # Calculate Spearman correlation coefficients
              fit <- cor( x = x, y = y, method = "spearman")

              # Return the coefficients as a data frame in order of coefficient size

              tbl <- data.frame(row.names = rownames(fit)[order(abs(fit), decreasing = TRUE)],
                                coefficient = fit[order(abs(fit), decreasing = TRUE)])

              return(tbl)
          })
#----------------------------------------------------------------------------------------------------
