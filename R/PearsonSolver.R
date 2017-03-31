#----------------------------------------------------------------------------------------------------
#' An S4 class to represent a Pearson solver
#'
#' @include Solver.R
#' @import methods
#' 
#' @name PearsonSolver-class
#' 

.PearsonSolver <- setClass ("PearsonSolver", contains = "Solver")
#----------------------------------------------------------------------------------------------------
#' Create a Solver class object using  Pearson correlation coefficients as the solver
#'
#' @param mtx.assay An assay matrix of gene expression data
#' @param quiet A logical denoting whether or not the solver should print output
#' 
#' @return A Solver class object with Pearson correlation coefficients as the solver
#'
#' @export
#' 
#' @examples
#' solver <- PearsonSolver()

PearsonSolver <- function(mtx.assay = matrix(), quiet=TRUE)
{
    obj <- .PearsonSolver(Solver(mtx.assay=mtx.assay, quiet=quiet))


    # Send a warning if there's a row of zeros
    if(!is.na(max(mtx.assay)) & any(rowSums(mtx.assay) == 0))
       warning("One or more gene has zero expression; this may yield 'NA' results and warnings when using Pearson correlations")

    obj

} #PearsonSolver, the constructor
#----------------------------------------------------------------------------------------------------
#' Run the Pearson Solver
#'
#' @rdname solve.Pearson
#' @aliases run.PearsonSolver solve.Pearson
#'
#' @description Given a TReNA object with Pearson as the solver, use the \code{\link{cor}} function to
#' estimate coefficients for each transcription factor as a perdictor of the target gene's expression level
#'
#' @param obj An object of class TReNA with "pearson" as the solver string
#' @param target.gene A designated target gene that should be part of the mtx.assay data
#' @param tfs The designated set of transcription factors that could be associated with the target gene.
#' @param tf.weights A set of weights on the transcription factors (default = rep(1, length(tfs)))
#' @param extraArgs Modifiers to the Pearson solver
#'
#' @return The set of Pearson Correlation Coefficients between each transcription factor and the target gene.
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

setMethod("run", "PearsonSolver",

          function (obj, target.gene, tfs, tf.weights=rep(1,length(tfs)), extraArgs=list()){

              # Check if target.gene is in the bottom 10% in mean expression; if so, send a warning              
              if(rowMeans(obj@mtx.assay)[target.gene] < stats::quantile(rowMeans(obj@mtx.assay), probs = 0.1)){                  
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

              # Calculate Pearson correlation coefficients
              fit <- stats::cor( x = x, y = y)

              # Return the coefficients as a data frame 
              tbl <- data.frame(row.names = rownames(fit)[order(abs(fit), decreasing = TRUE)],
                                coefficient = fit[order(abs(fit), decreasing = TRUE)])

              return(tbl)
          })
#----------------------------------------------------------------------------------------------------
