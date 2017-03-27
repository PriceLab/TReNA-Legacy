#------------------------------------------------------------------------------------------------------------------------
#' An S4 class to represent a LASSO P-Value solver
#'
#' @include Solver.R
#' @import lassopv
#' @import methods
#' 
#' @name LassoPVSolver-class

.LassoPVSolver <- setClass ("LassoPVSolver", contains="Solver")
#------------------------------------------------------------------------------------------------------------------------
#' Create a Solver class object using the LASSO P-Value solver
#'
#' @param mtx.assay An assay matrix of gene expression data
#' @param quiet A logical denoting whether or not the solver should print output
#' 
#' @return A Solver class object with LASSO P-Value as the solver
#'
#' @export
#' 
#' @examples
#' solver <- LassoPVSolver()

LassoPVSolver <- function(mtx.assay=matrix(), quiet=TRUE)
{
    obj <- .LassoPVSolver(Solver(mtx.assay=mtx.assay, quiet=quiet))

    # Send a warning if there's a row of zeros
    if(!is.na(max(mtx.assay)) & any(rowSums(mtx.assay) == 0))
       warning("One or more gene has zero expression; this may cause problems when using P-Value LASSO. You may want to try 'lasso' or 'ridge' instead.")

    obj

} # LassoPVSolver, the constructor
#------------------------------------------------------------------------------------------------------------------------
#' Get LassoPV Solver name
#'
#' @param obj An object of class LassoPVSolver
#' 
#' @return "LassoPVSolver"
#'
#' @export
#' 
#' @examples
#' solver <- LassoPVSolver()
#' getSolverName(solver)

setMethod("getSolverName", "LassoPVSolver",

  function (obj){
     return("LassoPVSolver")
     })

#----------------------------------------------------------------------------------------------------
#' Run the LASSO P-Value Solver
#'
#' @rdname solve.LassoPV
#' @aliases run.LassoPVSolver solve.LassoPV
#' 
#' @description Given a TReNA object with LASSO P-Value as the solver, use the \code{\link{lassopv}} function to estimate coefficients
#' for each transcription factor as a predictor of the target gene's expression level. 
#'
#' @usage
#' solve(obj, target.gene, tfs, tf.weights=rep(1,length(tfs)), extraArgs=list())
#' 
#' @param obj An object of class TReNA with "lassopv" as the solver string
#' @param target.gene A designated target gene that should be part of the mtx.assay data
#' @param tfs The designated set of transcription factors that could be associated with the target gene.
#' @param tf.weights A set of weights on the transcription factors (default = rep(1, length(tfs)))
#' @param extraArgs Modifiers to the Lasso P-Value solver
#'
#' @return A data frame containing the p-values for each transcription factor pertaining to the target gene
#' plus the Pearson correlations between each transcription factor and the target gene.
#'
#' @seealso \code{\link{lassopv}}
#'
#' @family solver methods
#' 
#' @examples
#' # Load included Alzheimer's data, create a TReNA object with Bayes Spike as solver, and solve
#' load(system.file(package="TReNA", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
#' trena <- TReNA(mtx.assay = mtx.sub, solver = "lassopv")
#' target.gene <- "MEF2C"
#' tfs <- setdiff(rownames(mtx.sub), target.gene)
#' tbl <- solve(trena, target.gene, tfs)

setMethod("run", "LassoPVSolver",

          function (obj, target.gene, tfs, tf.weights=rep(1,length(tfs)), extraArgs=list()){
              
              if(length(tfs) == 0)                  
                  return(data.frame())              

              # Check if target.gene is in the bottom 10% in mean expression; if so, send a warning         
              if(rowMeans(obj@mtx.assay)[target.gene] < stats::quantile(rowMeans(obj@mtx.assay), probs = 0.1)){                 
                  warning("Target gene mean expression is in the bottom 10% of all genes in the assay matrix")            
              }
              
              # we don't try to handle tf self-regulation           
              deleters <- grep(target.gene, tfs)
              if(length(deleters) > 0){                  
                  tfs <- tfs[-deleters]                  
                  if(!obj@quiet)                   
                      message(sprintf("LassoPVSolver removing target.gene from candidate regulators: %s", target.gene))                  
              }
              
              if( length(tfs) == 0 ) return( data.frame() )
              
              mtx <- obj@mtx.assay              
              stopifnot(target.gene %in% rownames(mtx))             
              stopifnot(all(tfs %in% rownames(mtx)))              
              features <- t(mtx[tfs,,drop=F ])              
              target <- as.numeric(mtx[target.gene,])

              # Run LASSO P-Value and return the P-values, ordered by increasing value
              fit <- lassopv(features, target)
              fit <- fit[order(fit, decreasing=FALSE)]

              # Add pearson correlations and make a data frame              
              correlations.of.betas.to.targetGene <- unlist(lapply(names(fit),
                                                                   function(x) stats::cor(mtx[x,], mtx[target.gene,])))
              tbl <- data.frame(row.names = names(fit),
                                p.values = fit,
                                gene.cor=correlations.of.betas.to.targetGene)
              return(tbl)
})
#----------------------------------------------------------------------------------------------------
#' Rescale LASSO P-Value Predictor Weights
#'
#' @aliases rescalePredictorWeights.LassoPVSolver
#'
#' @param obj An object of class LassoPVSolver
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

setMethod("rescalePredictorWeights", "LassoPVSolver",

   function (obj, rawValue.min, rawValue.max, rawValues){
      1 - ((rawValues-rawValue.min)/(rawValue.max-rawValue.min))
      })

#----------------------------------------------------------------------------------------------------

