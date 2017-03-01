#----------------------------------------------------------------------------------------------------
#' An S4 class to represent a solver
#'
#' @slot mtx.assay An assay matrix of gene expression data
#' @slot quiet A logical element indicating whether the solver should produce output
#' @slot state Environment variable

.Solver <- setClass ("Solver",
                     slots = c(mtx.assay="matrix",
                               quiet="logical",
                               state="environment")
                     )

#----------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#----------------------------------------------------------------------------------------------------
setGeneric("getSolverName",   signature="obj", function(obj, target.gene, tfs) standardGeneric ("getSolverName"))
setGeneric("getAssayData",    signature="obj", function(obj) standardGeneric ("getAssayData"))
setGeneric("run",             signature="obj", function(obj, target.gene, tfs, tf.weights, extraArgs=list()) standardGeneric ("run"))

setGeneric("rescalePredictorWeights",
                              signature="obj", function(obj, rawValue.min, rawValue.max, rawValues) standardGeneric ("rescalePredictorWeights"))
#----------------------------------------------------------------------------------------------------
#' Define an object of class Solver
#'
#' @param mtx.assay An assay matrix of gene expression data
#'
#' @return An object of the Solver class

Solver <- function(mtx.assay=matrix(), quiet=TRUE)
{

    # If a matrix is supplied, check the distribution to see if it's too big
    if(!is.na(max(mtx.assay)) & (max(mtx.assay) - min(mtx.assay)) > 1E4)
        warning("Matrix range exceeds 10,000 and may indicate skewed data; consider transforming your matrix.")

    env <- new.env(parent=emptyenv())
   .Solver(mtx.assay=mtx.assay, quiet=quiet, state=env)

} # TReNA, the constructor
#----------------------------------------------------------------------------------------------------
#' Get Assay Data from Solver
#'
#' @param obj An object of class Solver
#' 
#' @return The assay matrix of gene expression data associated with a Solver object

setMethod("getAssayData", "Solver",

   function (obj){
      obj@mtx.assay
      })

#----------------------------------------------------------------------------------------------------
