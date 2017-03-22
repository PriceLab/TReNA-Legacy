#----------------------------------------------------------------------------------------------------
#' @title Solver Class
#'
#' @name Solver-class
#' @rdname Solver-class
#'
#' @export
#'
#' @description
#' The Solver class is a generic class that governs the different solvers available in TReNA. A
#' Solver class object is constructed during creation of a TReNA object and resides within the
#' TReNA object. It is rarely called by itself; rather, interaction with a particular solver object
#' is achieved using the \code{\link{solve}} method on a TReNA object. 
#' 
#' @slot mtx.assay An assay matrix of gene expression data
#' @slot quiet A logical element indicating whether the Solver object should print output

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
#' @rdname Solver-class
#' 
#' @param mtx.assay An assay matrix of gene expression data
#' @param quiet A logical indicating whether or not the Solver object should print output
#'
#' @return An object of the Solver class
#'
#' @examples
#' mtx <- matrix(rnorm(10000), nrow = 100)
#' solver <- Solver(mtx)
#'
#' @seealso \code{\link{TReNA}}, \code{\link{solve}}

Solver <- function(mtx.assay=matrix(), quiet=TRUE)
{

    # If a matrix is supplied, check the distribution to see if it's too big
    if(!is.na(max(mtx.assay))){        
        mtx.ratio <- (max(mtx.assay) - quantile(mtx.assay,0.75))/(quantile(mtx.assay,0.75) - median(mtx.assay))        
        if(mtx.ratio > 1000){                    
            warning("Assay matrix may contain highly skewed data; consider transforming your matrix.")
            }
    }
    
    env <- new.env(parent=emptyenv())
   .Solver(mtx.assay=mtx.assay, quiet=quiet, state=env)

} # TReNA, the constructor
#----------------------------------------------------------------------------------------------------
#' @title Get Assay Data from Solver
#'
#' @description
#' Retrieve the assay matrix of gene expression data from a Solver object
#' 
#' @name getAssayData-methods
#' @aliases getAssayData
#' 
#' @param obj An object of class Solver
#' 
#' @return The assay matrix of gene expression data associated with a Solver object
#'
#' @examples
#'
#' # Create a Solver object using the included Alzheimer's data and retrieve the matrix
#' load(system.file(package="TReNA", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
#' solver <- Solver(mtx.sub)
#' mtx <- getAssayData(solver)

setMethod("getAssayData", "Solver",

   function (obj){
      obj@mtx.assay
      })

#----------------------------------------------------------------------------------------------------
