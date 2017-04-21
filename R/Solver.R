#----------------------------------------------------------------------------------------------------
#' @name Solver-class
#' @rdname Solver-class
#' @aliases Solver
#' 
#' @import methods

.Solver <- setClass ("Solver",
                     slots = c(mtx.assay="matrix",
                               quiet="logical",
                               state="environment")
                     )

#----------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#----------------------------------------------------------------------------------------------------
#' Retrieve the assay matrix of gene expression data from a Solver object
#' 
#' @rdname getAssayData
#' @aliases getAssayData
#' 
#' @param obj An object of class Solver
#'
#' @export
#' 
#' @return The assay matrix of gene expression data associated with a Solver object
#'
#' @examples
#' # Create a Solver object using the included Alzheimer's data and retrieve the matrix
#' load(system.file(package="TReNA", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
#' solver <- Solver(mtx.sub)
#' mtx <- getAssayData(solver)
#' 
setGeneric("getAssayData",    signature="obj", function(obj) standardGeneric ("getAssayData"))
#' @export
setGeneric("run",             signature="obj", function(obj, target.gene, tfs, tf.weights, extraArgs=list()) standardGeneric ("run"))
#' @export
setGeneric("rescalePredictorWeights",
                              signature="obj", function(obj, rawValue.min, rawValue.max, rawValues) standardGeneric ("rescalePredictorWeights"))
#----------------------------------------------------------------------------------------------------
#' Define an object of class Solver
#'
#' @description
#' The Solver class is a generic class that governs the different solvers available in TReNA. A
#' Solver class object is constructed during creation of a TReNA object and resides within the
#' TReNA object. It is rarely called by itself; rather, interaction with a particular solver object
#' is achieved using the \code{\link{solve}} method on a TReNA object.
#' 
#' @rdname Solver-class
#' 
#' @param mtx.assay An assay matrix of gene expression data
#' @param quiet A logical indicating whether or not the Solver object should print output
#'
#' @export
#' 
#' @return An object of the Solver class
#'
#' @examples
#' # Create a simple Solver object with default options
#' mtx <- matrix(rnorm(10000), nrow = 100)
#' solver <- Solver(mtx)
#'
#' @seealso \code{\link{getAssayData}}, \code{\link{TReNA}}, \code{\link{solve}}
#'
#' @family Solver class objects

Solver <- function(mtx.assay=matrix(), quiet=TRUE)
{
    # If a matrix is supplied, check the distribution to see if it's too big
    if(!is.na(max(mtx.assay))){        
        mtx.ratio <- (max(mtx.assay) - stats::quantile(mtx.assay,0.75))/(stats::quantile(mtx.assay,0.75) - stats::median(mtx.assay))        
        if(mtx.ratio > 1000){                    
            warning("Assay matrix may contain highly skewed data; consider transforming your matrix.")
            }
    }
    
    env <- new.env(parent=emptyenv())
   .Solver(mtx.assay=mtx.assay, quiet=quiet, state=env)

} # Solver, the constructor
#----------------------------------------------------------------------------------------------------
#' @describeIn Solver Retrieve the assay matrix of gene expression data
#'
#' @param obj An object of class Solver
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
#' Rescale the Predictor Weights
#'
#' Solvers such as LASSO penalize predictors on a scale of 1 (full weight) to infinity (zero weight).
#' With the \code{rescalePredictorWeights} method, incoming raw values can be scaled between a possibly
#' theoretical minimum and maximum value. 
#'
#' @rdname rescalePredictorWeights
#' @aliases rescalePredictorWeights
#'
#' @param obj An object of the Solver class
#' @param rawValue.min The minimum value of the raw expression values
#' @param rawValue.max The maximum value of the raw expression values
#' @param rawValues A matrix of raw expression values
#'
#' @export
#'
#' @return A matrix of the raw values re-scaled using the minimum and maximum values
#'
#' @examples
#' # Create a LassoSolver object using the included Alzheimer's data and rescale the predictors
#' load(system.file(package="TReNA", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
#' ls <- LassoSolver(mtx.sub)
#' raw.values <- c(241, 4739, 9854, 22215, 658334)
#' cooked.values <- rescalePredictorWeights(ls, rawValue.min = 1, rawValue.max = 1000000, raw.values)

setMethod("rescalePredictorWeights", "Solver",

          function(obj, rawValue.min, rawValue.max, rawValues){
              1 - ((rawValues-rawValue.min)/(rawValue.max-rawValue.min))           
          })
#----------------------------------------------------------------------------------------------------
