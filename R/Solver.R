#----------------------------------------------------------------------------------------------------
#' @name Solver-class
#' @rdname Solver-class
#' @aliases Solver
#' 
#' @import methods

.Solver <- setClass ("Solver",
                     slots = c(mtx.assay="matrix",
                               targetGene="character",
                               candidateRegulators="character",
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

#' Retrieve the target gene from a Solver object
#'
#' @rdname getTarget
#' @aliases getTarget
#'
#' @param obj An object of class Solver
#'
#' @return The target gene associated with a Solver object
#'
#' @examples
#' # Create a Solver object using the included Alzheimer's data and retrieve the target gene
#' load(system.file(package="TReNA", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
#' solver <- Solver(mtx.sub)
#' mtx <- getTarget(solver) 

#' @export
setGeneric("getTarget", signature = "obj", function(obj) standardGeneric("getTarget"))

#' Retrieve the target gene from a Solver object
#'
#' @rdname getTarget
#' @aliases getTarget
#'
#' @param obj An object of class Solver
#'
#' @return The target gene associated with a Solver object
#'
#' @examples
#' # Create a Solver object using the included Alzheimer's data and retrieve the regulators
#' load(system.file(package="TReNA", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
#' solver <- Solver(mtx.sub)
#' mtx <- getRegulators(solver) 

#' @export
setGeneric("getRegulators", signature = "obj", function(obj) standardGeneric("getRegulators"))

#' Show a Solver object
#'
#' @rdname show
#' @aliases show
#'
#' @param obj An object of class Solver
#'
#' @return The names of the top 10 candidate regulators, plus statistics on the assay matrix, the
#' target gene, and the total number of candidate regulators
#'
#' @examples
#' load(system.file(package="TReNA", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
#' tfs <- setdiff(rownames(mtx.sub),"MEF2C")
#' pearson.solver <- PearsonSolver(mtx.sub, "MEF2C", tfs)
#' show(pearson.solver)

setGeneric("show", signature = "obj", function(obj) standardGeneric("show"))
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

Solver <- function(mtx.assay=matrix(), targetGene, candidateRegulators, quiet=TRUE)
{
    # If a matrix is supplied, check the distribution to see if it's too big
    if(!is.na(max(mtx.assay))){
        mtx.ratio <- (max(mtx.assay) - stats::quantile(mtx.assay,0.75))/(stats::quantile(mtx.assay,0.75) - stats::median(mtx.assay))
        if(mtx.ratio > 1000){
            warning("Assay matrix may contain highly skewed data; consider transforming your matrix.")
            }
    }

    

    env <- new.env(parent=emptyenv())
    .Solver(mtx.assay=mtx.assay,
            targetGene = targetGene,
            candidateRegulators = candidateRegulators,
            quiet=quiet,
            state=env)

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
#' @describeIn Solver Retrieve the target gene for a Solver
#'
#' @param obj An object of class Solver
#' 
#' @examples
#'
#' # Create a Solver object using the included Alzheimer's data and retrieve the matrix
#' load(system.file(package="TReNA", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
#' solver <- Solver(mtx.sub)
#' mtx <- getTarget(solver)

setMethod("getTarget", "Solver",

   function (obj){
      obj@targetGene
      })
#----------------------------------------------------------------------------------------------------
#' @describeIn Solver Retrieve the candidate regulators for a Solver
#'
#' @param obj An object of class Solver
#' 
#' @examples
#'
#' # Create a Solver object using the included Alzheimer's data and retrieve the matrix
#' load(system.file(package="TReNA", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
#' solver <- Solver(mtx.sub)
#' mtx <- getRegulators(solver)

setMethod("getRegulators", "Solver",

   function (obj){
      obj@candidateRegulators
   })
#----------------------------------------------------------------------------------------------------
setMethod('show', 'Solver',

    function(obj) {
       regulator.count <- length(getRegulators(obj))
       if(regulator.count > 10){
          regulatorString <- paste(getRegulators(obj)[1:10], collapse=",")
          regulatorString <- sprintf("%s...", regulatorString);
          }
       else
          regulatorString <- paste(getRegulators(obj), collapse=",")

       msg = sprintf("RandomForestSolver with mtx.assay (%d, %d), targetGene %s, %d candidate regulators %s",
                     nrow(getAssayData(obj)), ncol(getAssayData(obj)),
                     getTarget(obj), regulator.count, regulatorString)
       cat (msg, '\n', sep='')
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
