#' @import methods
#'
#' @name CandidateFilter-class
#' @rdname CandidateFilter-class
#' @aliases CandidateFilter
#'
#' @slot mtx.assay An assay matrix of gene expression data
#' @slot quiet A logical denoting whether or not the CandidateFilter object should print output
#'
#' @family Filtering objects

.CandidateFilter <- setClass("CandidateFilter",
                    slots = c(quiet = "logical"
			      )
                    )
#----------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#----------------------------------------------------------------------------------------------------
#' Get candidate genes using a CandidateFilter object
#'
#' @rdname getCandidates
#' @aliases getCandidates
#'
#' @param obj An object of a CandidateFilter class
#' @param extraArgs A named list of extra arguments corresponding to the chosen filter
#'
#' @return A vector containing genes from the assay matrix that are selected by the filter
#'
#' @family getCandidate Methods
#' @export
setGeneric("getCandidates", signature="obj", function(obj, argsList) standardGeneric("getCandidates"))

#' Retrieve the assay matrix of gene expression data
#'
#' @aliases getFilterAssayData
#'
#' @param obj An object of a CandidateFilter class
#'
#' @return The assay matrix of gene expression data associated with a CandidateFilter object
#'
#' @examples
#'
#' # Create a CandidateFilter object using the included Alzheimer's data and retrieve the matrix
#' load(system.file(package="TReNA", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
#' my.filter <- CandidateFilter(mtx.sub)
#' mtx <- getFilterAssayData(my.filter)

#' @export
setGeneric("getFilterAssayData",    signature="obj", function(obj) standardGeneric ("getFilterAssayData"))
#----------------------------------------------------------------------------------------------------
#' CandidateFilter
#'
#' A CandidateFilter is an S4 class to represent a gene candidate filter. These filters can employ a variety of methods
#' to reduce the number of transcription factors used as predictors for solving a TReNA object.
#'
#' @rdname CandidateFilter-class
#' @aliases CandidateFilter
#'
#' @param mtx.assay An assay matrix of gene expression data
#' @param quiet A logical denoting whether or not the CandidateFilter object should print output
#'
#' @return An object of the Candidate filter class
#'
#' @seealso \code{\link{getCandidates}}, \code{\link{getFilterAssayData}}
#'
#' @export
#'
#' @examples
#' # Create an empty candidate filter
#' candidate.filter <- CandidateFilter(mtx.assay = matrix(), quiet=TRUE)

CandidateFilter <- function(quiet = TRUE)
{
    .CandidateFilter(quiet = quiet)

} # CandidateFilter, the constructor
#----------------------------------------------------------------------------------------------------
#' @describeIn CandidateFilter Retrieve the assay matrix of gene expression data
#'
#' @param obj An object of a CandidateFilter class
#' 
#' @examples
#'
#' # Create a CandidateFilter object using the included Alzheimer's data and retrieve the matrix
#' load(system.file(package="TReNA", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
#' my.filter <- CandidateFilter(mtx.sub)
#' mtx <- getFilterAssayData(my.filter)

setMethod("getFilterAssayData", "CandidateFilter",

   function (obj){
      obj@mtx.assay
      })
#----------------------------------------------------------------------------------------------------
