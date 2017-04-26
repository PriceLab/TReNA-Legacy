#' Class CandidateFilter
#'
#' @import methods
#'
#' @name CandidateFilter-class
#' @rdname CandidateFilter-class
#' @aliases CandidateFilter
#'
#' @slot mtx.assay An assay matrix of gene expression data
#' @slot quiet A logical denoting whether or not the CandidateFilter object should print output
#'
#' @return An object of the CandidateFilter class
#'
#' @seealso \code{\link{FootprintFilter}} \code{\link{VarianceFilter}} \code{\link{NullFilter}}

.CandidateFilter <- setClass("CandidateFilter",
                    slots = c(quiet = "logical"
			      )
                    )
#----------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#----------------------------------------------------------------------------------------------------
#' @export
setGeneric("getCandidates", signature="obj", function(obj, argsList) standardGeneric("getCandidates"))

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
#'
NULL
#----------------------------------------------------------------------------------------------------
