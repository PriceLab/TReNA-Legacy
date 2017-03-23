#' @title Class CandidateFilter
#'
#' @description
#' An S4 class to represent a gene candidate filter. These filters can employ a variety of methods
#' to reduce the number of transcription factors used as predictors for solving a TReNA object.
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
                    slots = c(mtx.assay = "matrix",
		              quiet = "logical"
			      )
)
#----------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#----------------------------------------------------------------------------------------------------
#' @export
setGeneric("getCandidates", signature="obj", function(obj,...) standardGeneric("getCandidates"))

#----------------------------------------------------------------------------------------------------
#' @title Class CandidateFilter
#' @name CandidateFilter-Class
#' @rdname CandidateFilter-Class
#'
#' @param mtx.assay An assay matrix of gene expression data
#' @param quiet A logical denoting whether or not the CandidateFilter object should print output
#'
#' @return An object of the Candidate filter class
#'
#' @examples
#' # Create an empty candidate filter
#' candidate.filter <- CandidateFilter(mtx.assay = matrix(), quiet=TRUE)

CandidateFilter <- function(mtx.assay = matrix(), quiet = TRUE)
{
    .CandidateFilter(mtx.assay = mtx.assay, quiet = quiet)

} # CandidateFilter, the constructor
#----------------------------------------------------------------------------------------------------
