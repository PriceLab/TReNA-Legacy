#----------------------------------------------------------------------------------------------------
#' An S4 class to represent a gene candidate filter
#'
#' @import methods
#' 
#' @name CandidateFilter-class
#' @param mtx.assay An assay matrix of gene expression data

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
#' @rdname CandidateFilter-Class
#'
#' @param mtx.assay An assay matrix of gene expression data
#'
#' @return An object of the Candidate filter class

CandidateFilter <- function(mtx.assay = matrix(), quiet = TRUE)
{
    .CandidateFilter(mtx.assay = mtx.assay, quiet = quiet)

} # CandidateFilter, the constructor
#----------------------------------------------------------------------------------------------------
#' Get candidate genes using the selected
