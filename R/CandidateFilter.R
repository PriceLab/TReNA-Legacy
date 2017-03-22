#----------------------------------------------------------------------------------------------------
#' An S4 class to represent a gene candidate filter
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
setGeneric("getCandidates", signature="obj", function(obj,...) standardGeneric("getCandidates"))

#----------------------------------------------------------------------------------------------------
#' Get candidate genes using the selected filter
#'
#' @name getCandidates-method
#' @rdname getCandidates
#'
#' @param mtx.assay An assay matrix of gene expression data
#'
#' @return A vector containing all genes selected as candidates by the filter

CandidateFilter <- function(mtx.assay = matrix(), quiet = TRUE)
{
    .CandidateFilter(mtx.assay = mtx.assay, quiet = quiet)

} # CandidateFilter, the constructor
#----------------------------------------------------------------------------------------------------
