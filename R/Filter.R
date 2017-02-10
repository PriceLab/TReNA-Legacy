#----------------------------------------------------------------------------------------------------
#' An S4 class to represent a filter
#'
#' @name Filter-class
#' @param mtx.assay An assay matrix of gene expression data

.Filter <- setClass("Filter",
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
#' @name getCandidates
#' @return A vector containing all genes selected as candidates by the filter

Filter <- function(mtx.assay = matrix(), quiet = TRUE)
{
    .Filter(mtx.assay = mtx.assay, quiet = quiet)

} # Filter, the constructor
#----------------------------------------------------------------------------------------------------
