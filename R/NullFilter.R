#' Apply a null filter
#'
#' @include CandidateFilter.R
#' @name NullFilter-class
#'
#' @param mtx.assay An assay matrix of gene expression data
#' @return An object of class NullFilter

#----------------------------------------------------------------------------------------------------
.NullFilter <- setClass ("NullFilter", contains = "CandidateFilter")

#----------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))

#----------------------------------------------------------------------------------------------------
#' Define an object of class Null Filter
#'
#' @param mtx.assay An assay matrix of gene expression data
#'
#' @return An object of the Null Filter class

NullFilter <- function(mtx.assay=matrix(), quiet = TRUE)
{
    .NullFilter(CandidateFilter(mtx.assay = mtx.assay, quiet = quiet))

} # NullFilter, the constructor
#----------------------------------------------------------------------------------------------------
#' Get candidate genes using the null filter
#'
#' @aliases getCandidates-NullFilter
#'
#' @param obj An object of class NullFilter
#' @return A vector containing all genes in the assay matrix

setMethod("getCandidates", "NullFilter",

    function(obj){
        # Simply return the genes
	genes <- rownames(obj@mtx.assay)
	return(genes)
	})
	
#----------------------------------------------------------------------------------------------------
