#' @title Create a NullFilter object
#'
#' @description
#' A NullFilter object allows for "filtering" of the genes in an assay matrix. Its associated
#' \code{getCandidates} method returns the list of transcription factors included in the assay matrix.
#' 
#' @include CandidateFilter.R
#' @import methods
#' 
#' @name NullFilter-class
#' @rdname NullFilter-class
#' @aliases NullFilter

#----------------------------------------------------------------------------------------------------
.NullFilter <- setClass ("NullFilter", contains = "CandidateFilter")

#----------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))

#----------------------------------------------------------------------------------------------------
#' @rdname NullFilter-class
#' 
#' @param mtx.assay An assay matrix of gene expression data
#' @param quiet A logical denoting whether or not the filter should print output
#' @return A CandidateFilter class object with null as the filtering method
#'
#' @export
#' 
#' @seealso \code{\link{getCandidates-NullFilter}}
#'
#' @family Filtering Objects
#'
#' @examples
#' load(system.file(package="TReNA", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
#' null.filter <- NullFilter(mtx.assay = mtx.sub)

NullFilter <- function(mtx.assay=matrix(), quiet = TRUE)
{
    .NullFilter(CandidateFilter(mtx.assay = mtx.assay, quiet = quiet))

} # NullFilter, the constructor
#----------------------------------------------------------------------------------------------------
#' Get candidate genes using the null filter
#'
#' @aliases getCandidates-NullFilter
#'
#' @usage
#' tfs <- getCandidates(null.filter)
#' 
#' @param obj An object of class NullFilter
#' @return A vector containing all genes in the assay matrix
#'
#' @seealso \code{\link{NullFilter}}
#' 
#' @family getCandidate Methods
#'
#' @examples
#'
#' # Using the included Alzheimer's data, return all transcription factors as candidates
#' load(system.file(package="TReNA", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
#' null.filter <- NullFilter(mtx.assay=mtx.sub)
#' tfs <- getCandidates(null.filter)

setMethod("getCandidates", "NullFilter",

    function(obj){
        # Simply return the genes
	genes <- rownames(obj@mtx.assay)
	return(genes)
	})
	
#----------------------------------------------------------------------------------------------------
