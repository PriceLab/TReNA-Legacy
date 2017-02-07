#' Apply a null filter
#'
#' @param mtx.assay An assay matrix of gene expression data
#' @return An object of class NullFilter

#----------------------------------------------------------------------------------------------------
.NullFilter <- setClass ("NullFilter",
                     slots = c(mtx.assay = "matrix",
		               quiet = "logical"
                     )
)

#----------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#----------------------------------------------------------------------------------------------------
setGeneric("getCandidates", signature="obj", function(obj) standardGeneric("getCandidates"))

#----------------------------------------------------------------------------------------------------
NullFilter <- function(mtx.assay=matrix(), quiet = TRUE)
{
    # Simply return the genes
    .NullFilter(mtx.assay = mtx.assay, quiet = quiet)

} # NullFilter, the constructor
#----------------------------------------------------------------------------------------------------
setMethod("getCandidates", "NullFilter",

    function(obj){
        # Simply return the genes
	genes <- rownames(obj@mtx.assay)
	return(genes)
	})
	
#----------------------------------------------------------------------------------------------------
