#' Filter based on gene expression variance
#'
#' @include Filter.R
#' @name VarianceFilter-class
#' @param mtx.assay An assay matrix of gene expression data

#----------------------------------------------------------------------------------------------------
.VarianceFilter <- setClass("VarianceFilter", contains = "Filter")

#----------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#----------------------------------------------------------------------------------------------------
VarianceFilter <- function(mtx.assay=matrix(), quiet=TRUE)
{
    .VarianceFilter(Filter(mtx.assay = mtx.assay, quiet = quiet))

} # VarianceFilter, the constructor
#----------------------------------------------------------------------------------------------------
#' Get candidate genes using the variance filter
#'
#' @aliases getCandidates-VarianceFilter
#' @param target.gene A designated target gene that should be part of the mtx.assay data
#'
#' @return A vector containing all genes with variances less than the target gene

setMethod("getCandidates", "VarianceFilter",

    function(obj,target.gene){
        # Designate the target genes and tfs
	tfs <- setdiff(rownames(obj@mtx.assay), target.gene)
	tf.mtx <- mtx.assay[-c(which(rownames(obj@mtx.assay) == target.gene)),]
	target.mtx <- mtx.assay[which(rownames(obj@mtx.assay) == target.gene),]

	# Find the variances
	tf.var <- apply(tf.mtx,1,var)
	target.var <- var(target.mtx)

	# Return only the genes with variances less than the target gene
	return(names(tf.var < target.var))
	}
)
#----------------------------------------------------------------------------------------------------
