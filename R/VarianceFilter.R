#' Filter based on gene expression variance
#'
#' @include CandidateFilter.R
#' @name VarianceFilter-class
#' @param mtx.assay An assay matrix of gene expression data

#----------------------------------------------------------------------------------------------------
.VarianceFilter <- setClass("VarianceFilter", contains = "CandidateFilter")

#----------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#----------------------------------------------------------------------------------------------------
VarianceFilter <- function(mtx.assay=matrix(), quiet=TRUE)
{
    .VarianceFilter(CandidateFilter(mtx.assay = mtx.assay, quiet = quiet))

} # VarianceFilter, the constructor
#----------------------------------------------------------------------------------------------------
#' Get candidate genes using the variance filter
#'
#' @aliases getCandidates-VarianceFilter
#'
#' @param obj An object of class VarianceFilter
#' @param target.gene A designated target gene that should be part of the mtx.assay data
#' @param var.size A user-specified percentage (0-1) of the target gene variance to use as a filter
#' (default = 0.5)
#'
#' @return A vector containing all genes with variances less than the target gene

setMethod("getCandidates", "VarianceFilter",

    function(obj,target.gene, var.size = 0.5){
        # Designate the target genes and tfs
	tfs <- setdiff(rownames(obj@mtx.assay), target.gene)
	tf.mtx <- obj@mtx.assay[-c(which(rownames(obj@mtx.assay) == target.gene)),]
	target.mtx <- obj@mtx.assay[which(rownames(obj@mtx.assay) == target.gene),]

	# Find the variances
	tf.var <- apply(tf.mtx,1,var)
	target.var <- var(target.mtx)

	# Return only the genes with variances within 50% of target gene variance
	return(names(which(tf.var > (1-var.size)*target.var & tf.var < (1+var.size)*target.var)))
	}
)
#----------------------------------------------------------------------------------------------------
