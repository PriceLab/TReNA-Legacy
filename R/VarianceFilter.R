#' @title Create a VarianceFilter object
#'
#' @description
#' A VarianceFilter object allows for filtering based on the variance of a target gene in relation to
#' other genes in the assay matrix. Using its associated \code{getCandidates} method, a VarianceFilter
#' object can be used to filter a list of possible transcription factors to those within a given
#' range of the variance of a supplied target gene.
#'
#' @include CandidateFilter.R
#' @import methods
#'
#' @name VarianceFilter-class
#' @rdname VarianceFilter-class
#' @aliases VarianceFilter

#----------------------------------------------------------------------------------------------------
.VarianceFilter <- setClass("VarianceFilter", contains = "CandidateFilter")

#----------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#----------------------------------------------------------------------------------------------------
#' @rdname VarianceFilter-class
#'
#' @param mtx.assay An assay matrix of gene expression data
#' @param quiet A logical denoting whether or not the solver should print output
#'
#' @return A CandidateFilter class object with variance as the filtering method
#'
#' @export
#' 
#' @seealso \code{\link{getCandidates-VarianceFilter}}, \code{\link{getFilterAssayData}}
#'
#' @return An object of the VarianceFilter class
#' 
#' @family Filtering Objects
#'
#' @examples
#' load(system.file(package="TReNA", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
#' variance.filter <- VarianceFilter(mtx.assay = mtx.sub)

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
#' @param argsList A named list containing two fields:
#' \itemize{
#' \item{"target.gene" A designated target gene that should be part of the mtx.assay data}
#' \item{"var.size" A user-specified percentage (0-1) of the target gene variance to use as a filter}
#' }
#' @seealso \code{\link{VarianceFilter}}
#'
#' @family getCandidate Methods
#'
#' @return A vector containing all genes with variances less than the target gene
#'
#' @examples
#'
#' # Using the included Alzheimer's dataset, filter out only those transcription factors with variance
#' # within 50% of the variance of MEF2C
#' load(system.file(package="TReNA", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
#' variance.filter <- VarianceFilter(mtx.assay = mtx.sub)
#'
#' target.gene <- "MEF2C"
#' tfs <- getCandidates(variance.filter, argsList = list("target.gene" = target.gene, "var.size" = 0.5))

setMethod("getCandidates", "VarianceFilter",

          function(obj,argsList){

              # Retrive the extra arguments
              target.gene <- argsList[["target.gene"]]
              var.size <- argsList[["var.size"]]

              # Designate the target genes and tfs              
              tfs <- setdiff(rownames(getFilterAssayData(obj)), target.gene)              
              tf.mtx <- getFilterAssayData(obj)[-c(which(rownames(getFilterAssayData(obj)) == target.gene)),]              
              target.mtx <- getFilterAssayData(obj)[which(rownames(getFilterAssayData(obj)) == target.gene),]
              
              # Find the variances              
              tf.var <- apply(tf.mtx,1,stats::var)              
              target.var <- stats::var(target.mtx)

              # Return only the genes with variances within the var.size of target gene variance
              var.idx <- which(tf.var > (1-var.size)*target.var & tf.var < (1+var.size)*target.var)
              tfs <- names(tf.var[var.idx])
              tf.vars <- tf.var[var.idx]

              return(list(tfs = tfs, tf.vars = tf.vars))
          }
)
#----------------------------------------------------------------------------------------------------
