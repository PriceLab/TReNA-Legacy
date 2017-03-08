#' Filter based on gene footprints
#'
#' @include CandidateFilter.R
#' @name FootprintFilter-class
#' @param mtx.assay An assay matrix of gene expression data

#----------------------------------------------------------------------------------------------------
.FootprintFilter <- setClass("FootprintFilter", contains = "CandidateFilter")

#----------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#----------------------------------------------------------------------------------------------------
FootprintFilter <- function(mtx.assay=matrix(), quiet=TRUE)
{
    .FootprintFilter(CandidateFilter(mtx.assay = mtx.assay, quiet = quiet))

} # FootprintFilter, the constructor
#----------------------------------------------------------------------------------------------------
#' Get candidate genes using the variance filter
#'
#' @aliases getCandidates-FootprintFilter
#'
#' @param obj An object of class FootprintFilter
#' @param target.gene A designated target gene that should be part of the mtx.assay data
#' @param genome.db.uri A connection to a genome database containing footprint information
#' @param project.db.uri A connection to a project database containing footprint information
#' @param size.upstream An integer denoting the distance upstream of the target gene to look for footprints
#' (default = 1000)
#' @param size.downstream An integer denoting the distance downstream of the target gene to look for footprints
#' (default = 1000)
#'
#' @return A vector containing all genes with variances less than the target gene

setMethod("getCandidates", "FootprintFilter",

    function(obj,target.gene, genome.db.uri, project.db.uri, size.upstream=1000, size.downstream=1000){
        # Create a FootprintFinder object and find the footprints
        fp <- FootprintFinder(genome.db.uri, project.db.uri, quiet=TRUE)
        tbl.fp <- getFootprintsForGene(fp, target.gene,
                                       size.upstream=size.upstream, size.downstream=size.downstream)

        # Convert the footprints to genes and close the database connection
        tbl.out <- mapMotifsToTFsMergeIntoTable(fp, tbl.fp)
        closeDatabaseConnections(fp)

	# Return the TFs
	return(tbl.out$tf)
	}
)
#----------------------------------------------------------------------------------------------------
