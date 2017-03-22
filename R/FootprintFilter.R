#' @title Create a FootprintFilter object 
#'
#' @description
#' A FootprintFilter object allows for filtering based on gene footprinting databases. Using its
#' associated \code{getCandidates} method and URIs for both a genome database and project database,
#' a FootprintFilter object can be used to filter a list of possible transcription factors to those
#' that match footprint motifs for a supplied target gene.
#'
#' @export
#' @include CandidateFilter.R
#' @name FootprintFilter-class
#' @rdname FootprintFilter-class
#' @aliases FootprintFilter

#----------------------------------------------------------------------------------------------------
.FootprintFilter <- setClass("FootprintFilter", contains = "CandidateFilter")

#----------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#----------------------------------------------------------------------------------------------------
#' @rdname FootprintFilter-class
#'
#' @param mtx.assay An assay matrix of gene expression data
#' @param quiet A logical denoting whether or not the filter should print output
#'
#' @seealso \code{\link{getCandidates-FootprintFilter}}
#'
#' @family Filtering Objects
#' 
#' @examples
#' load(system.file(package="TReNA", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
#' footprint.filter <- FootprintFilter(mtx.assay = mtx.sub)

FootprintFilter <- function(mtx.assay=matrix(), quiet=TRUE)
{
    .FootprintFilter(CandidateFilter(mtx.assay = mtx.assay, quiet = quiet))

} # FootprintFilter, the constructor
#----------------------------------------------------------------------------------------------------
#' Get candidate genes using the variance filter
#'
#' @aliases getCandidates-FootprintFilter
#'
#' @usage
#' footprint.filter <- FootprintFilter(mtx.assay, quiet)
#' tfs <- getCandidates(footprint.filter, target.gene, genome.db.uri, projct.db.uri)
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
#' @seealso \code{\link{FootprintFilter}}
#' 
#' @family getCandidate Methods
#' 
#' @return A vector containing all genes with variances less than the target gene
#'
#' @examples
#'
#' # Use footprint filter with the included SQLite database for MEF2C to filter candidates
#' # in the included Alzheimer's dataset
#' load(system.file(package="TReNA", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
#' footprint.filter <- FootprintFilter(mtx.assay = mtx.sub)
#'
#' target.gene <- "MEF2C"
#' db.address <- system.file(package="TReNA", "extdata")
#' genome.db.uri <- paste("sqlite:/",db.address,"genome.sub.db", sep = "/")
#' project.db.uri <- paste("sqlite:/",db.address,"project.sub.db", sep = "/")
#' 
#' tfs <- getCandidates(footprint.filter, target.gene, genome.db.uri, project.db.uri)
    

setMethod("getCandidates", "FootprintFilter",

    function(obj,target.gene, genome.db.uri, project.db.uri, size.upstream=1000, size.downstream=1000){
        # Create a FootprintFinder object and find the footprints
        fp <- FootprintFinder(genome.db.uri, project.db.uri, quiet=TRUE)
        tbl.fp <- getFootprintsForGene(fp, target.gene,
                                       size.upstream=size.upstream, size.downstream=size.downstream)

        # Convert the footprints to genes and close the database connection
        tbl.out <- mapMotifsToTFsMergeIntoTable(fp, tbl.fp)
        closeDatabaseConnections(fp)

        # Intersect the footprints with the rows in the matrix
        candidate.tfs <- intersect(tbl.out$tf, rownames(obj@mtx.assay))

	# Return the TFs
	return(candidate.tfs)
	}
)
#----------------------------------------------------------------------------------------------------
