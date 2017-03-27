#' @title Create a OpenChromFilter object 
#'
#' @description
#' An OpenChromFilter object allows for filtering based on a supplied chromosome and region. Using its
#' associated \code{getCandidates} method, a chromosome, and starting/ending locations for a region, 
#' an OpenChromFilter object can be used to filter a list of possible transcription factors to those
#' that match motifs within the supplied region
#'
#' @include CandidateFilter.R
#' @import methods
#' 
#' @name OpenChromFilter-class
#' @rdname OpenChromFilter-class
#' @aliases OpenChromFilter

#----------------------------------------------------------------------------------------------------
.OpenChromFilter <- setClass("OpenChromFilter", contains = "CandidateFilter")

#----------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#----------------------------------------------------------------------------------------------------
#' @rdname OpenChromFilter-class
#'
#' @param mtx.assay An assay matrix of gene expression data
#' @param quiet A logical denoting whether or not the filter should print output
#'
#' @seealso \code{\link{getCandidates-OpenChromFilter}}
#'
#' @export
#' 
#' @family Filtering Objects
#' 
#' @examples
#' load(system.file(package="TReNA", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
#' openchrom.filter <- OpenChromFilter(mtx.assay = mtx.sub)

OpenChromFilter <- function(mtx.assay=matrix(), quiet=TRUE)
{
    .OpenChromFilter(CandidateFilter(mtx.assay = mtx.assay, quiet = quiet))

} # OpenChromFilter, the constructor
#----------------------------------------------------------------------------------------------------
#' Get candidate genes using the open chromatin filter
#'
#' @aliases getCandidates-OpenChromFilter
#'
#' @usage
#' getCandidates(obj, extraArgs)
#' 
#' @param obj An object of class FootprintFilter
#' @param extraArgs
#' \itemize{
#' \item{"chromosome" A chromosome of interest that contains the regions to be used for filtering}
#' \item{"start" An integer denoting the starting point of the region of interest}
#' \item{"end" An integer denoting the ending point of the region of interest}
#' }
#'
#' @seealso \code{\link{OpenChromFilter}}
#' 
#' @family getCandidate Methods
#' 
#' @return A vector containing all genes with motifs in the supplied region
#'
#' @examples
#'
#' # Use open chromatin filter for MEF2C to filter candidates
#' # in the included Alzheimer's dataset
#' load(system.file(package="TReNA", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
#' openchrom.filter <- OpenChromFilter(mtx.assay = mtx.sub)
#' tfs <- getCandidates(openchrom.filter, extraArgs = list("chromosome" = "chr5",
#' "start" = 888936629, "end" = 898936629))

setMethod("getCandidates", "OpenChromFilter",

    function(obj,extraArgs){

        # Collect arguments from extraArgs
        chromosome <- extraArgs[["chromosome"]]
        start <- extraArgs[["start"]]
        end <- extraArgs[["end"]]

        # Connect to the UCSC MySQL hg38 database
        driver <- RMySQL::MySQL()
        host <- "genome-mysql.cse.ucsc.edu"
        user <- "genome"
        dbname <- "hg38"
        ucsc.db <- DBI::dbConnect(driver, user = user, host = host, dbname = dbname)

        # Pull out the regions corresponding to the region in the ENCODE
        query <- paste("select chrom, chromStart, chromEnd from wgEncodeRegDnaseClustered where",
                       sprintf("chrom = '%s'", chromosome),
                       sprintf("and chromStart >= %d", start),
                       sprintf("and chromEnd <= %d", end),
                       collapse = " ")
        regions <- DBI::dbGetQuery(ucsc.db, query)
        DBI::dbDisconnect(ucsc.db)

        # Convert the regions to motifs using FIMO
        
        # Convert the motifs to TFs and return those
        tbl.out <- regions ## Hack for now

        # Intersect the TFs with the rows in the matrix
        candidate.tfs <- intersect(tbl.out$tf, rownames(obj@mtx.assay))

	# Return the TFs
	return(candidate.tfs)
	}
)
#----------------------------------------------------------------------------------------------------
