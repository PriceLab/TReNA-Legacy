#' @title Create a FootprintFilter object
#'
#' @description
#' A FootprintFilter object allows for filtering based on gene footprinting databases. Using its
#' associated \code{getCandidates} method and URIs for both a genome database and project database,
#' a FootprintFilter object can be used to filter a list of possible transcription factors to those
#' that match footprint motifs for a supplied target gene.
#'
#' @include CandidateFilter.R
#' @import methods
#'
#' @name FootprintFilter-class
#' @rdname FootprintFilter-class
#' @aliases FootprintFilter

#----------------------------------------------------------------------------------------------------
.FootprintFilter <- setClass("FootprintFilter", contains = "CandidateFilter",
                             slots=list(genomeDB="character",
                                        footprintDB="character",
                                        regions="character",
                                        footprintFinder="FootprintFinder")
                             )

#----------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#----------------------------------------------------------------------------------------------------
#' @rdname FootprintFilter-class
#'
#' #' @param quiet A logical denoting whether or not the filter should print output
#'
#' @seealso \code{\link{getCandidates-FootprintFilter}}
#'
#' @export
#'
#' @family Filtering Objects
#'
#' @examples
#' load(system.file(package="TReNA", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
#' footprint.filter <- FootprintFilter()

FootprintFilter <- function(genomeDB, footprintDB, geneCenteredSpec=list(),
                            regionsSpec=list(), quiet=TRUE)
{
   fpFinder <- FootprintFinder(genomeDB, footprintDB, quiet=quiet)


   regions <- c();   # one or more chromLoc strings: "chr5:88903257-88905257"

   if(length(geneCenteredSpec) == 3){
      new.region <- with (geneCenteredSpec, getGenePromoterRegion(fpFinder, targetGene, tssUpstream, tssDownstream))
      new.region.chromLocString <- with(new.region, sprintf("%s:%d-%d", chr, start, end))
      regions <- c(regions, new.region.chromLocString)
      }

   if(length(regionsSpec) > 0)
      regions <- c(regions, regionsSpec)

   .FootprintFilter(CandidateFilter(quiet = quiet),
                    genomeDB=genomeDB,
                    footprintDB=footprintDB,
                    regions=regions,
                    footprintFinder=fpFinder)

} # FootprintFilter, the constructor
#----------------------------------------------------------------------------------------------------
#' Get candidate genes using the variance filter
#'
#' @aliases getCandidates-FootprintFilter
#'
#' @param obj An object of class FootprintFilter
#' @param extraArgs A named list containing 5 fields:
#' \itemize{
#' \item{"target.gene" A designated target gene that should be part of the mtx.assay data}
#' \item{"genome.db.uri" A connection to a genome database containing footprint information}
#' \item{"project.db.uri" A connection to a project database containing footprint information}
#' \item{"size.upstream" An integer denoting the distance upstream of the target gene to look for footprints}
#' \item{"size.downstream" An integer denoting the distance downstream of the target gene to look for footprints}
#' }
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
#' tfs <- getCandidates(footprint.filter, extraArgs = list("target.gene" = target.gene,
#' "genome.db.uri" = genome.db.uri, "project.db.uri" = project.db.uri,
#' "size.upstream" = 1000, "size.downstream" = 1000))


setMethod("getCandidates", "FootprintFilter",

     function(obj){

              # Create a FootprintFinder object and find the footprints
        fp <- FootprintFinder(obj@genomeDB, obj@footprintDB, quiet=TRUE)
        tbl.out <- data.frame()

        for(region in obj@regions){
           chromLoc <- .parseChromLocString(region)
           tbl.fp <- try(with(chromLoc, getFootprintsInRegion(fp, chrom, start, end)))
           if(!(class(tbl.fp) == "try-error")){
              tbl.out <- rbind(tbl.out, mapMotifsToTFsMergeIntoTable(fp, tbl.fp))
              }
           else{
             warning("FootprintFinder error with region %s", region)
             closeDatabaseConnections(fp)
             return(NULL)
              }
           } # for region

        closeDatabaseConnections(fp)
                # Intersect the footprints with the rows in the matrix
        candidate.tfs <- sort(unique(unlist(strsplit(tbl.out$tf, ";"))))
        return(list("tfs" = candidate.tfs, "tbl" = tbl.out))
        }) # getCandidates

#----------------------------------------------------------------------------------------------------
.parseChromLocString <- function(chromLocString)
{
   tokens.0 <- strsplit(chromLocString, ":", fixed=TRUE)[[1]]
   stopifnot(length(tokens.0) == 2)
   chrom <- tokens.0[1]

   tokens.1 <- strsplit(tokens.0[2], "-")[[1]]
   stopifnot(length(tokens.1) == 2)
   start <- as.integer(tokens.1[1])
   end <- as.integer(tokens.1[2])

   return(list(chrom=chrom, start=start, end=end))

} # parseChromLocString
#------------------------------------------------------------------------------------------------------------------------
