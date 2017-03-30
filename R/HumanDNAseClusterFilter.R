#' @title Create a HumanDNAseClusterFilter objects
#'
#' @description
#' An HumanDNAseClusterFilter object allows for filtering based on a supplied chromosome and region. Using its
#' associated \code{getCandidates} method, a chromosome, and starting/ending locations for a region,
#' an HumanDNAseClusterFilter object can be used to filter a list of possible transcription factors to those
#' that match motifs within the supplied region
#'
#' @include CandidateFilter.R
#' @import methods
#'
#' @name HumanDNAseClusterFilter-class
#' @rdname HumanDNAseClusterFilter-class
#' @aliases HumanDNAseClusterFilter

#----------------------------------------------------------------------------------------------------
.HumanDNAseClusterFilter <- setClass("HumanDNAseClusterFilter", contains = "CandidateFilter")

#----------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#----------------------------------------------------------------------------------------------------
#' @rdname HumanDNAseClusterFilter-class
#'
#' @param mtx.assay An assay matrix of gene expression data
#' @param quiet A logical denoting whether or not the filter should print output
#'
#' @seealso \code{\link{getCandidates-HumanDNAseClusterFilter}}
#'
#' @export
#'
#' @family Filtering Objects
#'
#' @examples
#' load(system.file(package="TReNA", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
#' openchrom.filter <- HumanDNAseClusterFilter(mtx.assay = mtx.sub)

HumanDNAseClusterFilter <- function(mtx.assay=matrix(), quiet=TRUE)
{
    .HumanDNAseClusterFilter(CandidateFilter(mtx.assay = mtx.assay, quiet = quiet))

} # HumanDNAseClusterFilter, the constructor
#----------------------------------------------------------------------------------------------------
#' Get candidate genes using the open chromatin filter
#'
#' @aliases getCandidates-HumanDNAseClusterFilter
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
#' @seealso \code{\link{HumanDNAseClusterFilter}}
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
#' openchrom.filter <- HumanDNAseClusterFilter(mtx.assay = mtx.sub)
#' tfs <- getCandidates(openchrom.filter, extraArgs = list("chromosome" = "chr5",
#' "start" = 888936629, "end" = 898936629))
#
# "DNASE I Hyptersensitivity Peak Clusters from ENCODE (95 cell types)"

setMethod("getCandidates", "HumanDNAseClusterFilter",

    function(obj, extraArgs){

        region.score.threshold <- 250;     # apparently what ucsc uses for this track
        motif.min.match.percentage <- 95;  # 95% of the maximum possible score, as calcluated by bioc:matchPWM
        stopifnot(all(c("chrom", "start", "end") %in% names(extraArgs)))

          # Collect arguments from extraArgs
        chrom <- extraArgs[["chrom"]]
        start <- extraArgs[["start"]]
        end <- extraArgs[["end"]]

        if("region.score.threshold" %in% names(extraArgs))
            region.score.threshold <- extraArgs$region.score.threshold

        if("motif.min.match.percentage" %in% names(extraArgs))
            motif.min.match.percentage <- extraArgs$motif.min.match.percentage

        tbl.regions <- .getRegions(chrom, start, end, region.score.threshold)
        seqs <- .getSequence(tbl.regions)
        tbl.motifs <- .getScoredMotifs(seqs, motif.min.match.percentage, obj@quiet)
        region.count <- length(seqs)
        tbl.summary <- data.frame()
        for(i in seq_len(region.count)){
            tbl.summary <- rbind(tbl.summary, cbind(tbl.motifs[[i]], tbl.regions[i,]))
            }
        colnames(tbl.summary) <- c("motif.start", "motif.end", "motif.width", "motif.score", "motif", "match",
                                   "strand", "chrom", "regionStart", "regionEnd", "regionScore", "sourceCount")
        tbl.summary <- tbl.summary[, c("chrom", "regionStart", "regionEnd", "regionScore", "sourceCount", "motif",
                                   "match", "motif.start", "motif.end", "motif.width", "motif.score", "strand")]
        tbl.mg <- read.table(system.file(package="TReNA", "extdata", "motifGenes.tsv"), sep="\t", as.is=TRUE,
                             header=TRUE)

        tfs.by.motif <- lapply(tbl.summary$motif, function(m) subset(tbl.mg, motif==m)$tf.gene)
        all.tfs <- sort(unique(unlist(tfs.by.motif)))
        tfs.by.motif.joined <- unlist(lapply(tfs.by.motif, function(m) paste(m, collapse=";")))
        tbl.summary$tf <- tfs.by.motif.joined
        result <- list(tbl=tbl.summary,
                       tfs=all.tfs)

    }) # getCandidates

#----------------------------------------------------------------------------------------------------
.getRegions <- function(chromosome, start, end, score.threshold=200)
{
   driver <- RMySQL::MySQL()
   host <- "genome-mysql.cse.ucsc.edu"
   user <- "genome"
   dbname <- "hg38"

   ucsc.db <- DBI::dbConnect(driver, user = user, host = host, dbname = dbname)

   # Pull out the regions corresponding to the region in the ENCODE
   query <- paste("select chrom, chromStart, chromEnd, score, sourceCount from wgEncodeRegDnaseClustered where",
                  sprintf("chrom = '%s'", chromosome),
                  sprintf("and chromStart >= %d", start),
                  sprintf("and chromEnd <= %d", end),
                  collapse = " ")

   suppressWarnings(  # MySQL returns unsigned integers.  do not wish to see the conversion warnings
     tbl.regions <- DBI::dbGetQuery(ucsc.db, query)
     )

   DBI::dbDisconnect(ucsc.db)

   subset(tbl.regions, score >= score.threshold)

} # .getRegions
#----------------------------------------------------------------------------------------------------
.getSequence <- function(tbl.regions)
{
   if(!exists("reference.genome")){   # move to constructor
      library(BSgenome.Hsapiens.UCSC.hg38)
      reference.genome <<- BSgenome.Hsapiens.UCSC.hg38
      }

   gr.regions <- with(tbl.regions, GRanges(seqnames=chrom, IRanges(start=chromStart, end=chromEnd)))
   seqs <- getSeq(reference.genome, gr.regions)
   as.character(seqs)

} # .getSequence
#----------------------------------------------------------------------------------------------------
.findMotifs <- function(sequence, pfms, min.match.percentage=95, quiet=TRUE)
{
   min.match.as.string <- sprintf("%02d%%", min.match.percentage)

   search <- function(motifName, mtx, seq){
      #browser()
      hits.fwd <- matchPWM(mtx, seq, with.score=TRUE, min.score=min.match.as.string)
      seq.revcomp <- as.character(reverseComplement(DNAString(seq)))
      hits.rev <- matchPWM(mtx, seq.revcomp, with.score=TRUE, min.score=min.match.as.string)
      tbl <- data.frame()
      if(length(hits.fwd) > 0){
          if(!quiet) printf("%d +", length(hits.fwd))
          match <- substring(as.character(subject(hits.fwd)), start(ranges(hits.fwd)), end(ranges(hits.fwd)))
          tbl <- data.frame(ranges(hits.fwd), score=mcols(hits.fwd)$score, motif=motifName, match=match, strand="+")
          }
      if(length(hits.rev) > 0){
          if(!quiet) printf("%d -", length(hits.rev))
          match <- substring(as.character(subject(hits.rev)), start(ranges(hits.rev)), end(ranges(hits.rev)))
          tbl.rev <- data.frame(ranges(hits.rev), score=mcols(hits.rev)$score, motif=motifName, match=match, strand="-")
          tbl <- rbind(tbl, tbl.rev)
          }
       return(tbl)
       }

    count <- length(pfms)
    xx <- lapply(1:count, function(i) search(names(pfms)[i], pfms[[i]], sequence))
    tbl.result <- do.call("rbind", xx)
    if(nrow(tbl.result) == 0){
      return(data.frame())
      }
    else{
      tbl.result$motif <- as.character(tbl.result$motif)
      #tbl.result$seq <- as.character(tbl.result$seq)
      tbl.result$match <- as.character(tbl.result$match)
      tbl.result$strand <- as.character(tbl.result$strand)
      return(tbl.result[order(tbl.result$score,decreasing=TRUE),])
      }

}  # .findMotifs
#------------------------------------------------------------------------------------------------------------------------
.getScoredMotifs <- function(seqs, min.match.percentage=95, quiet=TRUE)
{
     parseLine <- function(textOfLine) {
        # first delete the leading A, C, G or T.  then the square brackets.  then convert
        x <- substr(textOfLine, 2, nchar(textOfLine))
        x2 <- sub(" *\\[ *", "", x)
        x3 <- sub(" *\\] *", "", x2)
        counts <- as.integer(strsplit(x3, "\\s+", perl=TRUE)[[1]])
        return(counts)
        } # parseLine

    parseJasparPwm = function (lines) {
      stopifnot(length(lines)==5) # title line, one line for each base
      motif.name.raw = strsplit(lines[1], "\t")[[1]][1]
      motif.name <- gsub(">", "", motif.name.raw, fixed=TRUE)
        # expect 4 rows, and a number of columns we can discern from  the incoming text.
      a.counts <- parseLine(lines[[2]])
      c.counts <- parseLine(lines[[3]])
      g.counts <- parseLine(lines[[4]])
      t.counts <- parseLine(lines[[5]])
      stopifnot(length(a.counts) == length(c.counts))
      stopifnot(length(a.counts) == length(g.counts))
      stopifnot(length(a.counts) == length(t.counts))
      cols <- length(a.counts)
      mtx <- matrix (nrow=4, ncol=cols, dimnames=list(c('A','C','G','T'), as.character(1:cols)))
      mtx[1,] <- a.counts
      mtx[2,] <- c.counts
      mtx[3,] <- g.counts
      mtx[4,] <- t.counts
      return(list(title=motif.name, matrix=mtx))
      } # parsePwm

   readRawJasparMatrices = function (uri) {
     all.lines <- scan(uri, what=character(0), sep='\n', quiet=TRUE)
     title.lines <- grep ('^>', all.lines)
     title.line.count <- length (title.lines)
     max <- title.line.count - 1
     pwms <- list()
     for(i in 1:max){
       start.line <- title.lines [i]
       end.line <- title.lines [i+1] - 1
       x <- parseJasparPwm (all.lines [start.line:end.line])
       pwms[[i]] <- list(title=x$title, matrix=x$matrix)
       } # for i
     invisible (pwms)
     } # readRawJasparMatrices

   if(!exists("pfms")){
      uri <- "http://jaspar.genereg.net/html/DOWNLOAD/JASPAR_CORE/pfm/nonredundant/pfm_vertebrates.txt"
       x <- readRawJasparMatrices(uri)
         # normalize, so that a frequency sum of 1.0 is true across the 4 possible bases at each position
       pfms <<- lapply(x, function(e) apply(e$matrix, 2, function(col) col/sum(col)))
       names(pfms) <<- as.character(lapply(x, function(e) e$title))
       }

   result <- lapply(seqs, function(seq) .findMotifs(seq, pfms, min.match.percentage, quiet))
   result


} # .getScoredMotifs
#----------------------------------------------------------------------------------------------------
