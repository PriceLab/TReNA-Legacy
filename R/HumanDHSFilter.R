#----------------------------------------------------------------------------------------------------
.HumanDNAseClusterFilter <- setClass("HumanDNAseClusterFilter",
                                     contains="CandidateFilter",
                                     representation (pfms='list'))

#----------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#----------------------------------------------------------------------------------------------------
setGeneric('getEncodeRegulatoryTableNames', signature='obj', function(obj) standardGeneric ('getEncodeRegulatoryTableNames'))
setGeneric('getRegulatoryRegions', signature='obj',
           function(obj, tableName, chromosome, start, end, score.threshold=200, quiet=TRUE)
               standardGeneric ('getRegulatoryRegions'))
#----------------------------------------------------------------------------------------------------
HumanDNAseClusterFilter <- function(mtx.assay=matrix(), quiet=TRUE)
{
    uri <- "http://jaspar.genereg.net/html/DOWNLOAD/JASPAR_CORE/pfm/nonredundant/pfm_vertebrates.txt"
    x <- .readRawJasparMatrices(uri)
         # normalize, so that a frequency sum of 1.0 is true across the 4 possible bases at each position
    pfms <- lapply(x, function(e) apply(e$matrix, 2, function(col) col/sum(col)))
    names(pfms) <- as.character(lapply(x, function(e) e$title))

    .HumanDNAseClusterFilter(CandidateFilter(mtx.assay = mtx.assay, quiet = quiet), pfms=pfms)

} # HumanDNAseClusterFilter, the constructor
#----------------------------------------------------------------------------------------------------
setMethod("getEncodeRegulatoryTableNames", "HumanDNAseClusterFilter",

     function(obj){
       driver <- RMySQL::MySQL()
       host <- "genome-mysql.cse.ucsc.edu"
       user <- "genome"
       dbname <- "hg38"
       db <- DBI::dbConnect(driver, user = user, host = host, dbname = dbname)
       all.tableNames <- DBI::dbListTables(db);
          # manual check (13 apr 2017) shows that only wgEncodeReg Peak tabel
          # and the "wgEncodeRegDnaseClustered" table, have scored chromosomal regions in them
       tableNames <- grep("wgEncodeReg.*Peak$", DBI::dbListTables(db), value=TRUE)
       clusteredTable <- "wgEncodeRegDnaseClustered"
       if(clusteredTable %in% all.tableNames)
          tableNames <- c(clusteredTable, tableNames)
       lapply(dbListConnections(driver), DBI::dbDisconnect)
       invisible(tableNames)
       })

#----------------------------------------------------------------------------------------------------
setMethod("getCandidates", "HumanDNAseClusterFilter",

    function(obj, extraArgs){

          # Collect arguments from extraArgs
        stopifnot(all(c("chrom", "start", "end", "region.score.threshold", "motif.min.match.percentage",
                        "tableName") %in% names(extraArgs)))

        chrom <- extraArgs[["chrom"]]
        start <- extraArgs[["start"]]
        end <- extraArgs[["end"]]
        tableName <- extraArgs[["tableName"]]
        region.score.threshold <- extraArgs$region.score.threshold
        motif.min.match.percentage <- extraArgs$motif.min.match.percentage

        tbl.regions <- getRegulatoryRegions(obj, tableName, chrom, start, end)
        tbl.regions <- subset(tbl.regions, score >= region.score.threshold)
        seqs <- .getSequence(tbl.regions)
        tbl.motifs <- .getScoredMotifs(seqs, motif.min.match.percentage, obj@quiet)
        region.count <- length(seqs)
        tbl.summary <- data.frame()
        all.tfs <- c()
        for(i in seq_len(region.count)){
            #printf("HumanDNAseClusterFilter, line 100, combining tbl.motifs and tbl.regions into tbl.summary")
            #browser();
            if(nrow(tbl.motifs[[i]]) > 0)
               tbl.summary <- rbind(tbl.summary, cbind(tbl.motifs[[i]], tbl.regions[i,]))
            }
        if(nrow(tbl.summary) > 0){
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
           }
        list(tbl=tbl.summary,tfs=all.tfs)
	}) # getCandidates

#----------------------------------------------------------------------------------------------------
setMethod("getRegulatoryRegions", "HumanDNAseClusterFilter",

    function(obj, tableName, chromosome, start, end, score.threshold=200, quiet=TRUE) {

       driver <- RMySQL::MySQL()
       host <- "genome-mysql.cse.ucsc.edu"
       user <- "genome"
       dbname <- "hg38"

       db <- dbConnect(driver, user = user, host = host, dbname = dbname)

       #schema <- colnames(dbGetQuery(db, sprintf("select * from %s limit 1", tableName)))
       #suppressWarnings(dbGetQuery(db, sprintf("select * from %s limit 3", tableName)))

       #if(!all(c("chrom", "chromStart", "chromEnd") %in% schema)){
       #   printf("%s lacks chrom, start, end", tableName)
       #   lapply(dbListConnections(driver), dbDisconnect)
       #   return(data.frame())
       #   }

       #main.clause <- sprintf("select chrom, chromStart, chromEnd, score, sourceCount from %s where", tableName);
       main.clause <- sprintf("select * from %s where", tableName);

      # Pull out the regions corresponding to the region in ENCODE
       query <- paste(main.clause,
                      sprintf("chrom = '%s'", chromosome),
                      sprintf("and chromStart >= %d", start),
                      sprintf("and chromEnd <= %d", end),
                      collapse = " ")

         # handle the usual case first: a start:end region many times larger than a typical DHS region
       suppressWarnings(  # MySQL returns unsigned integers.  hide these unproblematic conversion warnings
         tbl.regions <- dbGetQuery(db, query)
          )

       if(!quiet)
           printf("%d DHS regions reported in %d bases, start:end unmodified", nrow(tbl.regions), 1 + end - start)

      # if no hits, then perhaps a very small region is requested, one which falls entirely within a DHS region
       if(nrow(tbl.regions) == 0) {
          extension <- 10000
          if(!quiet)
             printf("possible that start:end (%d) is small relative to DHS regions, extend by %d", (1 + end - start),
                   extension);
          query <- paste(main.clause,
                         sprintf("chrom = '%s'", chromosome),
                         sprintf("and chromStart >= %d", start - extension),
                         sprintf("and chromEnd   <= %d", end + extension),
                         collapse = " ")
           suppressWarnings(  # MySQL returns unsigned integers.  hide these unproblematic conversion warnings
              tbl.regionsExtended <- dbGetQuery(db, query)
              )
         if(nrow(tbl.regionsExtended) > 0) { # now find just the intersection of DHS and requested region
            if(!quiet)
               printf("tbl.regionsExtended: %d rows, now do intersect", nrow(tbl.regionsExtended));
             gr.regions <- with(tbl.regionsExtended, GRanges(seqnames=chromosome, IRanges(chromStart, chromEnd)))
             gr.target <- GRanges(seqnames=chromosome, IRanges(start, end))
             gr.intersect <- GenomicRanges::intersect(gr.target, gr.regions, ignore.strand=TRUE)
             if(length(gr.intersect) == 1){
                region.index <- subjectHits(findOverlaps(gr.target, gr.regions))
                tbl.regions <- tbl.regionsExtended[region.index,]
                tbl.regions$chromStart[1] <- max(tbl.regions$chromStart[1], start)
                tbl.regions$chromEnd[1] <- min(tbl.regions$chromEnd[1], end)
                } # one region from the extended query intersects with the requested start:end.
             } # small region query, within DHS region
          } # small region, extension search

   lapply(dbListConnections(driver), dbDisconnect)

   invisible(tbl.regions[, c("chrom", "chromStart", "chromEnd", "name", "score")])

   }) # getRegulatoryRegions

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
.matchForwardAndReverse <- function(sequence, pfm, motifName, min.match.percentage=95, quiet=TRUE)
{
   min.match.as.string <- sprintf("%02d%%", min.match.percentage)

   hits.fwd <- matchPWM(pfm, sequence, with.score=TRUE, min.score=min.match.as.string)
   seq.revcomp <- as.character(reverseComplement(DNAString(sequence)))
   hits.rev <- matchPWM(pfm, seq.revcomp, with.score=TRUE, min.score=min.match.as.string)

   tbl <- data.frame()
   if(length(hits.fwd) > 0){
      if(!quiet) printf("%d +", length(hits.fwd))
      match <- substring(as.character(subject(hits.fwd)), start(ranges(hits.fwd)), end(ranges(hits.fwd)))
      tbl <- data.frame(ranges(hits.fwd), score=mcols(hits.fwd)$score, motif=motifName, match=match, strand="+",
                        stringsAsFactors=FALSE)
      }

   if(length(hits.rev) > 0){
      if(!quiet) printf("%d -", length(hits.rev))
      match <- substring(as.character(subject(hits.rev)), start(ranges(hits.rev)), end(ranges(hits.rev)))
      tbl.rev <- data.frame(ranges(hits.rev), score=mcols(hits.rev)$score, motif=motifName, match=match, strand="-",
                            stringsAsFactors=FALSE)
         # transform the start/end so that they are forward-strand relative
      true.start <- 1 + nchar(sequence) - tbl.rev$end
      true.end   <- 1 + nchar(sequence) - tbl.rev$start
      tbl.rev$start <- true.start
      tbl.rev$end   <- true.end
      tbl <- rbind(tbl, tbl.rev)
      }

   tbl

} # .matchForwardAndReverse
#----------------------------------------------------------------------------------------------------
.findMotifs <- function(sequence, pfms, min.match.percentage=95, quiet=TRUE)
{
   min.match.as.string <- sprintf("%02d%%", min.match.percentage)

   search <- function(motifName, mtx, seq){
      hits.fwd <- matchPWM(mtx, seq, with.score=TRUE, min.score=min.match.as.string)
      seq.revcomp <- as.character(reverseComplement(DNAString(seq)))
      hits.rev <- matchPWM(mtx, seq.revcomp, with.score=TRUE, min.score=min.match.as.string)
      #browser()
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
    xx <- lapply(1:count, function(i) .matchForwardAndReverse(sequence, pfms[[i]], names(pfms)[i], min.match.percentage, quiet))
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
.parseLine <- function(textOfLine) {
   # first delete the leading A, C, G or T.  then the square brackets.  then convert
   x <- substr(textOfLine, 2, nchar(textOfLine))
   x2 <- sub(" *\\[ *", "", x)
   x3 <- sub(" *\\] *", "", x2)
   counts <- as.integer(strsplit(x3, "\\s+", perl=TRUE)[[1]])
   return(counts)
   } # parseLine

 .parseJasparPwm = function (lines) {
   stopifnot(length(lines)==5) # title line, one line for each base
   motif.name.raw = strsplit(lines[1], "\t")[[1]][1]
   motif.name <- gsub(">", "", motif.name.raw, fixed=TRUE)
     # expect 4 rows, and a number of columns we can discern from  the incoming text.
   a.counts <- .parseLine(lines[[2]])
   c.counts <- .parseLine(lines[[3]])
   g.counts <- .parseLine(lines[[4]])
   t.counts <- .parseLine(lines[[5]])
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

.readRawJasparMatrices = function (uri) {
  all.lines <- scan(uri, what=character(0), sep='\n', quiet=TRUE)
  title.lines <- grep ('^>', all.lines)
  title.line.count <- length (title.lines)
  max <- title.line.count - 1
  pwms <- list()
  for(i in 1:max){
    start.line <- title.lines [i]
    end.line <- title.lines [i+1] - 1
    x <- .parseJasparPwm (all.lines [start.line:end.line])
    pwms[[i]] <- list(title=x$title, matrix=x$matrix)
    } # for i
  invisible (pwms)
  } # readRawJasparMatrices


