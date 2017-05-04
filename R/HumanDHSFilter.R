#----------------------------------------------------------------------------------------------------
.HumanDHSFilter <- setClass("HumanDHSFilter",
                            contains="CandidateFilter",
                            slots=list(genomeName="character",
                                       genome="BSgenome",
                                       encodeTableName="character",
                                       fimoDB="DBIConnection",
                                       geneInfoDatabase.uri="character",   # access to gtf database
                                       geneCenteredSpec="list",
                                       regionsSpec="character",
                                       pfms="list",
                                       quiet="logical"))

#----------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#----------------------------------------------------------------------------------------------------
setGeneric("getEncodeRegulatoryTableNames", signature="obj", function(obj) standardGeneric ("getEncodeRegulatoryTableNames"))
setGeneric("getRegulatoryRegions", signature="obj",
           function(obj, encode.table.name, chromosome, start, end, score.threshold=200, quiet=TRUE)
               standardGeneric ("getRegulatoryRegions"))
setGeneric("getSequence", signature="obj", function(obj, tbl.regions) standardGeneric ("getSequence"))
setGeneric("geneSymbolToTSS", signature="obj", function(obj, geneSymbol) standardGeneric("geneSymbolToTSS"))
#----------------------------------------------------------------------------------------------------
HumanDHSFilter <- function(genomeName,
                           encodeTableName="wgEncodeRegDnaseClustered",
                           fimoDatabase.uri,
                           geneInfoDatabase.uri,
                           geneCenteredSpec=c(),
                           regionsSpec=c(),
                           quiet=TRUE)
{
   regions <- c();   # one or more chromLoc strings: "chr5:88903257-88905257"


    uri <- "http://jaspar.genereg.net/html/DOWNLOAD/JASPAR_CORE/pfm/nonredundant/pfm_vertebrates.txt"
    x <- .readRawJasparMatrices(uri)
         # normalize, so that a frequency sum of 1.0 is true across the 4 possible bases at each position
    pfms <- lapply(x, function(e) apply(e$matrix, 2, function(col) col/sum(col)))
    names(pfms) <- as.character(lapply(x, function(e) e$title))

    if(genomeName == "hg38"){
       library(BSgenome.Hsapiens.UCSC.hg38)
       reference.genome <- BSgenome.Hsapiens.UCSC.hg38
       }
    else if(genomeName == "hg19"){
       library(BSgenome.Hsapiens.UCSC.hg19)
       reference.genome <- BSgenome.Hsapiens.UCSC.hg19
       }
    else {
      stop(sprintf("HumanDHSFilter genome.name not in hg19, hg38: '%s'", genomeName))
       }


   fimo.db.info <- .parseDatabaseUri(fimoDatabase.uri)

   stopifnot(fimo.db.info$brand %in% c("postgres"))
   switch(fimo.db.info$brand,
          postgres={
             host <- fimo.db.info$host
             dbname <- fimo.db.info$name
             driver <- RPostgreSQL::PostgreSQL()
             fimo.db <- DBI::dbConnect(driver, user= "trena", password="trena", dbname=dbname, host=host)
             print(dbListTables(fimo.db))
             },
          printf("unrecognized fimo db protoocol '%s'", fimo.db.info$brand)
          )
   .HumanDHSFilter(CandidateFilter(quiet = quiet),
                   genomeName=genomeName,
                   fimoDB=fimo.db,
                   encodeTableName=encodeTableName,
                   geneInfoDatabase.uri=geneInfoDatabase.uri,
                   genome=reference.genome,
                   geneCenteredSpec=geneCenteredSpec,
                   regionsSpec=regionsSpec,
                   pfms=pfms,
                   quiet=quiet)

} # HumanDHSFilter, the constructor
#----------------------------------------------------------------------------------------------------
setMethod("getEncodeRegulatoryTableNames", "HumanDHSFilter",

     function(obj){
       driver <- RMySQL::MySQL()
       host <- "genome-mysql.cse.ucsc.edu"
       user <- "genome"
       dbname <- obj@genomeName
       db <- DBI::dbConnect(driver, user = user, host = host, dbname = dbname)
       all.tableNames <- DBI::dbListTables(db);
          # manual check (13 apr 2017) shows that only wgEncodeReg Peak tabel
          # and the "wgEncodeRegDnaseClustered" table, have scored chromosomal regions in them
       tableNames <- grep("wgEncodeReg.*Peak$", DBI::dbListTables(db), value=TRUE)
       clusteredTable <- switch(obj@genomeName,
                                hg19="wgEncodeRegDnaseClusteredV3",
                                hg38="wgEncodeRegDnaseClustered",
                                NA)
       if(clusteredTable %in% all.tableNames)
          tableNames <- c(clusteredTable, tableNames)
       lapply(dbListConnections(driver), DBI::dbDisconnect)
       invisible(tableNames)
       })

#----------------------------------------------------------------------------------------------------
setMethod("show", "HumanDHSFilter",

     function(object){
        s <- sprintf("HumanDHSFilter...")
        cat(s, sep="\n")
        })
#----------------------------------------------------------------------------------------------------
setMethod("geneSymbolToTSS", "HumanDHSFilter",

     function(obj, geneSymbol){
        geneInfo.db.info <- .parseDatabaseUri(obj@geneInfoDatabase.uri)
        host <- geneInfo.db.info$host
        dbname <- geneInfo.db.info$name
        driver <- RPostgreSQL::PostgreSQL()
        db.geneInfo <- DBI::dbConnect(driver, user= "trena", password="trena", dbname=dbname, host=host)
        print(dbListTables(db.geneInfo))
        #db.gtf <- dbConnect(PostgreSQL(), user= "trena", password="trena", dbname="gtf", host="whovian")
        query <- sprintf("select * from hg38human where moleculetype='gene' and gene_biotype='protein_coding' and gene_name='%s'",
                         obj@geneCenteredSpec$targetGene)
        tbl <- dbGetQuery(db.geneInfo, query);
        tss <- tbl$start[1];
        chrom <- tbl$chr
        DBI::dbDisconnect(db.geneInfo)
        list(chrom=chrom, tss=tss)
        })

#----------------------------------------------------------------------------------------------------
setMethod("getCandidates", "HumanDHSFilter",

    function(obj){

       regions <- c();   # one or more chromLoc strings: "chr5:88903257-88905257"

       if(length(obj@geneCenteredSpec) == 3){
          expected.fields <- c("targetGene", "tssUpstream", "tssDownstream")
          stopifnot(all(expected.fields %in% names(obj@geneCenteredSpec)))
          x <- geneSymbolToTSS(obj)
          start <- x$tss - obj@geneCenteredSpec$tssUpstream
          end <- x$tss + obj@geneCenteredSpec$tssDownstream
          new.region.chromLocString <- sprintf("%s:%d-%d", x$chrom, start, end)
          regions <- c(regions, new.region.chromLocString)
          }

       if(!(all(is.na(obj@regionsSpec)))){
          regions <- c(regions, obj@regionsSpec)
          }

       tbl.regions <- data.frame()

       for(region in regions){
           chromLoc <- .parseChromLocString(region)
           chrom <- chromLoc$chrom
           start <- chromLoc$start
           end   <- chromLoc$end
           encode.table.name <- obj@encodeTableName
           #region.score.threshold <- obj@region.score.threshold
           #motif.min.match.percentage <- obj@motif.min.match.percentage
           tbl.regions <- rbind(tbl.regions, getRegulatoryRegions(obj, encode.table.name, chrom, start, end))

           #browser()
           ##tbl.regions <- subset(tbl.regions, score >= region.score.threshold)
           #seqs <- getSequence(obj, tbl.regions)
           #tbl.motifs <- .getScoredMotifs(seqs, 75, obj@quiet)
           #region.count <- length(seqs)
           #tbl.summary <- data.frame()
           #all.tfs <- c()
           #for(i in seq_len(region.count)){
           #   #printf("HumanDHSFilter, line 100, combining tbl.motifs and tbl.regions into tbl.summary")
           #   #browser();
           #   if(nrow(tbl.motifs[[i]]) > 0)
           #      tbl.summary <- rbind(tbl.summary, cbind(tbl.motifs[[i]], tbl.regions[i,]))
           #   } # for i
           } # for region

         if(nrow(tbl.regions) > 0){
            start.pos <- min(tbl.regions$chromStart) - 100
            end.pos   <- max(tbl.regions$chromEnd) + 100
            fimo.chrom <- sub("chr", "", unique(tbl.regions$chrom))
            query <- sprintf("select * from fimo_hg38 where chrom='%s' and start >= %d and endpos <= %d",
                             fimo.chrom, start.pos, end.pos)
            tbl.fimo <- dbGetQuery(obj@fimoDB, query)
            gr.fimo <- GRanges(seqnames=paste("chr", tbl.fimo$chrom, sep=""), IRanges(start=tbl.fimo$start, end=tbl.fimo$endpos))
            gr.regions <- GRanges(seqnames=tbl.regions$chrom, IRanges(start=tbl.regions$chromStart, end=tbl.regions$chromEnd))
            tbl.ov <- as.data.frame(findOverlaps(gr.fimo, gr.regions))
            tbl.out <- cbind(tbl.fimo[tbl.ov$queryHits,], tbl.regions[tbl.ov$subjectHits,])
            tbl.mg <- read.table(system.file(package="TReNA", "extdata", "motifGenes.tsv"), sep="\t", as.is=TRUE, header=TRUE)
            tfs.by.motif <- lapply(tbl.out$motifname, function(m) subset(tbl.mg, motif==m)$tf.gene)
            all.tfs <- sort(unique(unlist(tfs.by.motif)))
            tfs.by.motif.joined <- unlist(lapply(tfs.by.motif, function(m) paste(m, collapse=";")))
            tbl.out$tf <- tfs.by.motif.joined
            }

         # if(nrow(tbl.summary) > 0){
         #  colnames(tbl.summary) <- c("motif.start", "motif.end", "motif.width", "motif.score",
         #                             "motif.maxScore", "motif.relativeScore", "motif", "match",
         #                             "strand", "chrom", "regionStart", "regionEnd",  "sourceCount", "regionScore")
#
#           colnames.in.order <- c("chrom", "regionStart", "regionEnd", "regionScore", "sourceCount",
#                                          "motif", "match", "strand", "motif.start", "motif.end", "motif.width",
#                                          "motif.score", "motif.maxScore", "motif.relativeScore")
#           tbl.summary <- tbl.summary[, colnames.in.order]
#           tbl.mg <- read.table(system.file(package="TReNA", "extdata", "motifGenes.tsv"), sep="\t", as.is=TRUE, header=TRUE)
#           tfs.by.motif <- lapply(tbl.summary$motif, function(m) subset(tbl.mg, motif==m)$tf.gene)
#           all.tfs <- sort(unique(unlist(tfs.by.motif)))
#           tfs.by.motif.joined <- unlist(lapply(tfs.by.motif, function(m) paste(m, collapse=";")))
#           tbl.summary$tf <- tfs.by.motif.joined
#           tbl.summary$motif.start <- -1 + tbl.summary$regionStart + tbl.summary$motif.start
#           tbl.summary$motif.end   <- -1 + tbl.summary$regionStart + tbl.summary$motif.end
#           }

        list(tbl=tbl.out,tfs=all.tfs)

	}) # getCandidates

#----------------------------------------------------------------------------------------------------
setMethod("getRegulatoryRegions", "HumanDHSFilter",

    function(obj, encode.table.name, chromosome, start, end, score.threshold=200) {

       driver <- RMySQL::MySQL()
       host <- "genome-mysql.cse.ucsc.edu"
       user <- "genome"
       dbname <- obj@genomeName

       db <- dbConnect(driver, user = user, host = host, dbname = dbname)

       #schema <- colnames(dbGetQuery(db, sprintf("select * from %s limit 1", tableName)))
       #suppressWarnings(dbGetQuery(db, sprintf("select * from %s limit 3", tableName)))

       #if(!all(c("chrom", "chromStart", "chromEnd") %in% schema)){
       #   printf("%s lacks chrom, start, end", tableName)
       #   lapply(dbListConnections(driver), dbDisconnect)
       #   return(data.frame())
       #   }

       main.clause <- sprintf("select * from %s where", encode.table.name);

      # Pull out the regions corresponding to the region in ENCODE
       query <- paste(main.clause,
                      sprintf("chrom = '%s'", chromosome),
                      sprintf("and chromStart >= %d", start),
                      sprintf("and chromEnd <= %d", end),
                      collapse = " ")

       printf("query: %s", query)

         # handle the usual case first: a start:end region many times larger than a typical DHS region
       suppressWarnings(  # MySQL returns unsigned integers.  hide these unproblematic conversion warnings
         tbl.regions <- dbGetQuery(db, query)
          )

       if(!obj@quiet)
           printf("%d DHS regions reported in %d bases, start:end unmodified", nrow(tbl.regions), 1 + end - start)

      # if no hits, then perhaps a very small region is requested, one which falls entirely within a DHS region
       if(nrow(tbl.regions) == 0) {
          extension <- 10000
          if(!obj@quiet)
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
         if(!obj@quiet) printf("query with extended region, %d rows", nrow(tbl.regionsExtended))
         if(nrow(tbl.regionsExtended) > 0) { # now find just the intersection of DHS and requested region
            if(!obj@quiet)
               printf("tbl.regionsExtended: %d rows, now have intersection", nrow(tbl.regionsExtended));
            gr.regions <- with(tbl.regionsExtended, GRanges(seqnames=chromosome, IRanges(chromStart, chromEnd)))
            gr.target <- GRanges(seqnames=chromosome, IRanges(start, end))
              # use range to collapse any multiple intersections down to that of the original target
            gr.intersect <- GenomicRanges::intersect(gr.target, gr.regions, ignore.strand=TRUE)
            if(!obj@quiet) printf("GenomicRanges intersections of extended region with original target: %d", length(gr.intersect))
            if(length(gr.intersect) >= 1){
               tbl.ov <- as.data.frame(findOverlaps(gr.intersect, gr.regions, type="any"))
               tbl.regions <- cbind(as.data.frame(gr.intersect[tbl.ov$queryHits]),
                                    tbl.regionsExtended[tbl.ov$subjectHits, c("name", "score")])
               colnames(tbl.regions) <- c("chrom", "chromStart", "chromEnd", "width", "strand", "name", "score")
               } # one or more region from the extended query intersects with the requested start:end.
            } # small region query, within DHS region
          } # small region, extension search

   lapply(dbListConnections(driver), dbDisconnect)

   tbl.regions$chrom <- as.character(tbl.regions$chrom)

       # the ucsc database call the  itemCount columns "name".  fix that
   if("name" %in% colnames(tbl.regions)){
      colnames(tbl.regions)[match("name", colnames(tbl.regions))] <- "count"
      }

   #tbl.regions$motif.start <- -1 + tbl.regions$regionStart + motif.start
   #tbl.regions$motif.end <- -1 + tbl.regions$regionEnd + motif.end
   invisible(tbl.regions[, c("chrom", "chromStart", "chromEnd",  "count", "score")])

   }) # getRegulatoryRegions

#----------------------------------------------------------------------------------------------------
setMethod("getSequence", "HumanDHSFilter",

   function(obj, tbl.regions){
     gr.regions <- with(tbl.regions, GRanges(seqnames=chrom, IRanges(start=chromStart, end=chromEnd)))
     seqs <- getSeq(obj@genome, gr.regions)
     as.character(seqs)
     })  # getSequence

#----------------------------------------------------------------------------------------------------
.matchForwardAndReverse <- function(sequence, pfm, motifName, min.match.percentage=95, quiet=TRUE)
{
   min.match.as.string <- sprintf("%02d%%", min.match.percentage)

   hits.fwd <- matchPWM(pfm, sequence, with.score=TRUE, min.score=min.match.as.string)
   hits.rev <- matchPWM(reverseComplement(pfm), sequence, with.score=TRUE, min.score=min.match.as.string)

   max.score <-  maxScore(pfm)
   tbl <- data.frame()
   if(length(hits.fwd) > 0){
      if(!quiet) printf("%d +", length(hits.fwd))
      match <- substring(as.character(subject(hits.fwd)), start(ranges(hits.fwd)), end(ranges(hits.fwd)))
      actual.score <- mcols(hits.fwd)$score
      relative.score <- actual.score/max.score
      tbl <- data.frame(ranges(hits.fwd),
                        score=mcols(hits.fwd)$score, maxScore=max.score,relativeScore=relative.score,
                        motif=motifName, match=match, strand="+",
                        stringsAsFactors=FALSE)
      } # hits.fwd

   if(length(hits.rev) > 0){
      if(!quiet) printf("%d -", length(hits.rev))
      match <- substring(as.character(subject(hits.rev)), start(ranges(hits.rev)), end(ranges(hits.rev)))
      actual.score <- mcols(hits.rev)$score
      relative.score <- actual.score/max.score
      tbl.rev <- data.frame(ranges(hits.rev),
                            score=mcols(hits.rev)$score, maxScore=max.score, relativeScore=relative.score,
                            motif=motifName, match=match, strand="-",
                            stringsAsFactors=FALSE)
         # transform the start/end so that they are forward-strand relative
      true.start <- 1 + nchar(sequence) - tbl.rev$end
      true.end   <- 1 + nchar(sequence) - tbl.rev$start
      tbl.rev$start <- true.start
      tbl.rev$end   <- true.end
      tbl <- rbind(tbl, tbl.rev)
      }

   #printf("returning match fwd/bwd tbl with %d rows", nrow(tbl))
   #print(table(tbl$strand))

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
   xx <- lapply(1:count, function(i) {
       .matchForwardAndReverse(sequence, pfms[[i]], names(pfms)[i], min.match.percentage, quiet)
       })

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
.parseLine <- function(textOfLine)
{
   # first delete the leading A, C, G or T.  then the square brackets.  then convert
   x <- substr(textOfLine, 2, nchar(textOfLine))
   x2 <- sub(" *\\[ *", "", x)
   x3 <- sub(" *\\] *", "", x2)
   counts <- as.integer(strsplit(x3, "\\s+", perl=TRUE)[[1]])

   return(counts)

} # .parseLine
#------------------------------------------------------------------------------------------------------------------------
.parseJasparPwm = function (lines)
{
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

} # .parsePwm
#------------------------------------------------------------------------------------------------------------------------
.readRawJasparMatrices = function (uri)
{
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

} # .readRawJasparMatrices
#------------------------------------------------------------------------------------------------------------------------
# TODO: move this duplicated code to the base class, or to R/utils.R
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

} # .parseChromLocString
#------------------------------------------------------------------------------------------------------------------------
# TODO: move this duplicated code to the base class, or to R/utils.R
.parseDatabaseUri <- function(database.uri)
{
   topLevel.tokens <- strsplit(database.uri, "://")[[1]]
   database.brand <- topLevel.tokens[1]
   secondLevel.tokens <- strsplit(topLevel.tokens[2], "/(?=[^/]+$)", perl = TRUE)[[1]]
   host <- secondLevel.tokens[1]
   database.name <- secondLevel.tokens[2]

   list(brand=database.brand, host=host, name=database.name)

} # .parseDatabaseUri
#----------------------------------------------------------------------------------------------------
geneSymbolToTSS
