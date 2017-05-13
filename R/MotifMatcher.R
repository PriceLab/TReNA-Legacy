.MotifMatcher <- setClass('MotifMatcher',
                          representation(name="character",
                                         genome="BSgenome",
                                         pfms="list",
                                         quiet="logical")
                        )
#------------------------------------------------------------------------------------------------------------------------
setGeneric("getPfms",            signature="obj", function(obj) standardGeneric ("getPfms"))
setGeneric("getSequence",        signature="obj", function(obj, tbl.regions, variants=NA_character_) standardGeneric ("getSequence"))
setGeneric("parseVariantString", signature="obj", function(obj, variantString) standardGeneric ("parseVariantString"))
setGeneric("findMatchesByChromosomalRegion", signature="obj",
           function(obj, tbl.regions, pwmMatchMinimumAsPercentage, variants=NA_character_)
              standardGeneric ("findMatchesByChromosomalRegion"))
#------------------------------------------------------------------------------------------------------------------------
MotifMatcher <- function(name=NA_character_,
                         genomeName="hg38",
                         pfms=list(),
                         quiet=TRUE)
{
   if(length(pfms) == 0){
      uri <- "http://jaspar.genereg.net/html/DOWNLOAD/JASPAR_CORE/pfm/nonredundant/pfm_vertebrates.txt"
      x <- .readRawJasparMatrices(uri)
         # normalize, so that a frequency sum of 1.0 is true across the 4 possible bases at each position
      pfms <- lapply(x, function(e) apply(e$matrix, 2, function(col) col/sum(col)))
      names(pfms) <- as.character(lapply(x, function(e) e$title))
      }

    if(genomeName == "hg38"){
       library(BSgenome.Hsapiens.UCSC.hg38)
       reference.genome <- BSgenome.Hsapiens.UCSC.hg38
       }
    else if(genomeName == "hg19"){
       library(BSgenome.Hsapiens.UCSC.hg19)
       reference.genome <- BSgenome.Hsapiens.UCSC.hg19
       }
    else {
      stop(sprintf("MotifMatch, genomeName not in hg19, hg38: '%s'", genomeName))
       }


   .MotifMatcher(name=name, genome=reference.genome,
                 pfms=pfms, quiet=quiet)

} # MotifMatcher constructor
#------------------------------------------------------------------------------------------------------------------------
setMethod("show", "MotifMatcher",

     function(object){
        s <- sprintf("MotifMatcher...")
        cat(s, sep="\n")
        })
#----------------------------------------------------------------------------------------------------
setMethod("findMatchesByChromosomalRegion", "MotifMatcher",

    function(obj, tbl.regions, pwmMatchMinimumAsPercentage, variants=NA_character_){
       #stopifnot(nrow(tbl.regions) == 1)   # simplify our manipulations.
       browser()
       tbl.regions <- getSequence(obj, tbl.regions, variants)
       colnames(tbl.regions) <- c("chrom", "chromStart", "chromEnd", "seq", "variant")   # preprend 'chrom' to start and end, to distinguish
       tbl.motifs <- .getScoredMotifs(tbl.regions$seq, pwmMatchMinimumAsPercentage, obj@quiet)
           # tbl.mg will soon come from MotifDb
       tbl.mg <- read.table(system.file(package="TReNA", "extdata", "motifGenes.tsv"), sep="\t", as.is=TRUE, header=TRUE)
       region.count <- nrow(tbl.regions)
       tbl.out <- data.frame()
       all.tfs <- c()
       for(i in seq_len(region.count)){
         if(nrow(tbl.motifs[[i]]) > 0)
            tbl.out <- rbind(tbl.out, cbind(tbl.motifs[[i]], tbl.regions[i,]))
            } # for i
       if(nrow(tbl.out) == 0)
          return(list(tbl=data.frame(), tfs=c()))

       tbl.out$start <- tbl.out$start + tbl.out$chromStart - 1;
       tbl.out$end <- tbl.out$end  + tbl.out$chromStart - 1;
             # change some column names
       colnames(tbl.out)[4] <- "motifscore"
       colnames(tbl.out)[2] <- "endpos"
       colnames(tbl.out)[7] <- "motifname"
       tbl.out <- tbl.out[, -c(3,5)] # get rid of width and maxScore columns
       desired.column.order <- c("motifname", "chrom", "start", "endpos", "strand", "motifscore", "relativeScore", "match",
                                 "chromStart", "chromEnd", "seq", "variant") #, "count", "score")
       tbl.out <- tbl.out[, desired.column.order]
       tbl.out <- tbl.out[order(tbl.out$motifscore, decreasing=TRUE),]
       tfs.by.motif <- lapply(tbl.out$motifname, function(m) subset(tbl.mg, motif==m)$tf.gene)
       all.tfs <- sort(unique(unlist(tfs.by.motif)))
       tfs.by.motif.joined <- unlist(lapply(tfs.by.motif, function(m) paste(m, collapse=";")))
       tbl.out$tf <- tfs.by.motif.joined
             # tentative, empirically obtained cutoff
          #upper.quartile.threshold <- fivenum(tbl.out$motifscore)[4]
          #tbl.out <- subset(tbl.out, motifscore >= upper.quartile.threshold)
       list(tbl=tbl.out, tfs=all.tfs)
       })

#----------------------------------------------------------------------------------------------------
#setMethod("findMatchesBySequence", "MotifMatcher",
#
#    function(obj, sequences){
#           # tbl.mg will soon come from MotifDb
#       if(length(sequences) == 0){
#          warning("zero sequences supplied, returning empty data.frame")
#          return(data.frame())
#          }
#       browser()
#       tbl.mg <- read.table(system.file(package="TReNA", "extdata", "motifGenes.tsv"), sep="\t", as.is=TRUE, header=TRUE)
#       tbl.motifs <- .getScoredMotifs(sequences, obj@pwmMatchMinimumAsPercentage, obj@quiet)
#       region.count <- length(sequences)
#       tbl.out <- data.frame()
#       all.tfs <- c()
#       for(i in seq_len(region.count)){
#         if(nrow(tbl.motifs[[i]]) > 0)
#            tbl.out <- rbind(tbl.out, cbind(tbl.motifs[[i]], tbl.regions[i,]))
#            } # for i
#      tbl.out$start <- tbl.out$start + tbl.out$chromStart - 1;
#      tbl.out$end <- tbl.out$end  + tbl.out$chromStart - 1;
#             # change some column names
#      colnames(tbl.out)[4] <- "motifscore"
#      colnames(tbl.out)[2] <- "endpos"
#      colnames(tbl.out)[7] <- "motifname"
#      tbl.out <- tbl.out[, -c(3,5)] # get rid of width and maxScore columns
#      desired.column.order <- c("motifname", "chrom", "start", "endpos", "strand", "motifscore", "relativeScore", "match",
#                                "chromStart", "chromEnd", "count", "score")
#      tbl.out <- tbl.out[, desired.column.order]
#      tbl.out <- tbl.out[order(tbl.out$motifscore, decreasing=TRUE),]
#      tfs.by.motif <- lapply(tbl.out$motifname, function(m) subset(tbl.mg, motif==m)$tf.gene)
#      all.bioc.tfs <- sort(unique(unlist(tfs.by.motif)))
#      tfs.by.motif.joined <- unlist(lapply(tfs.by.motif, function(m) paste(m, collapse=";")))
#      tbl.out$tf <- tfs.by.motif.joined
#             # tentative, empirically obtained cutoff
#          #upper.quartile.threshold <- fivenum(tbl.out$motifscore)[4]
#          #tbl.out <- subset(tbl.out, motifscore >= upper.quartile.threshold)
#       tbl.out
#       })
#

#------------------------------------------------------------------------------------------------------------------------
setMethod("getPfms", "MotifMatcher",

          function(obj){
             return(obj@pfms)
          })
#------------------------------------------------------------------------------------------------------------------------
.matchPwmForwardAndReverse <- function(sequence, pfm, motifName, min.match.percentage=95, quiet=TRUE)
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

} # .matchPwmForwardAndReverse
#----------------------------------------------------------------------------------------------------
# .findMotifs <- function(sequence, pfms, min.match.percentage=95, quiet=TRUE)
# {
#    min.match.as.string <- sprintf("%02d%%", min.match.percentage)
#
#    search <- function(motifName, mtx, seq){
#       hits.fwd <- matchPWM(mtx, seq, with.score=TRUE, min.score=min.match.as.string)
#       seq.revcomp <- as.character(reverseComplement(DNAString(seq)))
#       hits.rev <- matchPWM(mtx, seq.revcomp, with.score=TRUE, min.score=min.match.as.string)
#       #browser()
#       tbl <- data.frame()
#       if(length(hits.fwd) > 0){
#           if(!quiet) printf("%d +", length(hits.fwd))
#           match <- substring(as.character(subject(hits.fwd)), start(ranges(hits.fwd)), end(ranges(hits.fwd)))
#           tbl <- data.frame(ranges(hits.fwd), score=mcols(hits.fwd)$score, motif=motifName, match=match, strand="+")
#           }
#       if(length(hits.rev) > 0){
#           if(!quiet) printf("%d -", length(hits.rev))
#           match <- substring(as.character(subject(hits.rev)), start(ranges(hits.rev)), end(ranges(hits.rev)))
#           tbl.rev <- data.frame(ranges(hits.rev), score=mcols(hits.rev)$score, motif=motifName, match=match, strand="-")
#           tbl <- rbind(tbl, tbl.rev)
#           }
#        return(tbl)
#        }
#
#     count <- length(pfms)
#     xx <- lapply(1:count, function(i) {
#        .matchForwardAndReverse(sequence, pfms[[i]], names(pfms)[i], min.match.percentage, quiet)
#        })
#
#     tbl.result <- do.call("rbind", xx)
#     if(nrow(tbl.result) == 0){
#       return(data.frame())
#       }
#     else{
#       tbl.result$motif <- as.character(tbl.result$motif)
#       #tbl.result$seq <- as.character(tbl.result$seq)
#       tbl.result$match <- as.character(tbl.result$match)
#       tbl.result$strand <- as.character(tbl.result$strand)
#       return(tbl.result[order(tbl.result$score,decreasing=TRUE),])
#       }
#
# }  # .findMotifs
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
   if(is.null(result))
      result <- data.frame()

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
setMethod("getSequence", "MotifMatcher",

   function(obj, tbl.regions, variants=NA_character_){
        # either no variants are supplied, or a variant string is offered for each row
     #.distributeVariants(tbl.regions, variants)
     stopifnot(is.na(variants) || nrow(tbl.regions) == length(variants))
     gr.regions <- with(tbl.regions, GRanges(seqnames=chrom, IRanges(start=start, end=end)))
     seqs <- as.character(getSeq(obj@genome, gr.regions))
     tbl.regions$seq <- seqs
     variant.column <- rep("", nrow(tbl.regions))
     if(!is.na(variants)){
        i <- 0
        for(variant in variants){
           i <- i + 1
           variant.info <- parseVariantString(obj, variant)
           browser()
           variant.column[i] <- with(variant.info, sprintf("%s:%d:%s", chrom, loc, mut))
           gr.variants <- with(variant.info, GRanges(seqnames=chrom, IRanges(start=loc, end=loc)))
           tbl.overlaps <- as.data.frame(findOverlaps(gr.variants, gr.regions))
           for(r in 1:nrow(tbl.overlaps)){
              tbl.regions <- .injectSnp(tbl.regions, tbl.overlaps$queryHits[r], tbl.overlaps$subjectHits[r],
                                        variant.info)
           } # for r
          } # for variant
       } # !is.na(variants)
      tbl.regions$variant <- variant.column
      tbl.regions
     })  # getSequence

#----------------------------------------------------------------------------------------------------
.injectSnp <- function(tbl.regions, queryIndex, subjectIndex, variant.info)
{
   wt <- as.list(tbl.regions[subjectIndex,])
   stopifnot(variant.info$loc <= wt$end)
   stopifnot(variant.info$loc >= wt$start)
   stopifnot(wt$chrom == variant.info$chrom)

   offset <- 1 + variant.info$loc - wt$start
   mut.loc <- wt$start + offset
   wt.base <- substring(wt$seq, offset, offset)
   mutant.seq <- sprintf("%s%s%s", substr(wt$seq, 1, (offset-1)), variant.info$mut,
                                   substr(wt$seq, offset+1, nchar(wt$seq)))
      # create as many (initially identical) rows as there are mutant.seqs
   mutant.seq.count <- length(mutant.seq)
   tbl.regions.new <- as.data.frame(lapply(tbl.regions, rep, mutant.seq.count))
   tbl.regions.new$seq <- mutant.seq
   tbl.regions.new

} # .injextSnp
#----------------------------------------------------------------------------------------------------
setMethod("parseVariantString", "MotifMatcher",

    function(obj, variantString) {
       if(grepl("^rs", variantString)){
          require(SNPlocs.Hsapiens.dbSNP144.GRCh38)
          snp.info <- as.data.frame(snpsById(SNPlocs.Hsapiens.dbSNP144.GRCh38, variantString))[1,]
          chrom <- as.character(snp.info$seqnames)
          if(!grepl("chr", chrom))
             chrom <- sub("ch", "chr", chrom)
          start <- as.numeric(snp.info$pos)
          ambiguity.code <- snp.info$alleles_as_ambig
          elements.string <- IUPAC_CODE_MAP[[ambiguity.code]]
          elements <- strsplit(elements.string,'')[[1]]
          wt <- as.character(getSeq(obj@genome, chrom, start, start))
          mut <- setdiff(elements, wt)
          snp <- list(chrom=chrom, loc=start, wt=wt, mut=mut)
          variant.count <- length(snp$mut)
          result <- vector(mode="list", length=variant.count)
          for(v in 1:variant.count){
             result[[v]] <- snp
             result[[v]]$mut <- snp$mut[v]
             } # for v
          tbl.out <- do.call(rbind.data.frame, result)
          } # if rsid
       else{
          tokens <- strsplit(variantString, ":")[[1]]
          stopifnot(length(tokens) == 4)
          tbl.out <- data.frame(chrom=tokens[1], loc=as.numeric(tokens[2]), wt=tokens[3], mut=tokens[4])
          }
       tbl.out$chrom <- as.character(tbl.out$chrom)
       tbl.out$wt <- as.character(tbl.out$wt)
       tbl.out$mut <- as.character(tbl.out$mut)
       tbl.out
    })

#----------------------------------------------------------------------------------------------------
