.MotifMatcher <- setClass('MotifMatcher',
                       representation(name="character",
                                      sequence="character",
                                      pwmMatchMinimumAsPercentage="numeric",
                                      pfms="list")
                        )
#------------------------------------------------------------------------------------------------------------------------
setGeneric("findMatches",      signature="obj", function(obj) standardGeneric ("findMatches"))
setGeneric("getPfms",          signature="obj", function(obj) standardGeneric ("getPfms"))
#------------------------------------------------------------------------------------------------------------------------
MotifMatcher <- function(name=NA_character_,
                   sequence=NA_character_,
                   pwmMatchMinimumAsPercentage=90,
                   pfms=list())
{
   if(length(pfms) == 0){
      uri <- "http://jaspar.genereg.net/html/DOWNLOAD/JASPAR_CORE/pfm/nonredundant/pfm_vertebrates.txt"
      x <- .readRawJasparMatrices(uri)
         # normalize, so that a frequency sum of 1.0 is true across the 4 possible bases at each position
      pfms <- lapply(x, function(e) apply(e$matrix, 2, function(col) col/sum(col)))
      names(pfms) <- as.character(lapply(x, function(e) e$title))
      }

   .MotifMatcher(name=name, sequence=sequence, pwmMatchMinimumAsPercentage=pwmMatchMinimumAsPercentage, pfms=pfms)

} # MotifMatcher constructor
#------------------------------------------------------------------------------------------------------------------------
setMethod("findMatches", "MotifMatcher",

          function(obj){
             return(data.frame())
          })
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





