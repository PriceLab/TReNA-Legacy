#----------------------------------------------------------------------------------------------------
.HumanDHSFilter <- setClass("HumanDHSFilter",
                            contains="CandidateFilter",
                            slots=list(genomeName="character",
                                       genome="BSgenome",
                                       encodeTableName="character",
                                       pwmMatchPercentageThreshold="numeric",
                                       geneInfoDatabase.uri="character",   # access to gtf database
                                       geneCenteredSpec="list",
                                       regionsSpec="character",
                                       variants="character",
                                       pfms="list",
                                       quiet="logical"))

#----------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#----------------------------------------------------------------------------------------------------
setGeneric("getEncodeRegulatoryTableNames", signature="obj", function(obj) standardGeneric ("getEncodeRegulatoryTableNames"))
setGeneric("getRegulatoryRegions", signature="obj",
           function(obj, encode.table.name, chromosome, start, end, score.threshold=200, quiet=TRUE)
               standardGeneric ("getRegulatoryRegions"))
setGeneric("getSequence_tmp", signature="obj", function(obj, tbl.regions) standardGeneric ("getSequence_tmp"))
setGeneric("geneSymbolToTSS", signature="obj", function(obj, geneSymbol) standardGeneric("geneSymbolToTSS"))
#----------------------------------------------------------------------------------------------------
HumanDHSFilter <- function(genomeName,
                           encodeTableName="wgEncodeRegDnaseClustered",
                           #fimoDatabase.uri,
                           pwmMatchPercentageThreshold,
                           geneInfoDatabase.uri,
                           geneCenteredSpec=list(),
                           regionsSpec=list(),
                           variants=NA_character_,
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

   .HumanDHSFilter(CandidateFilter(quiet = quiet),
                   genomeName=genomeName,
                   #fimoDB=fimo.db,
                   pwmMatchPercentageThreshold=pwmMatchPercentageThreshold,
                   encodeTableName=encodeTableName,
                   geneInfoDatabase.uri=geneInfoDatabase.uri,
                   genome=reference.genome,
                   geneCenteredSpec=geneCenteredSpec,
                   regionsSpec=regionsSpec,
                   variants=variants,
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
        #print(dbListTables(db.geneInfo))
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

         # if present, extract a region from the geneCenteredSpec
       if(length(obj@geneCenteredSpec) == 3){
          expected.fields <- c("targetGene", "tssUpstream", "tssDownstream")
          stopifnot(all(expected.fields %in% names(obj@geneCenteredSpec)))
          x <- geneSymbolToTSS(obj)
          start <- x$tss - obj@geneCenteredSpec$tssUpstream
          end <- x$tss + obj@geneCenteredSpec$tssDownstream
          new.region.chromLocString <- sprintf("%s:%d-%d", x$chrom, start, end)
          regions <- c(regions, new.region.chromLocString)
          }

          # are there any explicit regions to include?
       if(!(all(is.na(obj@regionsSpec)))){
          regions <- c(regions, obj@regionsSpec)
          }

       printf("HumanDHSFilter::getCandidates");
       tbl.regions <- as.data.frame(do.call("rbind", lapply(regions, function(region) .parseChromLocString(region))))
       tbl.regions$chrom <- as.character(tbl.regions$chrom)
       tbl.regions$start <- as.integer(tbl.regions$start)
       tbl.regions$end <- as.integer(tbl.regions$end)
       tbl.dhs <- data.frame()
       for(r in 1:nrow(tbl.regions)){
          tbl.new <- getRegulatoryRegions(obj, obj@encodeTableName, tbl.regions$chrom[r], tbl.regions$start[r], tbl.regions$end[r])
          tbl.dhs <- rbind(tbl.dhs, tbl.new)
          }
       if(!obj@quiet)
          printf("found %d DHS regions in %d requested regions", nrow(tbl.dhs), nrow(tbl.regions))

       if(nrow(tbl.dhs) == 0){
          return(NA)
          }

       colnames(tbl.dhs) <- c("chrom", "start", "end", "count", "score")
       mm <- MotifMatcher(name="untitltedMotifMatcher", genomeName=obj@genomeName)
       x <- findMatchesByChromosomalRegion(mm, tbl.dhs,
                                           pwmMatchMinimumAsPercentage=obj@pwmMatchPercentageThreshold,
                                           variants=obj@variants)
       if(!obj@quiet)
          printf(" and %d motifs", nrow(x$tbl))

       colnames(x$tbl) <- c("motifName", "chrom", "motifStart", "motifEnd", "strand", "motifScore", "motifRelativeScore",
                            "match", "regulatoryRegionStart", "regualtoryRegionEnd", "regulatorySequence", "variant", "tfs")
       if(nchar(x$tbl$regulatorySequence[1]) > 20)
          x$tbl$regulatorySequence <- paste(substr(x$tbl$regulatorySequence, 1, 17), "...", sep="")
       x
    }) # getCandidates

#----------------------------------------------------------------------------------------------------
setMethod("getRegulatoryRegions", "HumanDHSFilter",

    function(obj, encode.table.name, chromosome, start, end, score.threshold=0) {

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

       if(!obj@quiet)
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
          if(!obj@quiet)
             printf("query with extension: %s", query)
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
setMethod("getSequence_tmp", "HumanDHSFilter",

   function(obj, tbl.regions){
     gr.regions <- with(tbl.regions, GRanges(seqnames=chrom, IRanges(start=chromStart, end=chromEnd)))
     seqs <- getSeq(obj@genome, gr.regions)
     as.character(seqs)
     })  # getSequence_tmp

#----------------------------------------------------------------------------------------------------
