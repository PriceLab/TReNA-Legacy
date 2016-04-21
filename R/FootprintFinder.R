printf <- function(...) print(noquote(sprintf(...)))
#------------------------------------------------------------------------------------------------------------------------
.FootprintFinder <- setClass("FootprintFinder",
                        slots = c(biomart="Mart",
                                  footprint.database.initializer="character",
                                  genes="character",
                                  tbl.locs="data.frame",
                                  state="environment",
                                  quiet="logical")
                        )
#------------------------------------------------------------------------------------------------------------------------
setGeneric("getGenes", signature="obj", function(obj) standardGeneric ("getGenes"))
setGeneric("getPromoterRegion", signature="obj", function(obj,  geneSymbol, size.upstream=1000, size.downstream=0)
                                                      standardGeneric("getPromoterRegion"))
setGeneric("getFootprints", signature="obj", function(obj,  geneSymbol, size.upstream=1000, size.downstream=0)
                                                      standardGeneric("getFootprints"))
#------------------------------------------------------------------------------------------------------------------------
FootprintFinder <- function(biomart, footprint.database.initializer, genes, quiet=TRUE)
{
   tokens <- strsplit(footprint.database.initializer, ":")[[1]]
   fp.db.protocol <- tokens[1]
   stopifnot(fp.db.protocol %in% c("sqlite"))
   fullpath <- tokens[2]
   stopifnot(file.exists(fullpath))

   tbl.locs <- .getGeneLocations (biomart, genes, quiet)

   state <- new.env(parent=emptyenv())
   state[["fp.db"]] == NULL

   if(fp.db.protocol == "sqlite"){
      if(!quiet) printf("FootprintFinder ctor opening sqlite connection to %s", fullpath)
      fp.db <- dbConnect(dbDriver("SQLite"), fullpath)
      state[["fp.db"]] <- fp.db
      } # sqlite

   if(is.null(state[["fp.db"]]))
       stop(sprintf("failed to connect to '%s'", footprint.database.initializer))

   .FootprintFinder(biomart=biomart,
               footprint.database.initializer=footprint.database.initializer,
               genes=genes,
               tbl.locs=tbl.locs,
               state=state,
               quiet=quiet)

} # FootprintFinder, the constructor
#------------------------------------------------------------------------------------------------------------------------
.getGeneLocations <- function(biomart, genes, quiet=TRUE)
{
    columns.desired <- c("entrezgene", "chromosome_name", "start_position", "end_position", "strand", "hgnc_symbol")
    if(!quiet)
       printf("biomart query for %d genes", length(genes))
    tbl.locs <- getBM(attributes=columns.desired, filters="hgnc_symbol", values=genes, mart=biomart)
    duplicates <- which(duplicated(tbl.locs$hgnc_symbol))
    if(length(duplicates) > 0)
       tbl.locs <- tbl.locs[-duplicates,]

    tbl.locs

} # .getGeneLocations
#------------------------------------------------------------------------------------------------------------------------
setMethod("getGenes", "FootprintFinder",

   function (obj){
      obj@genes
      })

#------------------------------------------------------------------------------------------------------------------------
setMethod("getPromoterRegion", "FootprintFinder",

   function(obj,  geneSymbol, size.upstream=1000, size.downstream=0){

      tbl.locs <- obj@tbl.locs

      stopifnot(geneSymbol %in% tbl.locs$hgnc_symbol)
      index <- grep(geneSymbol, tbl.locs$hgnc_symbol)[1]

      chrom <- sprintf("chr%s", tbl.locs$chromosome_name[index])
      start.orig <- tbl.locs$start_position[index]
      end.orig   <- tbl.locs$end_position[index]

      if(tbl.locs$strand[index] == -1){ # reverse (minus) strand.  TSS is at "end" position
         start.loc <- end.orig - size.downstream
         end.loc   <- end.orig + size.upstream
         }
      else{ #  forward (plus) strand.  TSS is at "start" position
        start.loc <- start.orig - size.upstream
        end.loc   <- start.orig + size.downstream
        }

     return(list(chr=chrom, start=start.loc, end=end.loc))
     })

#------------------------------------------------------------------------------------------------------------------------
setMethod("getFootprints", "FootprintFinder",

    function(obj,  geneSymbol, size.upstream=1000, size.downstream=0){
       stopifnot(length(geneSymbol) == 1)
       loc <- getPromoterRegion(obj, geneSymbol, size.upstream, size.downstream)
       query <- paste(c("select fp.chr, fp.mfpStart, fp.mfpEnd, fp.motifName, fp.pvalue, mg.motif, mg.tfs",
                        "from footprints fp",
                        "inner join motifsGenes mg",
                        "on fp.motifName=mg.motif",
                          sprintf("where fp.chr = '%s' and  fp.mfpStart > %d and fp.mfpEnd < %d",
                                  loc$chr, loc$start, loc$end)),
                        collapse=" ")
       dbGetQuery(obj@state[["fp.db"]], query)
       })

#----------------------------------------------------------------------------------------------------
