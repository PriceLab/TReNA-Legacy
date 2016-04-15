printf <- function(...) print(noquote(sprintf(...)))
#------------------------------------------------------------------------------------------------------------------------
.Footprints <- setClass("Footprints",
                        slots = c(biomart="Mart",
                                  footprint.database.initializer="character",
                                  genes="character",
                                  state="environment")
                        )
#------------------------------------------------------------------------------------------------------------------------
setGeneric("getGenes", signature="obj", function(obj) standardGeneric ("getGenes"))
setGeneric("getPromoterRegion", signature="obj", function(obj,  geneSymbol, size.upstream=1000, size.downstream=0)
                                                      standardGeneric("getPromoterRegion"))
setGeneric("getFootprints", signature="obj", function(obj,  geneSymbol, size.upstream=1000, size.downstream=0)
                                                      standardGeneric("getFootprints"))
#------------------------------------------------------------------------------------------------------------------------
Footprints <- function(biomart, footprint.database.initializer, genes=genes)
{
   tokens <- strsplit(footprint.database.initializer, ":")[[1]]
   fp.db.protocol <- tokens[1]
   stopifnot(fp.db.protocol %in% c("sqlite"))
   fullpath <- tokens[2]
   stopifnot(file.exists(fullpath))

   state <- new.env(parent=emptyenv())
   state[["fp.db"]] == NULL

   if(fp.db.protocol == "sqlite"){
      printf("opening sqlite connection to %s", fullpath)
      fp.db <- dbConnect(dbDriver("SQLite"), fullpath)
      state[["fp.db"]] <- fp.db
      } # sqlite

   if(is.null(state[["fp.db"]]))
       stop(sprintf("failed to connect to '%s'", footprint.database.initializer))

   .Footprints(biomart=biomart,
               footprint.database.initializer=footprint.database.initializer,
               genes=genes,
               state=state)

} # Footprints, the constructor
#------------------------------------------------------------------------------------------------------------------------
setMethod("getGenes", "Footprints",

   function (obj){
      obj@genes
      })

#------------------------------------------------------------------------------------------------------------------------
setMethod("getPromoterRegion", "Footprints",

   function(obj,  geneSymbol, size.upstream=1000, size.downstream=0){

      columns.desired <- c("entrezgene", "chromosome_name", "start_position", "end_position", "strand", "hgnc_symbol")
      tbl.loc <- getBM(attributes=columns.desired, filters="hgnc_symbol", values=geneSymbol, mart=obj@biomart)[1,]
      chrom <- sprintf("chr%s", tbl.loc$chromosome_name)
      start.orig <- tbl.loc$start_position
      end.orig   <- tbl.loc$end_position
      if(tbl.loc$strand == -1){ # reverse (minus) strand.  TSS is at "end" position
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
setMethod("getFootprints", "Footprints",

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

# result <- getBM(attributes=c("entrezgene", "chromosome_name", "start_position", "end_position", "strand", "hgnc_symbol"),
#                   filters="hgnc_symbol", values=geneSymbol, mart=ensembl.hg38)

