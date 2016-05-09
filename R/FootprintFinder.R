printf <- function(...) print(noquote(sprintf(...)))
#------------------------------------------------------------------------------------------------------------------------
# a FootprintFinder needs, at present, this sort of information
#
#   1) a complete table of genome features, for which our current model is ensembl's Homo_sapiens.GRCh38.84.chr.gtf
#   2) a footprint table, the output from cory's pipeline
#   3) a motif/TF map, with scores & etc
# 
#------------------------------------------------------------------------------------------------------------------------
.FootprintFinder <- setClass("FootprintFinder",
                        slots = c(database.uri="character",
                                  state="environment",
                                  quiet="logical")
                        )
#------------------------------------------------------------------------------------------------------------------------
setGeneric("getPromoterRegion", signature="obj", function(obj,  geneSymbol, size.upstream=1000, size.downstream=0)
                                                      standardGeneric("getPromoterRegion"))
setGeneric("getFootprints", signature="obj", function(obj,  geneSymbol, size.upstream=1000, size.downstream=0)
                                                      standardGeneric("getFootprints"))
#------------------------------------------------------------------------------------------------------------------------
FootprintFinder <- function(database.uri, quiet=TRUE)
{
   tokens <- strsplit(database.initializer, ":")[[1]]
   stopifnot(length(tokens) == 2)
   db.protocol <- tokens[1]
   db.description <- tokens[2]
   stopifnot(db.protocol %in% c("sqlite", "postgres"))
   tbl.locs <- .getGeneLocations (, genes, quiet)

   state <- new.env(parent=emptyenv())
   state[["db"]] == NULL
 
  if(!quiet) printf("FootprintFinder ctor opening database connection to %s", fullpath)

   if(db.protocol == "sqlite"){
      fullpath <- db.description
      stopifnot(file.exists(fullpath))
      db <- dbConnect(dbDriver("SQLite"), fullpath)
      state[["db"]] <- db
      } # sqlite

   if(db.protocol == "postgres"){ 
      postgres.tokens <- strsplit(db.description, "/")[[1]]
      stopifnot(length(postgres.tokens) == 3)
      db.host <- tokens[1]
      db.databaseName <- tokens[2]
      genome.assembly <- tokens[3]
      driver <- PostgreSQL()
      db.fp <- dbConnect(driver, user= "trena", password="trena", host=db.host)
      existing.databases <- dbGetQuery(db, "select datname from pg_database")[,1]      
      stopifnot(db.databaseName %in% existing.databases)
      dbDisconnect(db.fp)
      db.fp <- dbConnect(driver, user="trena", password="trena", dbname=db.databaseName, host=db.host)
      stopifnot("gtf" %in% existing.databases)
      db.gtf <- dbConnect(driver, user="trena", password="trena", dbname="gtf", host=db.host)
      stage[["db.fp"]] <- db.fp
      stage[["db.gtf"]] <- db.gtf
      }

   if(is.null(state[["db.fp"]]))
     stop(sprintf("failed to connect to footprint database specified in '%s'", db.description)

   if(is.null(state[["db.gtf"]]))
     stop(sprintf("failed to connect to gtf (genome) information specified in '%s'", db.description)


   .FootprintFinder(database.initializer=database.initializer,
                    state=state,
                    quiet=quiet)

} # FootprintFinder, the constructor
#------------------------------------------------------------------------------------------------------------------------
.getGeneLocations <- function(biomart, genes, quiet=TRUE)
{
    columns.desired <- c("entrezgene", "chromosome_name", "start_position", "end_position", "strand", "hgnc_symbol",
                         "ensembl_gene_id")
    if(!quiet)
       printf("biomart query for %d genes", length(genes))

        # identify ENSG ids.  lookup locations for them
    tbl.locs.ensg <- data.frame()
    tbl.locs.geneSymbol <- data.frame()

    ensg.ids <- grep("^ENSG", genes)

    if(length(ensg.ids) > 0){
       tbl.locs.ensg <- getBM(attributes=columns.desired, filters="ensembl_gene_id", values=genes, mart=biomart)
       }

    non.ensg.ids <- setdiff(genes, ensg.ids)
    if(length(non.ensg.ids) > 0){
       tbl.locs.geneSymbol <- getBM(attributes=columns.desired, filters="hgnc_symbol", values=genes, mart=biomart)
       duplicates <- which(duplicated(tbl.locs.geneSymbol$hgnc_symbol))
       if(length(duplicates) > 0)
          tbl.locs.geneSymbol <- tbl.locs.geneSymbol[-duplicates,]
       }

    rbind(tbl.locs.ensg, tbl.locs.geneSymbol)

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
      if(geneSymbol %in% tbl.locs$hgnc_symbol){
          index <- grep(geneSymbol, tbl.locs$hgnc_symbol)[1]
      } else if (geneSymbol %in% tbl.locs$ensembl_gene_id){
          index <- grep(geneSymbol, tbl.locs$ensembl_gene_id)
      } else {
          stop(sprintf("FootprintFinder:getPromoterRegion does not recognize geneSymbol '%s'",
                       geneSymbol))
      }

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
       query <- paste(c("select fp.chr, fp.mfpStart, fp.mfpEnd, fp.motifName, fp.pvalue, mg.motif, mg.tf",
                        "from footprints fp",
                        "inner join motifsGenes mg",
                        "on fp.motifName=mg.motif",
                          sprintf("where fp.chr = '%s' and  fp.mfpStart > %d and fp.mfpEnd < %d",
                                  loc$chr, loc$start, loc$end)),
                        collapse=" ")
       dbGetQuery(obj@state[["db"]], query)
       })

#----------------------------------------------------------------------------------------------------
