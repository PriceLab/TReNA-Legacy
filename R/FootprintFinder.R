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
                             slots = c(genome.db="DBIConnection",
                                       project.db="DBIConnection",
                                       quiet="logical")
                            )
#------------------------------------------------------------------------------------------------------------------------
setGeneric("getChromLoc", signature="obj", function(obj, name, biotype="protein_coding",
                                                    moleculetype="gene") standardGeneric("getChromLoc"))
setGeneric("getGenePromoterRegion", signature="obj", function(obj,  gene, size.upstream=1000, size.downstream=0,
                                                              biotype="protein_coding", moleculetype="gene")
                                                      standardGeneric("getGenePromoterRegion"))
setGeneric("getFootprintsForGene", signature="obj", function(obj,  gene, size.upstream=1000, size.downstream=0,
                                                                   biotype="protein_coding", moleculetype="gene")
                                                      standardGeneric("getFootprintsForGene"))
setGeneric("getFootprintsInRegion", signature="obj", function(obj, chromosome, start, end) standardGeneric("getFootprintsInRegion"))
setGeneric("getGtfGeneBioTypes", signature="obj", function(obj) standardGeneric("getGtfGeneBioTypes"))
setGeneric("getGtfMoleculeTypes", signature="obj", function(obj) standardGeneric("getGtfMoleculeTypes"))
setGeneric("closeDatabaseConnections", signature="obj", function(obj) standardGeneric("closeDatabaseConnections"))
setGeneric("getPromoterRegionsAllGenes",signature="obj", function(obj ,size.upstream=10000 , size.downstream=10000 ) standardGeneric("getPromoterRegionsAllGenes"))
#------------------------------------------------------------------------------------------------------------------------
.parseDatabaseUri <- function(database.uri)
{
   topLevel.tokens <- strsplit(database.uri, "://")[[1]]
   database.brand <- topLevel.tokens[1]
   secondLevel.tokens <- strsplit(topLevel.tokens[2], "/")[[1]]
   host <- secondLevel.tokens[1]
   database.name <- secondLevel.tokens[2]

   list(brand=database.brand, host=host, name=database.name)

} # .parseDatabaseUri
#------------------------------------------------------------------------------------------------------------------------
FootprintFinder <- function(genome.database.uri, project.database.uri, quiet=TRUE)
{
   genome.db.info <- .parseDatabaseUri(genome.database.uri)
   project.db.info <- .parseDatabaseUri(project.database.uri)
   stopifnot(genome.db.info$brand %in% c("postgres"))


      # open the genome database
   if(genome.db.info$brand == "postgres"){
      host <- genome.db.info$host
      dbname <- genome.db.info$name
      driver <- PostgreSQL()
      genome.db <- dbConnect(driver, user= "trena", password="trena", dbname=dbname, host=host)
      existing.databases <- dbGetQuery(genome.db, "select datname from pg_database")[,1]
      stopifnot(dbname %in% existing.databases)
      dbDisconnect(genome.db)
      genome.db <- dbConnect(driver, user="trena", password="trena", dbname=dbname, host=host)
      expected.tables <- c("gtf", "motifsgenes")
      stopifnot(all(expected.tables %in% dbListTables(genome.db)))
      if(!quiet){
         row.count <- dbGetQuery(genome.db, "select count(*) from gtf")[1,1]
         printf("%s: %d rows", sprintf("%s/gtf", genome.database.uri), row.count)
         row.count <- dbGetQuery(genome.db, "select count(*) from motifsgenes")[1,1]
         printf("%s: %d rows", sprintf("%s/motifsgenes", genome.database.uri), row.count)

         }
      } # if postgres

      # open the project database
   if(project.db.info$brand == "postgres"){
      host <- project.db.info$host
      dbname <- project.db.info$name
      driver <- PostgreSQL()
      project.db <- dbConnect(driver, user= "trena", password="trena", dbname=dbname, host=host)
      existing.databases <- dbGetQuery(project.db, "select datname from pg_database")[,1]
      stopifnot(dbname %in% existing.databases)
      dbDisconnect(project.db)
      project.db <- dbConnect(driver, user="trena", password="trena", dbname=dbname, host=host)
      expected.tables <- c("footprints")
      stopifnot(all(expected.tables %in% dbListTables(project.db)))
      if(!quiet){
         row.count <- dbGetQuery(project.db, "select count(*) from footprints")[1,1]
         printf("%s: %d rows", sprintf("%s/footprints", project.database.uri), row.count)
         }
     } # if postgres


   .FootprintFinder(genome.db=genome.db, project.db=project.db, quiet=quiet)

} # FootprintFinder, the constructor
#------------------------------------------------------------------------------------------------------------------------
setMethod("closeDatabaseConnections", "FootprintFinder",

     function(obj){
        if("DBIConnection" %in% is(obj@genome.db))
          dbDisconnect(obj@genome.db)
        if("DBIConnection" %in% is(obj@project.db))
          dbDisconnect(obj@project.db)
          })

#------------------------------------------------------------------------------------------------------------------------
setMethod("getGtfGeneBioTypes", "FootprintFinder",

     function(obj){
        sort(dbGetQuery(obj@genome.db, "select distinct gene_biotype from gtf")[,1])
        })

#------------------------------------------------------------------------------------------------------------------------
setMethod("getGtfMoleculeTypes", "FootprintFinder",

     function(obj){
        sort(dbGetQuery(obj@genome.db, "select distinct moleculetype from gtf")[,1])
        })

#------------------------------------------------------------------------------------------------------------------------
setMethod("getChromLoc", "FootprintFinder",

   function(obj, name, biotype="protein_coding", moleculetype="gene"){
      query <- paste("select gene_id, gene_name, chr, start, endpos, strand from gtf where ",
                     sprintf("(gene_name='%s' or gene_id='%s') ", name, name),
                     sprintf("and gene_biotype='%s' and moleculetype='%s'", biotype, moleculetype),
                     collapse=" ")
     dbGetQuery(obj@genome.db, query)
     })

#------------------------------------------------------------------------------------------------------------------------
setMethod("getGenePromoterRegion", "FootprintFinder",

   function(obj, gene, size.upstream=1000, size.downstream=0, biotype="protein_coding", moleculetype="gene"){

      tbl.loc <- getChromLoc(obj, gene, biotype=biotype, moleculetype=moleculetype)
      if(nrow(tbl.loc) != 1){
          warning(sprintf("no chromosomal location for %s (%s, %s)", gene, biotype, moleculetype))
          return(NA)
          }

      chrom <- tbl.loc$chr[1]
      start.orig <- tbl.loc$start[1]
      end.orig   <- tbl.loc$endpos[1]
      strand     <- tbl.loc$strand[1]

      if(strand == "-"){ # reverse (minus) strand.  TSS is at "end" position
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
setMethod("getFootprintsForGene", "FootprintFinder",

    function(obj,  gene, size.upstream=1000, size.downstream=0, biotype="protein_coding", moleculetype="gene"){
       stopifnot(length(gene) == 1)
       loc <- getGenePromoterRegion(obj, gene, size.upstream, size.downstream,
                                    biotype=biotype, moleculetype=moleculetype)
       if(!obj@quiet) print(loc)
       getFootprintsInRegion(obj, loc$chr, loc$start, loc$end)
       }) # getFootprintsForGene

#----------------------------------------------------------------------------------------------------
setMethod("getFootprintsInRegion", "FootprintFinder",

    function(obj, chromosome, start, end){
       query <- paste(c("select fp.chr, fp.mfpstart, fp.mfpend, fp.motifname, fp.pval, mg.motif, mg.tf_name, mg.tf_ensg",
                        "from footprints fp",
                        "inner join motifsgenes mg",
                        "on fp.motifName=mg.motif",
                          sprintf("where fp.chr = '%s' and  fp.mfpstart >= %d and fp.mfpend <= %d",
                                  chromosome, start, end)),
                        collapse=" ")
       if(!obj@quiet) print(query)
       dbGetQuery(obj@project.db, query)
       }) # getFootprintsInRegion

#----------------------------------------------------------------------------------------------------
setMethod("getPromoterRegionsAllGenes","FootprintFinder",

   function( obj , size.upstream=10000 , size.downstream=10000 ) {

   query <-
   paste( "select gene_name, gene_id, chr, start, endpos, strand from gtf where" ,
        "gene_biotype='protein_coding' and moleculetype='gene'" , sep=" " )

   genes = dbGetQuery( obj@genome.db , query )

   # function to get each transcript's TSS
   get_tss <-
   function( t ) {
      chrom <- genes$chr[t]
      start.orig <- genes$start[t]
      end.orig   <- genes$endpos[t]
      strand     <- genes$strand[t]

      if(strand == "-"){ # reverse (minus) strand.  TSS is at "end" position
         tss <- end.orig
         }
      else{ #  forward (plus) strand.  TSS is at "start" position
        tss <- start.orig
        }
     return( tss )
   }
   # apply get_tss to all transcripts
   tss = sapply( 1:nrow(genes) , get_tss )

   # assemble a bed file for the TSSs
   promoter_regions = unique( data.frame(
        chr = genes$chr ,
        start = tss - size.upstream , end = tss + size.downstream ,
        gene_name = genes$gene_name ,
        gene_id = genes$gene_id ))
   # GRanges obj
   gr = makeGRangesFromDataFrame( promoter_regions , keep.extra.columns = T )
   names(gr) = promoter_regions$gene_name
   return( gr )

}) # getPromoterRegionsAllGenes

















