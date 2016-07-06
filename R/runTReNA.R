#----------------------------------------------------------------------------------------------------
getFootprintsForTF <-
function( obj , tf ){
       query <-
paste( "select fp.chr, fp.mfpstart, fp.mfpend, fp.motifname, fp.pval, fp.score, mg.motif, mg.tf",
             "from footprints fp",
             "inner join motifsgenes mg",
             "on fp.motifName=mg.motif",
             paste("where mg.tf='",tf,"'" , sep="" ) ,  collapse = " " )

       if(!obj@quiet) print(query)
       dbGetQuery( obj@project.db , query )
} # getFootprintsForTF
#----------------------------------------------------------------------------------------------------
getGenePromoterRegions <-
function( obj , genelist , size.upstream = 10000 , size.downstream = 10000 ) {

  # note: this function runs very slowly. if the goal is to get the promoter regions for all genes,
  # use getPromoterRegionsAllGenes()

   promoter_regions = 
   lapply( 1:length(genelist) , function(x) {
      getGenePromoterRegion( 
         obj , 
         genelist[x] , 
         size.upstream = size.upstream , 
         size.downstream = size.downstream 
   )})
   chr = sapply( 1:length(promoter_regions) , function(x) promoter_regions[[x]]$chr )
   start = sapply( 1:length(promoter_regions) , function(x) promoter_regions[[x]]$start )
   end = sapply( 1:length(promoter_regions) , function(x) promoter_regions[[x]]$end )
   pr = data.frame( chr , start , end )
   promoter_regions.gr = makeGRangesFromDataFrame( pr )
   names(promoter_regions.gr) = genelist

   return( promoter_regions.gr )

} #getGenePromoterRegions
#----------------------------------------------------------------------------------------------------
getTfbsCountsPerGene <-
function( obj , tflist , promoter_regions , verbose = 1 ) {

   tfbs_counts_per_gene = list()
   for( x in 1:length(tflist) ) {
      if( verbose > 1 ) cat( "...Working on" , tflist[x] , "(" , x , "/" , length(x) ,  ")\n" )
      footprints.tf = getFootprintsForTF( obj , tflist[x] )
      if( nrow( footprints.tf ) == 0 ) {
         counts = rep( 0 , length(promoter_regions) )
         tfbs_counts_per_gene[[x]] = counts
      }
      colnames( footprints.tf )[1:3] = c("chr","start","end")
      footprints.tf.gr = makeGRangesFromDataFrame( footprints.tf )
      fp.reduced = reduce( footprints.tf.gr )
      counts = countOverlaps( promoter_regions , fp.reduced )
      tfbs_counts_per_gene[[x]] = counts
      if( verbose == 1 & x/10 == round(x/10) ) cat("...done" , x , "/" , length(tflist) ,"\n" )
   }
   tfbs_counts_per_gene = do.call( cbind , tfbs_counts_per_gene )
   colnames(tfbs_counts_per_gene) = tflist
   rownames(tfbs_counts_per_gene) = names(promoter_regions)
   return( tfbs_counts_per_gene )

} #getTfbsCountsPerGene
#----------------------------------------------------------------------------------------------------
getTfbsCountsPerGeneMC <-
function( tflist , promoter_regions , verbose = 1 , cores = 3 , genome="hg38" , tissue="lymphoblast") {

   library(doParallel)
   registerDoParallel( cores = cores )

   tfbs_counts_per_gene =
   foreach( x=1:length(tflist) ) %dopar% {
      genome.db.uri <- paste("postgres://whovian/" , genome , sep="" )
      project.db.uri <-  paste("postgres://whovian/" , tissue , sep="" )
      obj <- FootprintFinder(genome.db.uri, project.db.uri, quiet=TRUE)
      if( verbose > 1 ) cat( "...Working on" , tflist[x] , "(" , x , "/" , length(x) ,  ")\n" )
      footprints.tf = getFootprintsForTF( obj , tflist[x] )
      closeDatabaseConnections( obj )

      if( nrow( footprints.tf ) == 0 ) {
         counts = rep( 0 , length(promoter_regions) )
         return(counts)
      }
      colnames( footprints.tf )[1:3] = c("chr","start","end")
      footprints.tf.gr = makeGRangesFromDataFrame( footprints.tf )
      fp.reduced = reduce( footprints.tf.gr )
      counts = countOverlaps( promoter_regions , fp.reduced )
      if( verbose == 1 & x/10 == round(x/10) ) cat("...done" , x , "/" , length(tflist) ,"\n" )
      return( counts )
   }
   tfbs_counts_per_gene = do.call( cbind , tfbs_counts_per_gene )
   colnames(tfbs_counts_per_gene) = tflist
   rownames(tfbs_counts_per_gene) = names(promoter_regions)
   return( tfbs_counts_per_gene )

} #getTfbsCountsPerGene
#----------------------------------------------------------------------------------------------------
makeTfbsCountsTbl <- 
function( genome = "hg38" , tissue = "lymphoblast" ,  genelist = NULL , tflist = NULL , 
   biotype = "protein_coding" , moleculetype = "gene" ,
   size.upstream = 10000 , size.downstream = 10000 , cores = 1 , verbose = 1 ) {

   genome.db.uri <- paste("postgres://whovian/" , genome , sep="" )
   project.db.uri <-  paste("postgres://whovian/" , tissue , sep="" )
   obj <- FootprintFinder(genome.db.uri, project.db.uri, quiet=TRUE)

   if( is.null( genelist ) ) {
      if( verbose >= 1 ) 
          cat("no gene list is given. using all gene_ids from obj@genome.db\n")
      query = paste( "select gene_name from gtf",
         sprintf("where gene_biotype='%s' and moleculetype='%s'", biotype, moleculetype ) ,
         collapse = " " )
      genelist = dbGetQuery( obj@genome.db , query )
      genelist = genelist[,1]
   }

   if( is.null( tflist ) ) {
      if( verbose >= 1 )
          cat("no tf list is given. using all tfs from obj@project.db\n")
      query = "select distinct tf from motifsgenes"
      tflist = dbGetQuery( obj@project.db , query )[,1]
   }

   if( verbose >= 1 ) cat("fetching promoter regions for all genes\n")

   promoter_regions = 
      getPromoterRegionsAllGenes( obj )

   closeDatabaseConnections(obj)

   if( verbose >= 1 ) cat("counting binding sites for each tf in each region\n")

   tfbs_counts_per_gene = 
      getTfbsCountsPerGeneMC(
          tflist = tflist , promoter_regions = promoter_regions ,
          cores = cores , genome=genome , tissue=tissue )

   return( tfbs_counts_per_gene )


} #makeTfbsCountsTbl
#----------------------------------------------------------------------------------------------------
runTReNA <- function(gene, mtx, candidate.tfs)
{
  trena.lasso <- TReNA(mtx, solver="lasso")
  trena.bs    <- TReNA(mtx, solver="bayesSpike")
  trena.rf    <- TReNA(mtx, solver="randomForest")

  tfs <- intersect(candidate.tfs, rownames(mtx))
  printf("%8s: %d tf candidates", gene, length(tfs))
  tbl.lasso <- solve(trena.lasso, gene, tfs, extraArgs = list(alpha = 1.0))
  tbl.lasso <- subset(tbl.lasso, abs(beta) > 0.01)
  tbl.bs <- solve(trena.bs, gene, tfs)
  if(nrow(tbl.bs) > 0)
    tbl.bs <- subset(tbl.bs, pval < 0.05)

  suppressWarnings(
    rf.out <- solve(trena.rf, gene, tfs)
    )

  tbl.rf = rf.out$edges
  tbl.rf <-subset(tbl.rf, IncNodePurity >= fivenum(tbl.rf$IncNodePurity)[4])
  tbl.rf <- tbl.rf[order(tbl.rf$IncNodePurity, decreasing=TRUE),,drop=FALSE]

  tbl.lasso$gene <- rownames(tbl.lasso)
  tbl.lasso$method <- "lasso"
  tbl.lasso$score <- tbl.lasso$beta

  tbl.bs$gene <- rownames(tbl.bs)
  tbl.bs$method <- "bayesSpike"
  tbl.bs$score <- tbl.bs$beta

  tbl.rf$gene <- rownames(tbl.rf)
  tbl.rf$method <- "randomForest"
  tbl.rf$score <- tbl.rf$IncNodePurity

  tbl.all <- rbind.fill(tbl.lasso, tbl.bs, tbl.rf)
  tbl.all[order(abs(tbl.all$gene.cor), decreasing=TRUE),]

} # runTReNA
#----------------------------------------------------------------------------------------------------







