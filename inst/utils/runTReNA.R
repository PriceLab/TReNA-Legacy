#----------------------------------------------------------------------------------------------------
#' Get Footprints for a TF
#' @description Get all footprints from the project database for a specified transcription factor
#'
#' @param obj An object of class FootprintFilter (??)
#' @param tf A transcription factor

getFootprintsForTF <-
function( obj , tf ){

       motifsgenes = dbReadTable( obj@project.db , "motifsgenes" )
       fields = colnames(motifsgenes)
       
       if( "tf_name" %in% fields ) {

       query <-
paste( "select fp.chr, fp.mfpstart, fp.mfpend, fp.motifname, fp.pval, fp.score, mg.motif, mg.tf_name",
             "from footprints fp",
             "inner join motifsgenes mg",
             "on fp.motifName=mg.motif",
             paste("where mg.tf_name='",tf,"'" , sep="" ) ,  collapse = " " )

       }

       if( "tf" %in% fields ) {

       query <-
paste( "select fp.chr, fp.mfpstart, fp.mfpend, fp.motifname, fp.pval, fp.score, mg.motif, mg.tf",
             "from footprints fp",
             "inner join motifsgenes mg",
             "on fp.motifName=mg.motif",
             paste("where mg.tf='",tf,"'" , sep="" ) ,  collapse = " " )

       }

       if(!obj@quiet) print(query)
       dbGetQuery( obj@project.db , query )
} # getFootprintsForTF
#----------------------------------------------------------------------------------------------------
#' Get Gene Promoter Regions
#' @description Get the promoter regions for a list of genes. Region sizes can be tuned by specifying
#' the positions upstream and downstream of the gene that bound the region. 
#'
#' @param obj An object of class FootprintFilter(??)
#' @param genelist A gene list
#' @param size.upstream Number of basepairs to grab upstream of each gene (default = 10000)
#' @param size.downstream Number of basepairs to grab downstream of each gene (default = 10000)

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
#' Get Transcription Factor Binding Site Counts Per Promoter
#'
#' @description For a FootprintFilter object, get the transcription factor binding site counts for each promoter
#'
#' @param obj An object of class FootprintFilter
#' @param tflist A list of transcription factors
#' @param promotor_regions A list of promoter regions
#' @param verbose A logical indicating whether the function should produce verbose output (default = 1)
#'
#' @return A numeric vector containing transcription factor binding site counts per promoter 


getTfbsCountsPerPromoter <-
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

} #getTfbsCountsPerPromoter
#----------------------------------------------------------------------------------------------------
#' Get Transcription Factor Binding Site Counts Per Promoter with Clustering
#'
#' @description Get the transription factor binding site counts for each promoter and cluster
#'
#' @import doParallel
#' 
#' @param tflist A list of transcription factors
#' @param promoter_regions A list of promoter regions
#' @param verbose A logical indicating whether the function should produce verbose output (default = 1)
#' @param cores The number of computational cores to be devoted to the calculation (default = 5)
#' @param genome.db.uri A web address to the genome database
#' @param project.db.uri A web address to the project database
#'
#' @return A numeric vector containing transcription factor binding site counts per gene


getTfbsCountsPerPromoterMC <-
function( tflist , promoter_regions , verbose = 1 , cores = 5 , genome.db.uri , project.db.uri ) {

   cl = makeForkCluster( cores )
   registerDoParallel( cl )

   tfbs_counts_per_gene =
   foreach( x=1:length(tflist) ) %dopar% {
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
   stopCluster( cl )
   return( tfbs_counts_per_gene )

} #getTfbsCountsPerPromoterMC
#----------------------------------------------------------------------------------------------------
#' Get Hi-C Enhancer Regions
#'
#' Connect to the SQL database of hg38 enhancers and select only those regions where method = Hi-C
#' 
#' @return A table of enhancers where the method is Hi-C

getHiCEnhancerRegions <- function() {  

   enhancerdb <- 
      dbConnect(PostgreSQL(), user= "trena", password="trena", dbname="enhancers.hg38", host="whovian")
   query <- "select * from enhancers where method='Hi-C'"
   tbl = dbGetQuery( enhancerdb, query )
   return( tbl )
}
#----------------------------------------------------------------------------------------------------
#' Get DNase Enhancer Regions
#'
#' Connect to the SQL database of hg38 enhancers and select only those regions where method = DNase-DNase
#'
#' @return A table of enhancers where the method is DNase-DNase

getDnaseEnhancerRegions <- function() {   

   enhancerdb <- 
      dbConnect(PostgreSQL(), user= "trena", password="trena", dbname="enhancers.hg38", host="whovian")
   query <- "select * from enhancers where method='DNase-DNase'"
   tbl = dbGetQuery( enhancerdb, query )
   return( tbl )
}
#----------------------------------------------------------------------------------------------------
#' Get Gene List from Gtf
#'
#' Query the SQL gtf table to get a list of genes matching the supplied "biotype" and "moleculetype"
#'
#' @param genome.db.uri A web address to a genome database
#' @param project.db.uri A web address to a project database
#' @param biotype A string specifying the biotype of interest (default = "protein_coding") **What are the options?
#' @param moleculetype A string specifying the moleculetype of interest (default = "gene") **What are the options?
#' @param use_gene_ids A logical character indicating whether or not to use gene IDs as opposed to gene names (default = T)
#'
#' @return A list of genes resulting from the query

getGenelistFromGtf <- function( genome.db.uri , project.db.uri , 
   biotype="protein_coding" , moleculetype="gene" , use_gene_ids = T ) {

      obj <- FootprintFinder(genome.db.uri, project.db.uri, quiet=TRUE)

      if( use_gene_ids == T ) {
         query = paste( "select gene_id from gtf",
            sprintf("where gene_biotype='%s' and moleculetype='%s'", biotype, moleculetype ) ,
            collapse = " " )
      }
      if( use_gene_ids == F ) {
         query = paste( "select gene_name from gtf",
            sprintf("where gene_biotype='%s' and moleculetype='%s'", biotype, moleculetype ) ,
            collapse = " " )
      }
      genelist = dbGetQuery( obj@genome.db , query )
      genelist = genelist[,1]
      closeDatabaseConnections(obj)
      return( genelist )

} # getGenelistFromGtf
#----------------------------------------------------------------------------------------------------
#' Get Transcription Factor Binding Site Counts in Enhancers
#'
#' @import doParallel

getTfbsCountsInEnhancers <-
function( tflist = NULL , genelist = NULL , enhancertype = "Hi-C" , verbose = 2 , 
   cores = 5 , genome.db.uri , project.db.uri , use_gene_ids = T ,
   biotype = "protein_coding" , moleculetype="gene" ) {


   if( enhancertype == "Hi-C" )
      enhancers = getHiCEnhancerRegions()
   if( enhancertype == "DNase" )
      enhancers = getDnaseEnhancerRegions()

   if( use_gene_ids == T )
      enhancers = enhancers[,c("chr2","start2","end2","geneid")]
   if( use_gene_ids == F )
      enhancers = enhancers[,c("chr2","start2","end2","genename")]
   enhancers = unique( enhancers )
   colnames(enhancers) = c("chr","start","end","gene")
   enhancers = makeGRangesFromDataFrame( enhancers , keep.extra.columns = T )

   if( is.null( genelist )) {
      if( verbose >= 1 )
          cat( "no gene list is given, using all genes from obj@genome.db\n" )
      genelist = getGenelistFromGtf( genome.db.uri=genome.db.uri , project.db.uri=project.db.uri , use_gene_ids = use_gene_ids )
   }

   if( is.null( tflist ) ) {
      if( verbose >= 1 )
          cat("no tf list is given. using all tfs from obj@project.db\n")
      obj <- FootprintFinder(genome.db.uri, project.db.uri, quiet=TRUE)
      motifsgenes = dbReadTable( obj@project.db , "motifsgenes" )
      if( "tf_name" %in% colnames(motifsgenes) )
         tflist = unique( motifsgenes$tf_name )
      if( "tf" %in% colnames(motifsgenes) )
         tflist = unique( motifsgenes$tf )
      closeDatabaseConnections(obj)
   }

   if( verbose >= 1 ) cat("counting TFBSs per enhancer\n")


   cl <- makeForkCluster(cores)
   registerDoParallel( cl )

   tfbs_counts_per_enhancer =
   foreach( x=1:length(tflist) ) %dopar% {
      obj <- FootprintFinder(genome.db.uri, project.db.uri, quiet=TRUE)
      footprints.tf = getFootprintsForTF( obj , tflist[x] )
      closeDatabaseConnections( obj )

      if( nrow( footprints.tf ) == 0 ) {
         counts = rep( 0 , length(enhancers) )
         return(counts)
      }
      colnames( footprints.tf )[1:3] = c("chr","start","end")
      footprints.tf.gr = makeGRangesFromDataFrame( footprints.tf )
      fp.reduced = reduce( footprints.tf.gr )
      counts = countOverlaps( enhancers , fp.reduced )
      if( verbose > 1 & x/10 == round(x/10) ) cat("...done" , x , "/" , length(tflist) ,"\n" )
      return( counts )
   }
   stopCluster(cl)
   tfbs_counts_per_enhancer = do.call( cbind , tfbs_counts_per_enhancer )
   colnames(tfbs_counts_per_enhancer) = tflist

   enhancertargets = enhancers$gene

   if( verbose >= 1 ) cat("counting TFBSs per gene\n")

   cl = makeForkCluster( cores )
   registerDoParallel( cl )

   tfbs_counts_per_gene =
   foreach( gene=genelist ) %dopar% {
      enhancers_per_gene = which( enhancertargets == gene )
      if( length(enhancers_per_gene) == 0 ) return( rep(0,length(tflist)) )
      tmp = tfbs_counts_per_enhancer[ which( enhancertargets == gene ) , , drop = F ]
      counts = colSums( tmp )
      return(counts)
   }
   tfbs_counts_per_gene = do.call( rbind , tfbs_counts_per_gene )
   colnames(tfbs_counts_per_gene) = tflist
   rownames(tfbs_counts_per_gene) = genelist

   stopCluster( cl )

   if( use_gene_ids == T ) {
      cat("converting TF names to ENSG ids\n")

      obj <- FootprintFinder(genome.db.uri, project.db.uri, quiet=TRUE)
      query = paste( "select gene_id, gene_name from gtf",
         sprintf("where gene_biotype='%s' and moleculetype='%s'", biotype, moleculetype ) ,
         collapse = " " )
      gene_ids_names = dbGetQuery( obj@genome.db , query )
      matchRows = match( colnames(tfbs_counts_per_gene) , gene_ids_names$gene_name )
      tf_ids = gene_ids_names[ matchRows , "gene_id" ]
      colnames(tfbs_counts_per_gene) = tf_ids
      na_tf_ids = is.na( tf_ids )
      tfbs_counts_per_gene = tfbs_counts_per_gene[ , which(na_tf_ids ==F) ]
      
   }

   return( tfbs_counts_per_gene )

} # getTFbsCountsInEnhancers
#----------------------------------------------------------------------------------------------------
getTfbsCountsInPromoters <- 
function( genome.db.uri , project.db.uri ,  genelist = NULL , tflist = NULL , 
   biotype = "protein_coding" , moleculetype = "gene" , use_gene_ids = T ,
   size.upstream = 10000 , size.downstream = 10000 , cores = 1 , verbose = 1 ) {

   if( is.null( genelist )) {
      if( verbose >= 1 )
          cat( "no gene list is given, using all genes from obj@genome.db\n" )
      genelist = getGenelistFromGtf( genome.db.uri=genome.db.uri , project.db.uri=project.db.uri , use_gene_ids = use_gene_ids )
   }

   if( is.null( tflist ) ) {
      if( verbose >= 1 )
          cat("no tf list is given. using all tfs from obj@project.db\n")
      obj <- FootprintFinder(genome.db.uri, project.db.uri, quiet=TRUE)
      motifsgenes = dbReadTable( obj@project.db , "motifsgenes" )
      if( "tf_name" %in% colnames(motifsgenes) )
         tflist = unique( motifsgenes$tf_name )
      if( "tf" %in% colnames(motifsgenes) )
         tflist = unique( motifsgenes$tf )
      closeDatabaseConnections(obj)
   }

   if( verbose >= 1 ) cat("fetching promoter regions for all genes\n")

      obj <- FootprintFinder(genome.db.uri, project.db.uri, quiet=TRUE)
      promoter_regions = 
         getPromoterRegionsAllGenes( obj , use_gene_ids = use_gene_ids )
      closeDatabaseConnections(obj)

   if( verbose >= 1 ) cat("counting binding sites for each tf in each region\n")

   tfbs_counts_per_gene = 
      getTfbsCountsPerPromoterMC(
          tflist = tflist , promoter_regions = promoter_regions ,
          cores = cores , genome.db.uri=genome.db.uri , project.db.uri=project.db.uri )

   if( use_gene_ids == T ) {
      cat("converting TF names to ENSG ids\n")

      obj <- FootprintFinder(genome.db.uri, project.db.uri, quiet=TRUE)
      query = paste( "select gene_id, gene_name from gtf",
         sprintf("where gene_biotype='%s' and moleculetype='%s'", biotype, moleculetype ) ,
         collapse = " " )
      gene_ids_names = dbGetQuery( obj@genome.db , query )
      matchRows = match( colnames(tfbs_counts_per_gene) , gene_ids_names$gene_name )
      tf_ids = gene_ids_names[ matchRows , "gene_id" ]
      colnames(tfbs_counts_per_gene) = tf_ids
      na_tf_ids = is.na( tf_ids )
      tfbs_counts_per_gene = tfbs_counts_per_gene[ , which(na_tf_ids == F ) ]
      
   }

   return( tfbs_counts_per_gene )

} #getTfbsCountsInPromoters
#----------------------------------------------------------------------------------------------------
runAllSolversForGene <- function(gene, mtx, candidate.tfs)
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

} # runAllSolversForGene
#----------------------------------------------------------------------------------------------------
#' Make a TRN from promoter counts and expression
#'
#' @import doParallel

makeTrnFromPromoterCountsAndExpression <- 
function( counts , expr , method = "lasso" , alpha = 0.5 ,
   cores = NULL , center_and_scale = T , candidate_regulator_method = "quantile" ,
   tfbs.quantile.thresh = 0.75 ){   

   stopifnot( method %in% c("lasso","randomForest") )

   # prepare expression data
   expr = as.matrix(expr)
   if( center_and_scale == T ) {
      gene.mean = rowMeans(expr)
      gene.sd = apply( expr , 1 , stats::sd )
      expr = ( expr - gene.mean ) / gene.sd
   }

   cat("selecting genes for which we have both expression data and tfbs counts\n")
   isec.tfs = intersect( colnames(counts) , rownames(expr) )
   isec.genes = intersect( rownames(counts) , rownames(expr) )
   stopifnot( length(isec.genes) > 0 && length(isec.tfs) > 0 )
   tfbs = counts[ isec.genes , isec.tfs ]
   expr = expr[ isec.genes , ]

   stopifnot( all( rownames(tfbs) == rownames(expr) ) )

   # set up multi-core operation
   if( is.null(cores) ) cores = round( detectCores() / 3 )
   cl = makeForkCluster( cores )
   registerDoParallel( cl )

   # select candidate regulators for each gene based on tfbs

   if( candidate_regulator_method == "quantile" ) {
      if( is.null( tfbs.quantile.thresh ) )  tfbs.quantile.thresh = 0.75
      tfbs.quantile = apply( tfbs , 2 , quantile , probs = tfbs.quantile.thresh )
      candidate_regulators = t( t(tfbs) > tfbs.quantile )
   }
   if( candidate_regulator_method == "all" ) {
      candidate_regulators = tfbs
   }



   if( method == "lasso" ) {

      # set up the TReNA object
      trena = TReNA(mtx.assay=expr,solver=method)

      cat("determining an appropriate value for the panalty parameter lambda\n")
      fit.cv =
      foreach( target.gene=sample(rownames(expr),100) ) %dopar% {
         tfs = names(which( candidate_regulators[target.gene,] > 0 ))
         fit = solve(trena,target.gene,tfs,extraArgs=list(alpha=alpha,keep.metrics=T))
      }

      lambda = do.call( c ,
         lapply(1:length(fit.cv), function(i) fit.cv[[i]]$lambda))
      lambda.median = median(lambda,na.rm=T)
 
      cat("fiting model for all genes using the median lambda from fit.cv\n")

      fit2 =
      foreach( target.gene=rownames(expr) ) %dopar% {
         # tfs = names(which(candidate_regulators[target.gene,]==T))
         tfs = names(which( candidate_regulators[target.gene,] > 0 ))
         fit = solve(trena,target.gene,tfs,
             extraArgs=list(alpha=alpha,lambda=lambda.median,keep.metrics=T))
         if( length(fit) > 0 ) {
            if( nrow(fit$mtx.beta) > 0 ) {
               fit$mtx.beta$target = target.gene
               fit$mtx.beta$tf = rownames(fit$mtx.beta)
            }
         }
         return( fit )
      }

      r2 = do.call( c ,
         lapply(1:length(fit2), function(i) fit2[[i]]$r2))
      n.nonzero = do.call( c ,
         lapply(1:length(fit2), function(i) nrow(fit2[[i]]$mtx.beta)))
      trn = do.call( rbind ,
         lapply(1:length(fit2), function(i) fit2[[i]]$mtx.beta))

   }

   if( method == "randomForest" ) {

      trena <- TReNA(mtx.assay=enorm, solver="randomForest", quiet=FALSE)
   
      trn0 =
      foreach( target.gene=rownames(expr) ) %dopar% {
         tfs <- names(which( candidate_regulators[target.gene,] > 0 ))
         result <- solve(trena, target.gene, tfs)
         if( is.null(result) == F ) {
            result$edges$target = target.gene
            result$edges$tf = rownames(result$edges)
         }
         return(result)
      }

      r2 = do.call( c ,
         lapply(1:length(trn0), function(i) trn0[[i]]$r2))
      trn = do.call( rbind ,
         lapply(1:length(trn0), function(i) trn0[[i]]$edges))

   }

   stopCluster( cl )

   return( list( trn = trn , r2 = r2 ))

} # makeTrnFromPromoterCountsAndExpression






