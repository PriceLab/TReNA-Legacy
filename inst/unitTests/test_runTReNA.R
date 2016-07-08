library(TReNA)
library(RUnit)
#----------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_getFootprintsForTF()
   test_getGenePromoterRegions()
   test_getTfbsCountsPerPromoter()
   test_getTfbsCountsInPromoters()
   test_makeTrnFromPromoterCountsAndExpression()
   test_makeTrnFromPromoterCountsAndExpression.useAllTFs()
} # runTests
#----------------------------------------------------------------------------------------------------
test_getFootprintsForTF = function()
{
   printf("--- test_getFootprintsForTF")

   genome.db.uri <- "postgres://whovian/hg38"
   project.db.uri <-  "postgres://whovian/lymphoblast"
   fp <- FootprintFinder(genome.db.uri, project.db.uri, quiet=TRUE)

   tf = "RXRA"

   footprints = getFootprintsForTF( obj = fp , tf = tf )

   checkTrue( nrow(footprints) > 10000 )


} #test_getFootprintsForTF
#----------------------------------------------------------------------------------------------------
test_getGenePromoterRegions = function(quiet=F)
{
   if( quiet==F ) printf("--- test_getGenePromoterRegions")

   genome.db.uri <- "postgres://whovian/hg38"
   project.db.uri <-  "postgres://whovian/lymphoblast"
   fp <- FootprintFinder(genome.db.uri, project.db.uri, quiet=TRUE)

   # get gene_name field from gtf
      query = paste( "select gene_name from gtf",
         sprintf("where gene_biotype='%s' and moleculetype='%s'", "protein_coding", "gene" ) ,
         collapse = " " )
      genelist = dbGetQuery( fp@genome.db , query )
      genelist = genelist[,1]
      genelist_sample = sample( genelist , 5 )

  promoter_regions = getGenePromoterRegions( fp , genelist_sample )

  checkTrue( length( promoter_regions ) == 5 )
  checkTrue( all( width(ranges(promoter_regions)) == 20001 ))

  invisible( promoter_regions )

} #test_getGenePromoterRegions
#----------------------------------------------------------------------------------------------------
test_getTfbsCountsPerPromoter <- function()
{
   printf("--- test_getTfbsCountsPerPromoter")

   genome.db.uri <- "postgres://whovian/hg38"
   project.db.uri <-  "postgres://whovian/lymphoblast"
   fp <- FootprintFinder(genome.db.uri, project.db.uri, quiet=TRUE)

   promoter_regions = test_getGenePromoterRegions(quiet=T)

   # get TF names
      query = "select distinct tf from motifsgenes"
      tflist = dbGetQuery( fp@project.db , query )[,1]
   tflist = sample( tflist , 10 )

   tfbs_counts = getTfbsCountsPerPromoter( fp , tflist , promoter_regions = promoter_regions )

   checkTrue( sum( tfbs_counts ) > 0 )
   checkTrue( all( colnames(tfbs_counts) == tflist ))

} # test_getTfbsCountsPerPromoter
#----------------------------------------------------------------------------------------------------
test_getTfbsCountsInPromoters <- function()
{
   printf("--- test_getTfbsCountsInPromoters")

   tfbs_counts = getTfbsCountsInPromoters( genome.db.uri = "postgres://whovian/hg38" , 
      project.db.uri = "postgres://whovian/lymphoblast" , 
      tflist = c("RXRA","NR3C3","SATB2","EMX2","SP1","SP2") , cores = 2 , verbose = 2 )


   checkTrue( ncol(tfbs_counts) == 5 )
   checkTrue( max( apply( tfbs_counts , 2 , max ) ) < 100 )
   checkTrue( colnames(tfbs_counts)[1] == "ENSG00000186350" )

} # test_getTfbsCountsInPromoters
#----------------------------------------------------------------------------------------------------
test_getTfbsCountsInEnhancers <- function()
{
   printf("--- test_getTfbsCountsInPromoters")

   tfbs_counts = getTfbsCountsInEnhancers( genome.db.uri = "postgres://whovian/hg38" ,
      project.db.uri = "postgres://whovian/lymphoblast" , 
      tflist = c("RXRA","NR3C3","SATB2","EMX2","SP1","SP2") , cores = 2 )

   checkTrue( any( colSums( tfbs_counts ) > 0 ))

}
#----------------------------------------------------------------------------------------------------
test_makeTrnFromPromoterCountsAndExpression <- function()
{
   printf("--- test_makeTrnFromPromoterCountsAndExpression")
 

   print(load(system.file(package="TReNA", "extdata/promoter_tfbs_counts.gene_ids.hg38.lymphoblast.RData")))

   print(load(system.file(package="TReNA","extdata/GSE37772.expr.RData")))

   trn = makeTrnFromPromoterCountsAndExpression( 
      counts = promoter_counts , expr = expr2 , method = "lasso" )

   edges = trn$trn
   r2 = trn$r2
   
   checkTrue( nrow( edges ) > 100000 )
   checkTrue( median(r2 , na.rm = T ) > 0.2 )
} # test_makeTrnFromPromoterCountsAndExpression
#----------------------------------------------------------------------------------------------------
test_makeTrnFromPromoterCountsAndExpression.useAllTFs <- function()
{
   printf("--- test_makeTrnFromPromoterCountsAndExpression")


   print(load(system.file(package="TReNA", "extdata/promoter_tfbs_counts.gene_ids.hg38.lymphoblast.RData")))

   print(load(system.file(package="TReNA","extdata/GSE37772.expr.RData")))
 
   trn2 = makeTrnFromPromoterCountsAndExpression(
      counts = promoter_counts , expr = expr2 , method = "lasso" , candidate_regulator_method = "all" )

   edges = trn2$trn
   r2 = trn2$r2

   checkTrue( nrow( edges ) > 100000 )
   checkTrue( median(r2 , na.rm = T ) > 0.2 )

   kIn = table( edges$target )
   kOut = table( edges$tf )

   checkTrue( median( kIn ) > 10 & median( kIn ) < 25 )
   checkTrue( max( kOut ) < 3000 )


} # test_makeTrnFromPromoterCountsAndExpression.useAllTFs
#----------------------------------------------------------------------------------------------------
if(!interactive()) runTests()









