library(TReNA)
library(RUnit)
#----------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_getFootprintsForTF()
   test_getGenePromoterRegions()
   test_getTfbsCountsPerGene()
   test_makeTfbsCountsTbl()
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
test_getTfbsCountsPerGene <- function()
{
   printf("--- test_getTfbsCountsPerGene")

   genome.db.uri <- "postgres://whovian/hg38"
   project.db.uri <-  "postgres://whovian/lymphoblast"
   fp <- FootprintFinder(genome.db.uri, project.db.uri, quiet=TRUE)

   promoter_regions = test_getGenePromoterRegions(quiet=T)

   # get TF names
      query = "select distinct tf from motifsgenes"
      tflist = dbGetQuery( fp@project.db , query )[,1]
   tflist = sample( tflist , 10 )

   tfbs_counts = getTfbsCountsPerGene( fp , tflist , promoter_regions = promoter_regions )

   checkTrue( sum( tfbs_counts ) > 0 )
   checkTrue( all( colnames(tfbs_counts) == tflist ))

} # test_getTfbsCountsPerGene
#----------------------------------------------------------------------------------------------------
test_makeTfbsCountsTbl <- function()
{
   printf("--- test_makeTfbsCountsTbl")

   tfbs_counts = makeTfbsCountsTbl( genome="hg38" , tissue="lymphoblast" , 
      tflist = c("RXRA","SATB2","EMX2","SP1","SP2") , cores = 2 , verbose = 2 )

   checkTrue( max( apply( tfbs_counts , 2 , max ) ) < 100 )

}
#----------------------------------------------------------------------------------------------------
if(!interactive()) runTests()









