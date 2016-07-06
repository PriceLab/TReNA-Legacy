library(TReNA)

# 0. specify the tissue and genome for analysis
genome = "hg38"
tissue = "lymphoblast"
outdir = "/proj/price1/sament/TReNA/inst/extdata/"
cores = 5
verbose = 2


# 1, get counts of binding sites for each TF proximal to each gene (by default +/- 10kb from the TSS)

promoter_counts = getTfbsCountsInPromoters( 
   genome=genome , tissue=tissue , # specify the genome build and tissue source for footprints
   size.upstream = 10000 , size.downstream = 10000 , # define the size of the window around the TSS
   cores = cores , verbose = verbose ) # specify the number of cores for parallelization

save( promoter_counts , 
   file = paste( outdir , "promoter_tfbs_counts." , genome , "." , tissue , ".RData" , sep="" ))
  
# 2, get counts of binding sites for each TF in each gene's enhancers

# 2a. use enhancer-promoter loops from Hi-C (Rao et al. 2015)

enhancer_counts_hic = getTfbsCountsInEnhancers(
   genome=genome , tissue=tissue , # specify the genome build and tissue source for footprints
   enhancertype = "Hi-C" , # source of enhancer-promoter loops (Hi-C vs. DNase-Dnase correlation)
   cores = cores , verbose = verbose )

save( enhancer_counts_hic , 
   file = paste( outdir , "enhancer_tfbs_counts_hic." , genome , "." , tissue , ".RData" , sep="" ))

# 2b. use enhancer-promoter correlations from DNase-seq (Thurman et al. 2012)

enhancer_counts_dnase = getTfbsCountsInEnhancers(
   genome=genome , tissue=tissue , # specify the genome build and tissue source for footprints
   enhancertype = "DNase" , # source of enhancer-promoter loops (Hi-C vs. DNase-Dnase correlation)
   cores = cores , verbose = verbose ) 
 
save( enhancer_counts_dnase , 
   file = paste( outdir , "enhancer_tfbs_counts_dnase." , genome , "." , tissue , ".RData" , sep="" ))




