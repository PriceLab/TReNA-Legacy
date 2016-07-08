library(TReNA)

# 0. set options
genome.db.uri = "postgres://whovian/hg38"
project.db.uri = "postgres://whovian/lymphoblast"
out = "/proj/price1/sament/lymphoblast_trn/hg38.lymphoblast"
cores = 5
method = "lasso"
print(load(system.file(package="TReNA","/extdata/GSE37772.expr.RData")))
expr = expr2
rm( expr2 )


# 1, get counts of binding sites for each TF proximal to each gene (by default +/- 10kb from the TSS)

promoter_counts = getTfbsCountsInPromoters( 
   genome.db.uri=genome.db.uri , project.db.uri=project.db.uri , 
   size.upstream = 10000 , size.downstream = 10000 , # define the size of the window around the TSS
   cores = cores ) # specify the number of cores for parallelization

save( promoter_counts , 
   file = paste( out , ".promoter_tfbs_counts.gene_ids.RData" , sep="" ))
  
# 2, get counts of binding sites for each TF in each gene's enhancers

# 2a. use enhancer-promoter loops from Hi-C (Rao et al. 2015)

enhancer_counts_hic = getTfbsCountsInEnhancers(
   genome.db.uri=genome.db.uri , project.db.uri=project.db.uri , 
   enhancertype = "Hi-C" , # source of enhancer-promoter loops (Hi-C vs. DNase-Dnase correlation)
   cores = cores )

save( enhancer_counts_hic , 
   file = paste( out , ".enhancer_tfbs_counts_hic.RData" , sep="" ))

# 2b. use enhancer-promoter correlations from DNase-seq (Thurman et al. 2012)

enhancer_counts_dnase = getTfbsCountsInEnhancers(
   genome.db.uri=genome.db.uri , project.db.uri=project.db.uri ,   
   enhancertype = "DNase" , # source of enhancer-promoter loops (Hi-C vs. DNase-Dnase correlation)
   cores = cores ) 
 
save( enhancer_counts_dnase , 
   file = paste( out , ".enhancer_tfbs_counts_dnase.gene_ids.RData" , sep="" ))









