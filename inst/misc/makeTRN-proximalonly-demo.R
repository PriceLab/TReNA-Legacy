library(TReNA)

# 0. set options
genome = "hg38"
tissue = "lymphoblast"
outdir = "/proj/price1/sament/TReNA/inst/extdata/"
cores = 5
verbose = 2
method = "lasso"
print(load("/proj/price1/sament/TReNA/inst/extdata/GSE37772.expr.RData"))
expr = expr2
rm( expr2 )
# 1, get counts of binding sites for each TF proximal to each gene (by default +/- 10kb from the TSS)

promoter_counts = getTfbsCountsInPromoters( 
   genome=genome , tissue=tissue , # specify the genome build and tissue source for footprints
   size.upstream = 10000 , size.downstream = 10000 , # define the size of the window around the TSS
   cores = cores , verbose = verbose ) # specify the number of cores for parallelization

save( promoter_counts , 
   file = paste( outdir , "promoter_tfbs_counts.gene_ids." , genome , "." , tissue , ".RData" , sep="" ))
  
# 2. build TRN model by integrating TFBS counts with expression data

trn = makeTrnFromPromoterCountsAndExpression( counts = promoter_counts , expr = expr , 
   method = method , cores = cores )

save( trn , 
   file = paste( outdir , "trn." , genome , "." , tissue , ".RData" , sep = "" ))


