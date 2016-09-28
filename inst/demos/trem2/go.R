library(RPostgreSQL)
library(TReNA)
#------------------------------------------------------------------------------------------------------------------------
source("../../../../BDDS/trenadb/utils.R")
#------------------------------------------------------------------------------------------------------------------------
if(!exists("db.trena"))
   db.trena <- dbConnect(PostgreSQL(), user="trena", password="trena", dbname="trena", host="whovian")

if(!exists("db.gtf"))
   db.gtf <- dbConnect(PostgreSQL(), user= "trena", password="trena", dbname="gtf", host="whovian")

if(!exists("trem2")){
   tbl.tmp <- dbGetQuery(db.gtf, "select * from hg38human where gene_name='TREM2' and moleculetype='gene'")
   trem2 <- list(chrom=tbl.tmp[1, "chr"], start=tbl.tmp[1, "start"])
   }

if(!exists("db.hint"))
   db.hint <- dbConnect(PostgreSQL(), user="trena", password="trena", dbname="hint", host="whovian")

if(!exists("db.wellington"))
   db.wellington <- dbConnect(PostgreSQL(), user="trena", password="trena", dbname="wellington", host="whovian")

if(!exists("tbl.genesmotifs"))
    tbl.genesmotifs <- dbGetQuery(db.trena, "select * from tfmotifs")

#------------------------------------------------------------------------------------------------------------------------
chrom <- trem2$chrom
start <- trem2$start - 1000
end   <- trem2$start + 1000
target.gene <- "TREM2"

tbl.h <- createHintTable(chrom, start, end)
motifs <- unique(tbl.h$motif.h)
candidate.tfs <- sort(unique(subset(tbl.genesmotifs, motif %in% tbl.h$motif)$gene))

if(!exists("mtx.tcx"))   # "mtx.tcx"     "mtx.tcx.ctl" "mtx.tcx.ad" 
   print(load("~/github/Private_Cory_Data/inst/extdata/prepped.tcx.matrices.RData"))


matrices <- list(all=mtx.tcx, ad=mtx.tcx.ad, ctl=mtx.tcx.ctl)

results <- lapply(matrices, function(mtx){
   stopifnot(target.gene %in% rownames(mtx))
   candidate.regulators <- intersect(candidate.tfs, rownames(mtx))
   genes.of.interest <- c(target.gene, candidate.regulators)
   mtx.sub <- mtx[genes.of.interest,]
   mtx.adjusted <- asinh(mtx.sub)
   trena <- TReNA(mtx.assay=mtx.adjusted, solver="lasso", quiet=FALSE)
   solve(trena, target.gene, candidate.regulators, extraArgs =list(alpha=0.1, lambda=NULL))
   })



    
