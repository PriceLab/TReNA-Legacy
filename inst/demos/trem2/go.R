library(RPostgreSQL)
library(TReNA)
library (RColorBrewer)
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
target.gene <- "TREM2"
shoulder <- 1000
chrom <- trem2$chrom
start <- trem2$start - shoulder
end   <- trem2$start + shoulder

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
   tbl <- solve(trena, target.gene, candidate.regulators, extraArgs=list(alpha=0.1, lambda=NULL))
   })

results2 <- lapply(names(results), function(name){
    column.names <- colnames(results[[name]])
    new.names <- paste(name, column.names, sep=".")
    colnames(results[[name]]) <- new.names
    results[[name]]
    })

names(results2) <- names(results)

all.tfs <- sort(unique(c(rownames(results2[[1]]),
                         rownames(results2[[2]]),
                         rownames(results2[[3]]))))
all.conditions <- c(colnames(results2[[1]]),
                    colnames(results2[[2]]),
                    colnames(results2[[3]]))
m <- matrix(0, nrow=length(all.tfs), ncol=length(all.conditions), dimnames=list(all.tfs, all.conditions))
for(i in 1:length(results2)){
    mtx <- as.matrix(results2[[i]])
    row.names <- rownames(mtx)
    column.names <- colnames(mtx)
    print(row.names)
    print(column.names)
    print(mtx[row.names, column.names])
    print(m[row.names, column.names])
    m[row.names, column.names] <- mtx
    }

m.beta <- m[, grep("beta", colnames(m))]
m.beta.noAll <- m.beta[, -(grep("all", colnames(m.beta)))]
heatmap.2(m.beta, margins=c(10,10), trace='none')
heatmap.2(m.beta.noAll, margins=c(10,10), trace='none', cexCol=2)
