# runTReNA-demo.R: a simple way to run all current solvers, naively combining results into a
# score for each TF candidate meeting threshold.
#
library(TReNA)
library(plyr)

stopifnot(packageVersion("TReNA") >= '0.99.40')
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
# our demo gene is MEF2C, for which a small rna expression matrix is included
# in the TReNA package:
#
#           gene_id gene_name  chr    start   endpos strand
# 1 ENSG00000081189     MEF2C chr5 88717117 88904257      -
#
goi <- "MEF2C"
start <- 88717117 - 1000
end   <- 88904257 + 1000
chrom <- "chr5"
print(load(system.file(package="TReNA", "extdata", "ampAD.154genes.mef2cTFs.278samples.RData")))
stopifnot(dim(mtx.sub) == c(154, 278))

   # just keep rows with some variance across samples
sd <- apply(mtx.sub, 1, sd)
deleters <- which(sd < 1)
if(length(deleters) > 0)
  mtx.sub <- mtx.sub[-deleters,]
mtx.tmp <- mtx.sub - min(mtx.sub) + 0.001
mtx.sub.log2 <- log2(mtx.tmp)
fivenum(mtx.sub.log2)
genes <- rownames(mtx.sub.log2)

   # no TF self-loops, but target must be in the expression matrix
stopifnot("MEF2C" %in% genes)
tf.candidates <- genes[-grep("MEF2C", genes)]
tbl.all <- runTReNA(goi, mtx.sub.log2, tf.candidates)
print(head(tbl.all, n=10))
   # compare the methods
print(head(table(tbl.all$gene, tbl.all$method)))
