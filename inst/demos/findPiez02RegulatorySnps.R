library(PrivateCoryData)
library(TReNA)
library(GenomicRanges)

stopifnot(packageVersion("PrivateCoryData") >= "1.0.18")
stopifnot(packageVersion("TReNA") >= '0.99.33')
stopifnot(packageVersion("GenomicRanges") >= '1.22.4')

#library(igvR)
#library(snpFoot)
#stopifnot(packageVersion("igvR") >= '0.99.23')
#stopifnot(packageVersion("snpFoot") >= '1.0.45')

goi <- "LINC01444"   # not enough expression variance to use
goi <- "PIEZO2"
genome.db.uri <- "postgres://whovian/hg38"
project.db.uri <-  "postgres://whovian/wholeBrain"
fp <- FootprintFinder(genome.db.uri, project.db.uri, quiet=TRUE)
# getGtfGeneBioTypes(fp) to see range of values
# getGtfMoleculeTypes(fp) to see range of values
# make sure goi is recognized.  as a side effect (and needed here) we also learn the ENSG for LINC01444
print(getChromLoc(fp, goi, biotype="protein_coding", moleculetype="gene"))
#           gene_id gene_name   chr    start   endpos strand
# 1 ENSG00000264301 LINC01444 chr18 14969001 14970468      -
goi <- "ENSG00000154864"
tbl.fp <- getFootprintsForGene(fp, goi, size.upstream=3000,  size.downstream=1000, biotype="protein_coding")

candidate.tfs.ensg <- sort(unique(tbl.fp$tf_ensg))  #  355

print(load(system.file(package="PrivateCoryData", "extdata", "ampADMayo.64253genes.278samples.RData")))
  # just keep rows with some variance across samples
sd <- apply(mtx, 1, sd)
deleters <- which(sd < 1)
if(length(deleters) > 0)
    mtx.sub <- mtx[-deleters,]

mtx.tmp <- mtx.sub - min(mtx.sub) + 0.001
mtx.sub.log2 <- log2(mtx.tmp)
fivenum(mtx.sub.log2)  # [1] -9.9657843  0.8107618  3.6262014  5.4345771 10.0052973
trena <- TReNA(mtx.assay=mtx.sub.log2, solver="lasso", quiet=FALSE)

stopifnot(goi %in% rownames(mtx.sub.log2))

filtered.candidate.tfs <- intersect(candidate.tfs.ensg, rownames(mtx.sub.log2))  # 143/355

tbl.lasso <- solve(trena, goi, filtered.candidate.tfs)
tbl.lasso <- subset(tbl.lasso, abs(beta) > 0.1)   # 14 tfs in model

# quick hack to get gene names
query <- sprintf("select gene_id, gene_name from gtf where gene_id in ('%s')", paste(rownames(tbl.lasso), collapse="','"))
tbl.ids <- unique(dbGetQuery(fp@genome.db, query))
order <- match(rownames(tbl.lasso), tbl.ids$gene_id)
name <- tbl.ids$gene_name[order]
tbl.lasso$name <- name

# now extract only those footprints (20 out of 3693) associated with tfs which have good lasso betas,
# predicting the expression of the target gene.
# tbl.fpoi:  "footprints of interest"
tbl.fpoi <- subset(tbl.fp, tf_ensg %in% rownames(tbl.lasso))

# intersect these 20 footprints with snps in eblue's bed file, using GRanges

full.path <- system.file(package="PrivateCoryData", "extdata", "elizabethBlue", "chr18_CU0070F.shared.hg38.bed")
stopifnot(file.exists(full.path))
tbl.bed <- read.table(full.path, sep="\t", as.is=TRUE)
colnames(tbl.bed) <- c("chr", "start", "end", "signature");

gr.fp <- with(tbl.fpoi, GRanges(seqnames=chr, IRanges(start=mfpstart, end=mfpend)))           # 20 ranges
padding <- 10
gr.snp <- with(tbl.bed, GRanges(seqnames=chr, IRanges(start=start-padding, end=end+padding)))  # 1366 ranges
tbl.overlaps <- as.data.frame(findOverlaps(gr.fp, gr.snp))
# create a new "snp" column to add to tbl.fpoi, with initial value for each of -1
snp <- rep(-1, length(gr.fp))
# now look for snps which fall in any of the footprints, and record the location of that snp,
# overwiting the -1 with the snp's signature
snp[tbl.overlaps$queryHits] <- tbl.bed[tbl.overlaps$subjectHits, "signature"]

tbl.fpoi$snp <- snp
print(subset(tbl.fpoi, snp != -1))
