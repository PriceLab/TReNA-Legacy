library(TReNA)

db <- dbConnect(PostgreSQL(), host="whovian", user="trena", password="trena", dbname="hg38")
target.gene <- "MEF2C"
biotype <- "protein_coding"
moleculetype <- "gene"
query <- paste("select gene_name, chr, start, endpos, strand from gtf where ",
                 sprintf("gene_name='%s' ", target.gene),
                 sprintf("and gene_biotype='%s' and moleculetype='%s'", biotype, moleculetype),
                 collapse=" ")


tbl.target <- dbGetQuery(db, query)
anchor.loc <- tbl.target[1, "start"]   # note: on reverse strand, so start > endpos
anchor.chrom <- tbl.target[1, "chr"]

genome.db.uri <- "postgres://whovian/hg38"
project.db.uri <-  "postgres://whovian/wholeBrain"
fp <- FootprintFinder(genome.db.uri, project.db.uri, quiet=TRUE)

shoulder <- 10^4
chrom <- anchor.chrom
start.loc <- anchor.loc - shoulder
end.loc <- anchor.loc  + shoulder
tbl.fp <- getFootprintsInRegion(fp, chrom, start.loc, end.loc)
dim(tbl.fp)
candidate.tfs <- sort(unique(tbl.fp$tf))
length(candidate.tfs)  # 433

print(load(system.file(package="TReNA", "extdata/ampAD.154genes.mef2cTFs.278samples.RData")))  # mtx.xub
mtx.tmp <- mtx.sub - min(mtx.sub) + 0.001
mtx.log2 <- log2(mtx.tmp)
fivenum(mtx.log2)  # [1] -9.9657843  0.8107618  3.6262014  5.4345771 10.0052973

final.candidate.tfs <- intersect(candidate.tfs, rownames(mtx.log2))
trena <- TReNA(mtx.assay=mtx.log2, solver="lasso", quiet=FALSE)

tfs <- setdiff(final.candidate.tfs, target.gene)   # in case target is itself a TF with a footprint in region
tbl <- solve(trena, target.gene, tfs)
print(head(tbl))

