library(TReNA)

db <- dbConnect(PostgreSQL(), host="whovian", user="trena", password="trena", dbname="hg38")

# find the chrom loc for piez02

target.gene <- "PIEZO2"
biotype <- "protein_coding"
moleculetype <- "gene"
query <- paste("select gene_name, chr, start, endpos, strand from gtf where ",
                 sprintf("gene_name='%s' ", target.gene),
                 sprintf("and gene_biotype='%s' and moleculetype='%s'", biotype, moleculetype),
                 collapse=" ")

tbl.target <- dbGetQuery(db, query)
#  gene_name   chr    start   endpos strand
#     PIEZO2 chr18 10666483 11148762      -

anchor.loc <- tbl.target[1, "start"]   # note: on reverse strand, so start > endpos
anchor.chrom <- tbl.target[1, "chr"]

# find all genes withing 1MB of this anchor.loc

shoulder <- 10^6
query <- paste("select gene_name, chr, start, endpos, strand, gene_biotype, moleculetype from gtf ",
               sprintf("where start >= %d ", anchor.loc - shoulder),
               sprintf("and chr = '%s' ", anchor.chrom),
               sprintf("and endpos <=  %d ", anchor.loc + shoulder),
               sprintf("and gene_biotype in ('protein_coding', 'lincRNA') "),
               sprintf("and moleculetype = 'gene'"),
               collapse=" ")

tbl.region <- dbGetQuery(db, query)   # 6 protein-coding, 14 lincRNA genes for shoulder of 10^6

genome.db.uri <- "postgres://whovian/hg38"
project.db.uri <-  "postgres://whovian/wholeBrain"
fp <- FootprintFinder(genome.db.uri, project.db.uri, quiet=TRUE)

shoulder <- 10^3
# create an empty list, then fill it with data.frames, one for each gene  in the specified
# region around PIEZO2

tbls.all <- list()

for(row in 1:nrow(tbl.region)){
   chrom <- tbl.region$chr[row]
   gene <- tbl.region$gene_name[row]
   printf("about to query fp for region of gene %s", gene)
   start.loc <- tbl.region$start[row] - shoulder
   end.loc   <- tbl.region$start[row] + shoulder
   tbl.fp <- getFootprintsInRegion(fp, chrom, start.loc, end.loc)
   if(nrow(tbl.fp) > 0)
      tbl.fp$target <- gene
   tbls.all[[gene]] <- tbl.fp
   } # for gene

# collapse the list into a single data.frame
tbl <- do.call('rbind', tbls.all)

# ponder the distribution of footprints in the vicinity of each gene
print(table(tbl$target))

# create a list of lists, one per gene in the nominated region, in
# which each sublist is candidate transcription factors for that gene

target.genes <- unique(tbl$target)
tfs.per.target <- lapply(target.genes, function(target.gene) {tfs <- unique(subset(tbl, target==target.gene)$tf)})
names(tfs.per.target) <- target.genes
print(sapply(tfs.per.target, length))

