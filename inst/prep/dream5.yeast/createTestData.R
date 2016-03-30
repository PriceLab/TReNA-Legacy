library(org.Sc.sgd.db)
library(GO.db)
library(RUnit)

mtx <- as.matrix(read.table("net4_expression_data.tsv", sep="\t",  header=TRUE, as.is=TRUE, row.names=NULL))  # 536 x 5950
tbl.ids <- read.table("net4_gene_ids.tsv", sep="\t", header=FALSE, as.is=TRUE)  # 5950 x 2
colnames(tbl.ids) <- c("G", "orf")
tfs <- read.table("net4_transcription_factors.tsv", sep="\t", header=FALSE, as.is=TRUE)$V1

  # leave mtx names in Gn style for now
  #indices <- match(colnames(mtx), tbl.ids$V1)
  #colnames(mtx) <- tbl.ids$V2[indices]
  #tbl.geneNames <- select(org.Sc.sgd.db, keys=colnames(mtx), keytype="ORF", columns <- c("ORF", "GENENAME"))
  #length(which(is.na(tbl.geneNames$GENENAME)))
  #new.names <- tbl.geneNames$GENENAME
  #nas <- which(is.na(new.names))
  #new.names[nas] <- tbl.geneNames$ORF[nas]
  #colnames(mtx) <- new.names

# prepare a lookup table Gn -> gene symbol
db <- org.Sc.sgd.db
keytypes(db)  # learn the columns
tbl.orfSymbol <- select(db, keys=tbl.ids$orf, keytype="ORF", columns <- c("ORF", "GENENAME"))
tbl.merged <- merge(tbl.ids, tbl.orfSymbol, by.x="orf", by.y="ORF")
tbl.ids <- tbl.merged[, c("G", "orf", "GENENAME")]
colnames(tbl.ids) <- c("G", "orf", "symbol")
#deleters <- grep("decoy", tbl.ids$orf)
#if(length(deleters) > 0)
#   tbl.ids <- tbl.ids[-deleters,]
NAs <- which(is.na(tbl.ids$symbol))
if(length(NAs) > 0)
    tbl.ids$symbol[NAs] <- tbl.ids$orf[NAs]

# now replace mtx colnames
new.colnames <- tbl.ids[match(colnames(mtx), tbl.ids$G), "symbol"]
colnames(mtx) <- new.colnames
new.rownames <- sprintf("sample.%03d", 1:nrow(mtx))
rownames(mtx) <- new.rownames
mtx <- t(mtx)

dim(mtx)# [1] 5667  536
tbl.gold <- read.table("yeastGoldStandard.tsv", sep="\t", as.is=TRUE)
colnames(tbl.gold) <- c("TF", "target", "score")
table(tbl.gold$score) #      0: 223262   1: 3940
tbl.gold <- subset(tbl.gold, score==1)
all(tbl.gold$TF %in% tbl.ids$G) # [1] TRUE
all(tbl.gold$target %in% tbl.ids$G) # [1] TRUE
tbl.gold$TF <- tbl.ids$symbol[match(tbl.gold$TF, tbl.ids$G)]
tbl.gold$target <- tbl.ids$symbol[match(tbl.gold$target, tbl.ids$G)]

# spot check the top tfs, make sure they have proper GO terms
top.tfs <- names(tail(sort(table(tbl.gold$TF))))
tbl.go <- select(db, keys=top.tfs, keytype="GENENAME", columns="GO")
tbl.gobp <- subset(tbl.go, ONTOLOGY=="BP")
tbl.goTerms <- select(GO.db, keys=tbl.gobp$GO, keytype="GOID", columns=c("GOID", "TERM"))
tbl.gobp <- merge(tbl.gobp, tbl.goTerms, by.x="GO", by.y="GOID")
checkEquals(as.list(tail(sort(table(tbl.gobp$TERM)), n=2)),
               list(`regulation of transcription, DNA-templated`=36,
                    `transcription, DNA-templated` = 64))

tfs <- read.table("net4_transcription_factors.tsv", sep="\t", as.is=TRUE)$V1  # 333

# add a correlations column to tbl.gold:  it's a weak signal, approximatly normal
# fivenum(correlations)  # [1] -0.626543305 -0.003879309  0.167823346  0.333898289  0.914131294

correlations <- vector("numeric", nrow(tbl.gold))
correlations <- sapply(1:nrow(tbl.gold), function(row) {
    tf <- tbl.gold$TF[row];
    target <- tbl.gold$target[row];
    cor(mtx[tf,], mtx[target,])
    })

tbl.gold$cor <- correlations

serializedDataFile <- "../../extdata/dream5.net4.yeast.RData"
printf("writing mtx (%d,%d) tbl.gold (%d,%d) and tbl.ids(%d,%d) to %s",
       nrow(mtx), ncol(mtx), nrow(tbl.gold), ncol(tbl.gold), nrow(tbl.ids), ncol(tbl.ids),
       serializedDataFile)
save(mtx.expr=mtx, gold=tbl.gold, id.mapper=tbl.ids, file=serializedDataFile)
