#----------------------------------------------------------------------------------------------------
library(TReNA)
library(plyr)
library(limma)

#----------------------------------------------------------------------------------------------------
assess_methodsAgainstDistributions <- function(mtx.sub, target.gene, tfs)
{    
    #fivenum(mtx.sub)# 0.000000    1.753137   12.346965   43.247467 1027.765854
    
    # Transform with log2 
    mtx.tmp <- mtx.sub - min(mtx.sub) + 0.001
    mtx.log2 <- log2(mtx.tmp)
    #fivenum(mtx.log2)  # [1] -9.9657843  0.8107618  3.6262014  5.4345771 10.005297

    # Transform sub with asinh
    mtx.asinh <- asinh(mtx.sub)
    #fivenum(mtx.asinh)  # [1] 0.000000 1.327453 3.208193 4.460219 7.628290

    # Transform via VOOM transformation
    mtx.voom <- voom(mtx.sub)$E
    #fivenum(mtx.voom) 

    printf("--- Testing LASSO")

    trena <- TReNA(mtx.assay=mtx.sub, solver="lasso", quiet=FALSE)
    lasso1 <- solve(trena, target.gene, tfs)
    lasso1 <- data.frame(gene = rownames(lasso1),
                         lasso.as.is = lasso1$beta)
   
    
    trena <- TReNA(mtx.assay=mtx.log2, solver="lasso", quiet=FALSE)
    lasso2 <- solve(trena, target.gene, tfs)
    lasso2 <- data.frame(lasso.log2 = lasso2$beta,
                         gene = rownames(lasso2))
    
    trena <- TReNA(mtx.assay=mtx.asinh, solver="lasso", quiet=FALSE)
    lasso3 <- solve(trena, target.gene, tfs)
    lasso3 <- data.frame(lasso.asinh = lasso3$beta,
                         gene = rownames(lasso3))

    trena <- TReNA(mtx.assay=mtx.voom, solver="lasso", quiet=FALSE)
    lasso4 <- solve(trena, target.gene, tfs)
    lasso4 <- data.frame(lasso.voom = lasso4$beta,
                         gene = rownames(lasso4))


    lasso1$gene <- as.character(lasso1$gene)
    lasso2$gene <- as.character(lasso2$gene)
    lasso3$gene <- as.character(lasso3$gene)
    lasso4$gene <- as.character(lasso4$gene)

    # Grab the top 10 genes from each
    lasso1.top10 <- lasso1$gene[order(abs(lasso1$lasso.as.is), decreasing=TRUE)][1:10]
    lasso2.top10 <- lasso2$gene[order(abs(lasso2$lasso.log2), decreasing=TRUE)][1:10]
    lasso3.top10 <- lasso3$gene[order(abs(lasso3$lasso.asinh), decreasing=TRUE)][1:10]
    lasso4.top10 <- lasso4$gene[order(abs(lasso4$lasso.voom), decreasing = TRUE)][1:10]

    printf("--- Testing Bayes Spike")

    trena <- TReNA(mtx.assay=mtx.sub, solver="bayesSpike", quiet=FALSE)
    bs1 <- solve(trena, target.gene, tfs)
    bs1 <- data.frame(bs.as.is = bs1$beta,
                         gene = rownames(bs1))

    trena <- TReNA(mtx.assay=mtx.log2, solver="bayesSpike", quiet=FALSE)
    bs2 <- solve(trena, target.gene, tfs)
    bs2 <- data.frame(bs.log2 = bs2$beta,
                         gene = rownames(bs2))

    trena <- TReNA(mtx.assay=mtx.asinh, solver="bayesSpike", quiet=FALSE)
    bs3 <- solve(trena, target.gene, tfs)
    bs3 <- data.frame(bs.asinh = bs3$beta,
                         gene = rownames(bs3))

    trena <- TReNA(mtx.assay=mtx.voom, solver="bayesSpike", quiet=FALSE)
    bs4 <- solve(trena, target.gene, tfs)
    bs4 <- data.frame(bs.voom = bs4$beta,
                         gene = rownames(bs4))

    bs1$gene <- as.character(bs1$gene)
    bs2$gene <- as.character(bs2$gene)
    bs3$gene <- as.character(bs3$gene)
    bs4$gene <- as.character(bs4$gene)
    
    # Grab the top 10 genes from each
    bs1.top10 <- bs1$gene[order(abs(bs1$bs.as.is), decreasing=TRUE)][1:10]
    bs2.top10 <- bs2$gene[order(abs(bs2$bs.log2), decreasing=TRUE)][1:10]
    bs3.top10 <- bs3$gene[order(abs(bs3$bs.asinh), decreasing=TRUE)][1:10]
    bs4.top10 <- bs4$gene[order(abs(bs4$bs.voom), decreasing = TRUE)][1:10]

    printf("--- Testing Random Forest")

    trena <- TReNA(mtx.assay=mtx.sub, solver="randomForest", quiet=FALSE)
    rf.result <- solve(trena, target.gene, tfs)
    rf1 <- data.frame(rf.as.is = rf.result$edges$IncNodePurity,
                      gene = rownames(rf.result$edges))

    trena <- TReNA(mtx.assay=mtx.log2, solver="randomForest", quiet=FALSE)
    rf.result.2 <- solve(trena, target.gene, tfs)
    rf2 <- data.frame(rf.log2 = rf.result.2$edges$IncNodePurity,
                      gene = rownames(rf.result.2$edges))

    trena <- TReNA(mtx.assay=mtx.asinh, solver="randomForest", quiet=FALSE)
    rf.result.3 <- solve(trena, target.gene, tfs)
    rf3 <- data.frame(rf.asinh = rf.result.3$edges$IncNodePurity,
                      gene = rownames(rf.result.3$edges))

    trena <- TReNA(mtx.assay=mtx.voom, solver="randomForest", quiet=FALSE)
    rf.result.4 <- solve(trena, target.gene, tfs)
    rf4 <- data.frame(rf.voom = rf.result.4$edges$IncNodePurity,
                      gene = rownames(rf.result.4$edges))

    rf1$gene <- as.character(rf1$gene)
    rf2$gene <- as.character(rf2$gene)
    rf3$gene <- as.character(rf3$gene)
    rf4$gene <- as.character(rf4$gene)

    # Grab the top 10 genes from each
    rf1.top10 <- rf1$gene[order(abs(rf1$rf.as.is), decreasing=TRUE)][1:10]
    rf2.top10 <- rf2$gene[order(abs(rf2$rf.log2), decreasing=TRUE)][1:10]
    rf3.top10 <- rf3$gene[order(abs(rf3$rf.asinh), decreasing=TRUE)][1:10]    
    rf4.top10 <- rf4$gene[order(abs(rf4$rf.voom), decreasing=TRUE)][1:10]

    # Take the union of all the genes
    all.genes <- unique(c(lasso1.top10,
                       lasso2.top10,
                       lasso3.top10,
                       lasso4.top10,
                       bs1.top10,
                       bs2.top10,
                       bs3.top10,
                       bs4.top10,
                       rf1.top10,
                       rf2.top10,
                       rf3.top10,
                       rf4.top10)) # Returns 41 genes

    # Pull out the specified genes and assemble a table
    lasso1.sub <- subset(lasso1, lasso1$gene %in% all.genes)
    lasso2.sub <- subset(lasso2, lasso2$gene %in% all.genes)
    lasso3.sub <- subset(lasso3, lasso3$gene %in% all.genes)
    lasso4.sub <- subset(lasso4, lasso4$gene %in% all.genes)
    bs1.sub <- subset(bs1, bs1$gene %in% all.genes)
    bs2.sub <- subset(bs2, bs2$gene %in% all.genes)
    bs3.sub <- subset(bs3, bs3$gene %in% all.genes)
    bs4.sub <- subset(bs4, bs4$gene %in% all.genes)
    rf1.sub <- subset(rf1, rf1$gene %in% all.genes)
    rf2.sub <- subset(rf2, rf2$gene %in% all.genes)
    rf3.sub <- subset(rf3, rf3$gene %in% all.genes)
    rf4.sub <- subset(rf4, rf4$gene %in% all.genes)

    # Join it all in a table (requires plyr package)
    tbl.all <- join_all(list(lasso1.sub,
                             lasso2.sub,
                             lasso3.sub,
                             lasso4.sub,
                             bs1.sub,
                             bs2.sub,
                             bs3.sub,
                             bs4.sub,
                             rf1.sub,
                             rf2.sub,
                             rf3.sub,
                             rf4.sub),
                        by = 'gene', type = 'full')

    # Replace NA's in LASSO columns with 0s
    tbl.all[is.na(tbl.all)] <- 0

    # Add the gene correlation column
    tbl.all$gene.cor <- NA
    for(gene in tbl.all$gene){
        tbl.all$gene.cor[[which(tbl.all$gene == gene)]] <-
            rf.result$edges$gene.cor[[which(rownames(rf.result$edges) == gene)]]
    }

    # Order the rows and return it 
    tbl.all <- tbl.all[order(abs(tbl.all$gene.cor), decreasing = TRUE),]
    invisible(tbl.all)

} #assess_methodsAgainstDistributions
#----------------------------------------------------------------------------------------------------
assess_ampAD154AllSolversAndDistributions <- function(){

    # Load the matrix and create the transformed versions
    load(system.file(package="TReNA", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
    target.gene <- "MEF2C"
    tfs <- setdiff(rownames(mtx.sub), "MEF2C")
   
    tbl.all <- assess_methodsAgainstDistributions(mtx.sub,target.gene,tfs)


} #assess_ampAD154AllSolversAndDistributions
#----------------------------------------------------------------------------------------------------
transform_ampADFromRaw <- function(){

    library(edgeR)
    file.path <- "./inst/extdata/MayoRNAseq_RNAseq_TCX_geneCounts2.tsv"
    transposedCounts <- read.table(file.path, header = TRUE, check.names = FALSE, as.is = TRUE, sep = "\t", stringsAsFactors = TRUE)
    browser()
    rownames(transposedCounts) <- transposedCounts$ensembl_id
    
    # Make it into a matrix of counts and drop the ENSEMBL IDs
    mtx <- as.matrix(transposedCounts[, -1])

    # Make a DGEList object
    expr <- DGEList(mtx, group = rep(1, ncol(mtx)))

    # Calculate the normalization factors for TMM
    normFactors <- calcNormFactors(expr, method = ("TMM"))

    # Compute counts per million (returns a matrix)
    normalizedCpm <- cpm(normFactors)
    
    # Organize as a data frame w/ENSEMBL IDs and data values
    ### Where are the IDs to map to ENSEMBL? 
    ids <- as.data.frame(row.names(normalizedCpm))
    colnames(ids) <-c("ensembl_id")
    row.names(normalizedCpm) <- NULL
    df <- as.data.frame(normalizedCpm)
    newDF <- cbind(ids, df)

    # Return the normalized data frame
    invisible(newDF)
    
} #transform_ampADFromRaw
#----------------------------------------------------------------------------------------------------
