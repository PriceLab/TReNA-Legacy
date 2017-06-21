library(TReNA)
library(RUnit)
#----------------------------------------------------------------------------------------------------
# the vrk2 promoter snp
# chr2:57907313-57907333
sequence <- "ACCAGCATGCAAATTAGACAA"
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_basicConstructor()

} # runTests
#------------------------------------------------------------------------------------------------------------------------
test_basicConstructor <- function(reuse=FALSE)
{
   viz <- TReNA.Viz()
   viz

} # test_basicConstructor
#----------------------------------------------------------------------------------------------------
test_addDemoGraph <- function()
{
   printf("--- test_addDemoGraph")
   g <- graphNEL(edgemode="directed")
   g <- addNode(c("A", "B"), g)
   g <- addEdge("A", "B", g)

   tv <- TReNA.Viz(portRange=12011)
   TReNA:::addGraph(tv, g)
   layout(tv, "grid")

   tv

} # test_addDemoGraph
#----------------------------------------------------------------------------------------------------
test_addBedTrack <- function()
{
   printf("--- test_addBedTrack")

   tv <- TReNA.Viz(portRange=11011:11051)

   tbl.bed <- data.frame(chrom=rep("chr18", 3),
                         start=c(26860880, 26860974, 26861152),
                         end=c(26860887, 26860983, 26861161),
                         name=c("MA0808.1", "MA0090.2", "MA0809.1"),
                         score=c(0.8234864, 0.8221845, 0.8307060),
                         stringsAsFactors=FALSE)
    addBedTrackFromDataFrame(tv, sprintf("test_%d", round(sample(1:100, 1))), tbl.bed)
    addBedTrackFromHostedFile(tv,
                              trackName="brain HINT",
                              uri="http://pshannon.systemsbiology.net/annotations/brain_hint.bed.gz",
                              index.uri="http://pshannon.systemsbiology.net/annotations/brain_hint.bed.gz.tbi",
                              displayMode="SQUISHED")
    addBedTrackFromHostedFile(tv,
                             trackName="EncodeDHSclustered",
                             uri="http://pshannon.systemsbiology.net/annotations/dhsClusters_hg38.bed.gz",
                             index.uri="http://pshannon.systemsbiology.net/annotations/dhsClusters_hg38.bed.gz.tbi",
                             displayMode="SQUISHED")
    addBedTrackFromHostedFile(tv,
                             trackName="AQP4 snps",
                             uri="http://pshannon.systemsbiology.net/annotations/aqp4-chr18-snps.bed",
                             index.uri=NA,
                             displayMode="SQUISHED")


    tv

} # test_addBedTrack
#----------------------------------------------------------------------------------------------------
test_.graphToJSON <- function()
{
   printf("--- test_.graphToJSON")

   g <- graphNEL(edgemode='directed')

   nodeDataDefaults(g, attr='label') <- 'default node label'
   nodeDataDefaults(g, attr='type') <- 'default node label'
   nodeDataDefaults(g, attr='count') <- 0
   nodeDataDefaults(g, attr='score') <- 0.0

   edgeDataDefaults(g, attr='edgeType') <- 'undefined'
   edgeDataDefaults(g, attr='count') <- 0
   edgeDataDefaults(g, attr='score') <- 0.0

   g <- graph::addNode('A', g)
   g <- graph::addNode('B', g)
   g <- graph::addNode('C', g)

   all.nodes <- nodes(g)
   nodeData(g, c('A', 'B', 'C'), 'type') <- c('t_one', 't_two', 't_three')
   nodeData(g, all.nodes, 'label') <- all.nodes
   nodeData(g, all.nodes, 'score') <- runif(3, 0, 1)
   nodeData(g, all.nodes, 'count') <- sample(1:10, 3)

     # now add 1 edge
   g <- graph::addEdge('A', 'B', g)
   g <- graph::addEdge('B', 'C', g)
   g <- graph::addEdge('C', 'A', g)

   edgeData(g, 'A', 'B', 'edgeType') <- 'et_one'
   edgeData(g, 'B', 'C', 'edgeType') <- 'et_two'

   edgeData(g, 'A', 'B', 'score') <- runif(1, 0, 1)
   edgeData(g, 'C', 'A', 'score') <- runif(1, 0, 1)

   edgeData(g, 'B', 'C', 'count') <- sample(1:10, 1)
   edgeData(g, 'C', 'A', 'count') <- sample(1:10, 1)

   g.json <- TReNA:::.graphToJSON(g)
      # a pretty good check, but only an assist to the real test, which is to load this the
      # browser with cy.json(JSON.parse(<g.json string>))

   tbl <- fromJSON(g.json)[[1]][[1]]
   checkEquals(nrow(tbl), length(nodes(g)) + length(edgeNames(g)))
   checkEquals(colnames(tbl), c("id", "label", "type", "count", "score", "source", "target", "edgeType"))
   checkEquals(unlist(lapply(tbl, class), use.names=FALSE),
               c("character", "character", "character", "integer", "numeric", "character", "character", "character"))
   checkEquals(tbl$type[1:3], c("t_one", "t_two", "t_three"))
   checkEquals(tbl$edgeType[4:6], c("et_one", "et_two", "undefined"))


} # test_.graphToJSON
#----------------------------------------------------------------------------------------------------
test_geneRegulatoryModelToGraph <- function()
{
   printf("--- test_geneRegulatoryModelToGraph")

   load(system.file(package="TReNA", "extdata", "aqp4-model-regRegions.RData"))
   stopifnot(exists("tbl.gm"))
   rows <- nrow(tbl.gm)
   tbl.grm <- data.frame(tf=rownames(tbl.gm),
                         pearson=tbl.gm$gene.cor,
                         spearman=rep(NA_real_, rows),
                         betaLasso=rep(NA_real_, rows),
                         randomForest=tbl.gm$IncNodePurity,
                         pcaMax=rep(NA_real_, rows),
                         concordance=rep(NA_real_, rows),
                         stringsAsFactors=FALSE)
   stopifnot(exists("tbl.reg"))
   keepers <- unique(unlist(lapply(tbl.grm$tf, function(tf) grep(tf, tbl.reg$tf))))
   tbl.reg <- tbl.reg[keepers,]
   target.gene <- "AQP4"
   aqp4.tss <- 26865884   # minus strand, leftward positions are negative (downstream), rightward are positive (upstream)
   tbl.reg$distance.from.tss <- tbl.reg$motifStart - aqp4.tss

   tfs.expanded <- strsplit(tbl.reg$tf, ";", fixed=TRUE)
   tfs.winnowed <- vector(mode="character", length=nrow(tbl.reg))

   for(i in 1:length(tfs.expanded)){
      full.set <- tfs.expanded[[i]]
      reduced.set <- intersect(full.set, tbl.grm$tf)
      if(all(is.na(reduced.set))){
         #printf("no match to %s", tbl.reg$tf[i])
         reduced.set <- "NO_INTERSECTION"
         }
      else{
        if(length(reduced.set) > 1){
           #printf("multiple intersection: %s", tbl.reg$tf[i])
           reduced.set <- paste(reduced.set, collapse=";")
           } # if > 1
        } # not na
      tfs.winnowed[[i]] <- reduced.set
      } # for i

   tbl.reg$tf <- tfs.winnowed
   failures <- which(tfs.winnowed == "NO_INTERSECTION")
   if(length(failures) > 0)
      tbl.reg <- tbl.reg[-failures,]

   set.seed(17)
   #sample.edges <- unique(unlist(lapply(tbl.grm$tf, function(tf) head(grep(tf, tbl.reg$tf)))))[1:100]
   #tbl.reg <- tbl.reg[1:100,]
   tv <- TReNA.Viz(portRange=11011:11051)

   g <- geneRegulatoryModelToGraph(tv, target.gene, tbl.grm, tbl.reg)
   checkEquals(sort(nodes(g)),
                   c("COL1A1","COL1A1.fp.upstream.01505.L11.MA0471.1","COL1A1.fp.upstream.01523.L10.MA0599.1",
                     "COL1A1.fp.upstream.01523.L11.MA0079.3","COL1A1.fp.upstream.01523.L15.MA0516.1",
                     "COL1A1.fp.upstream.01530.L12.GLI1..3.p2","COL1A1.fp.upstream.01532.L9.ZIC1..3.p2",
                     "E2F7","GLI1","GLI2","GLI3","KLF12","KLF14","KLF2","KLF3","SP1","ZIC1"))

   checkEquals(sort(edgeNames(g)),
               c("COL1A1.fp.upstream.01505.L11.MA0471.1~COL1A1",
                 "COL1A1.fp.upstream.01523.L10.MA0599.1~COL1A1",
                 "COL1A1.fp.upstream.01523.L11.MA0079.3~COL1A1",
                 "COL1A1.fp.upstream.01523.L15.MA0516.1~COL1A1",
                 "COL1A1.fp.upstream.01530.L12.GLI1..3.p2~COL1A1",
                 "COL1A1.fp.upstream.01532.L9.ZIC1..3.p2~COL1A1",
                 "E2F7~COL1A1.fp.upstream.01505.L11.MA0471.1",
                 "GLI1~COL1A1.fp.upstream.01530.L12.GLI1..3.p2",
                 "GLI2~COL1A1.fp.upstream.01530.L12.GLI1..3.p2",
                 "GLI3~COL1A1.fp.upstream.01530.L12.GLI1..3.p2",
                 "KLF12~COL1A1.fp.upstream.01523.L10.MA0599.1",
                 "KLF12~COL1A1.fp.upstream.01523.L11.MA0079.3",
                 "KLF12~COL1A1.fp.upstream.01523.L15.MA0516.1",
                 "KLF14~COL1A1.fp.upstream.01523.L10.MA0599.1",
                 "KLF14~COL1A1.fp.upstream.01523.L11.MA0079.3",
                 "KLF14~COL1A1.fp.upstream.01523.L15.MA0516.1",
                 "KLF2~COL1A1.fp.upstream.01523.L10.MA0599.1",
                 "KLF2~COL1A1.fp.upstream.01523.L11.MA0079.3",
                 "KLF2~COL1A1.fp.upstream.01523.L15.MA0516.1",
                 "KLF3~COL1A1.fp.upstream.01523.L10.MA0599.1",
                 "KLF3~COL1A1.fp.upstream.01523.L11.MA0079.3",
                 "KLF3~COL1A1.fp.upstream.01523.L15.MA0516.1",
                 "SP1~COL1A1.fp.upstream.01523.L10.MA0599.1",
                 "SP1~COL1A1.fp.upstream.01523.L11.MA0079.3",
                 "SP1~COL1A1.fp.upstream.01523.L15.MA0516.1",
                 "ZIC1~COL1A1.fp.upstream.01532.L9.ZIC1..3.p2"))

       # --- select one footprint node, check its attributes
    checkEquals(nodeData(g, "COL1A1.fp.upstream.01523.L15.MA0516.1", attr="type")[[1]], "footprint")
    checkEquals(nodeData(g, "COL1A1.fp.upstream.01523.L15.MA0516.1", attr="label")[[1]], "MA0516.1")
    checkEquals(nodeData(g, "COL1A1.fp.upstream.01523.L15.MA0516.1", attr="motif")[[1]], "MA0516.1")
    print(checkEquals(nodeData(g, "COL1A1.fp.upstream.01523.L15.MA0516.1", attr="distance")[[1]], 1523))

      # null attributes, which should never happen (but just did, by assigning
      # the noa "pearson" from tbl.gm$perason (note mispelling)) are not
      # return by the nodeData function.  detect them, therefore, by counting lengths

    nodeCount <- length(nodes(g))
    checkEquals(length(unlist(nodeData(g, attr="type"), use.names=FALSE)), nodeCount)
    checkEquals(length(unlist(nodeData(g, attr="label"), use.names=FALSE)), nodeCount)
    checkEquals(length(unlist(nodeData(g, attr="distance"), use.names=FALSE)), nodeCount)
    checkEquals(length(unlist(nodeData(g, attr="pearson"), use.names=FALSE)), nodeCount)
    checkEquals(length(unlist(nodeData(g, attr="randomForest"), use.names=FALSE)), nodeCount)
    checkEquals(length(unlist(nodeData(g, attr="betaLasso"), use.names=FALSE)), nodeCount)
    checkEquals(length(unlist(nodeData(g, attr="motif"), use.names=FALSE)), nodeCount)
    checkEquals(length(unlist(nodeData(g, attr="xPos"), use.names=FALSE)), nodeCount)
    checkEquals(length(unlist(nodeData(g, attr="yPos"), use.names=FALSE)), nodeCount)

    edgeCount <- length(edgeNames(g))
    checkEquals(length(unlist(edgeData(g, attr="edgeType"), use.names=FALSE)), edgeCount)

    invisible(g)

} # test_geneModelToGraph
#------------------------------------------------------------------------------------------------------------------------
validModelData <- function(models)
{
   required.geneModelColumnNames <- c("tf", "pearson", "spearman", "betaLasso", "randomForest", "pcaMax", "concordance")
   required.regulatoryRegionsColumnNames <- c("motifName", "chrom", "motifStart", "motifEnd", "strand",
                                              "motifScore", "motifRelativeScore", "match",
                                              "distance.from.tss",
                                              "chromStart", "chromEnd", "seq", "status", "tf", "regionName")
   valid <- TRUE;  # be optimistic
   model.names <- names(models)
   checkTrue(length(model.names) > 0)
      # for simplicity, and especially for use in graph node attributes, we want model names with NO spaces
   noSpacesInNames <- all(grepl(" ", model.names) == FALSE)
   if(!noSpacesInNames){
      valid <- FALSE
      warning("model names should contain no spaces")
      }

   for(model.name in names(models)){
      stopifnot(all(names(models[[model.name]]) %in% c("tbl.model", "tbl.regulatoryRegions")))
      tbl.model <- models[[model.name]]$tbl.model

      missing.in.model <- setdiff(required.geneModelColumnNames, colnames(tbl.model))
      valid.0 <- length(missing.in.model) == 0;
      if(!valid.0)
         warning(sprintf("missing columns in tbl.model for %s: %s", model.name, paste(missing.in.model, collapse=", ")))
      valid <- valid & valid.0
      missing.in.regions <- setdiff(required.geneModelColumnNames, colnames(tbl.model))
      valid.1 <- length(missing.in.regions) == 0;
      if(!valid.1)
         warning(sprintf("missing columns in tbl.model for %s: %s", model.name, paste(missing.in.regions, collapse=", ")))
      valid <- valid & valid.1
      } # for model.name

    valid

} # validModelData
#------------------------------------------------------------------------------------------------------------------------
test_validModelData <- function()
{
   printf("--- test_validModelData")

   tbl.model.test <- data.frame(tf=c("a"),
                                pearson=c(0),
                                spearman=c(0),
                                betaLasso=c(0),
                                randomForest=c(0),
                                pcaMax=c(0),
                                concordance=c(0),
                                stringsAsFactors=FALSE)

   tbl.regions.test <- data.frame(motifName=c("a"),
                                  chrom=c("a"),
                                  motifStart=c(0),
                                  motifEnd=c(0),
                                  strand=c("+"),
                                  motifScore=c(0),
                                  motifRelativeScore=c(0),
                                  match=c("a"),
                                  distance.from.tss=c(0),
                                  chromStart=c(0),
                                  chromEnd=c(0),
                                  seq=c(0),
                                  status=c("wt"),
                                  tf=c("a"),
                                  regionName=c("xyz"),
                                  stringsAsFactors=FALSE
                                  )

   models <- list(wt=list(tbl.model=tbl.model.test,  tbl.regulatoryRegions=tbl.regions.test),
                  mut=list(tbl.model=tbl.model.test, tbl.regulatoryRegions=tbl.regions.test),
                  m3=list(tbl.model=tbl.model.test,  tbl.regulatoryRegions=tbl.regions.test))

   checkTrue(validModelData(models))

   names(models)[1] <- "wild type"
   suppressWarnings(checkTrue(!validModelData(models)))

} # test_validModelData
#------------------------------------------------------------------------------------------------------------------------
buildMultiModelGraph <- function(target.gene, models)
{
   stopifnot(validModelData(models))

   g <- graphNEL(edgemode = "directed")
   model.names <- names(models)

   node.attribute.specs <- list(type="undefined",
                                label="default node label",
                                distance=0,
                                pearson=0,
                                randomForest=0,
                                pcaMax=0,
                                concordance=0,
                                betaLasso=0,
                                motif="",
                                xPos=0,
                                yPos=0)
   edge.attribute.spec <- list(edgeType="undefined")
   attribute.classes <- c("", model.names)  # "" (no prefix) is the currently displayed set of attibutes

      # create current version of these attributes, and then
      # per-model versions, which get mapped to current
      # in response to user's interactive choice on the cyjs user interface
      # the "current version" is, e.g., "distance".
      # per-model ("wt" and "mut" versions) become "wt.distance" and "mut.distance"
      # and are used by copying e.g. all wt.xxx attributes into the current (non-prefixed)
      # attribute, upon which the cyjs style is defined

   for(class.name in attribute.classes){
      class.name.prefix <- class.name  # with possible "." appended, permits standard and model-specific attributes
      if(nchar(class.name) > 0)
         class.name.prefix <- sprintf("%s.", class.name)
      noa.names.without.prefix <- names(node.attribute.specs)
      noa.names <- sprintf("%s%s", class.name.prefix, noa.names.without.prefix)
      noa.count <- length(node.attribute.specs)
      for(i in 1:noa.count){
         nodeDataDefaults(g, attr=noa.names[i]) <- node.attribute.specs[[noa.names.without.prefix[i]]]
         }
      } # for class

   edgeDataDefaults(g, attr = "edgeType") <- "undefined"

   tfs <- c()
   regulatoryRegions <- c()

   for(model in models){  # collect all the tf and regulatory region nodes
     tbl.model <- model$tbl.model
     tfs <- unique(c(tfs, tbl.model$tf))
     tbl.reg <- model$tbl.regulatoryRegions
     regulatoryRegions <- unique(c(regulatoryRegions, tbl.reg$regionName))
     } # for model

   all.nodes <- unique(c(target.gene, tfs, regulatoryRegions))
   g <- addNode(all.nodes, g)

   nodeData(g, target.gene, "type") <- "targetGene"
   nodeData(g, tfs, "type")         <- "TF"
   nodeData(g, regulatoryRegions, "type")  <- "regulatoryRegion"
   nodeData(g, all.nodes, "label")  <- all.nodes

      # add edges, edge attribute, and the constant attributes for all of the regulatoryRegion nodes

   for(model in models){
      tfs <- model$tbl.regulatoryRegions$tf
      regRegions <- model$tbl.regulatoryRegions$regionName
      suppressWarnings(g <- addEdge(tfs, regRegions, g))
      edgeData(g,  tfs, regRegions, "edgeType") <- "bindsTo"
      suppressWarnings(g <- addEdge(regRegions, target.gene, g))
      edgeData(g, regRegions, target.gene, "edgeType") <- "regulatorySiteFor"
      nodeData(g, tbl.reg$regionName, "label") <- tbl.reg$motifName
      nodeData(g, tbl.reg$regionName, "distance") <- tbl.reg$distance.from.tss
      nodeData(g, tbl.reg$regionName, "motif") <- tbl.reg$motifName
      } # for model

      # now copy in the first model's tf node data

   tbl.model <- models[[1]]$tbl.model
   nodeData(g, tbl.model$tf, attr="randomForest") <- tbl.model$randomForest
   nodeData(g, tbl.model$tf, attr="pearson") <- tbl.model$pearson

     # now copy in each of the model's tf node data in turn
   model.names <- names(models)
   for(model.name in model.names){
      tbl.model <- models[[model.name]]$tbl.model
      noa.name <- sprintf("%s.%s", model.name, "randomForest")
      nodeData(g,  tbl.model$tf, attr=noa.name) <- tbl.model$randomForest
      noa.name <- sprintf("%s.%s", model.name, "pearson")
      nodeData(g,  tbl.model$tf, attr=noa.name) <- tbl.model$pearson
      } # for model.name

   xyz <- 99

#    browser();
#    xyz <- 99
#
#    nodeData(g, tfs, "pearson") <- tbl.gm$pearson
#    nodeData(g, tfs, "betaLasso") <- tbl.gm$betaLasso
#    nodeData(g, tfs, "randomForest") <- tbl.gm$randomForest
#    nodeData(g, tfs, "pcaMax") <- tbl.gm$pcaMax
#    nodeData(g, tfs, "concordance") <- tbl.gm$concordance
#
#    #browser()
#    g <- addEdge(tbl.reg$tf, tbl.reg$regionName, g)
#    edgeData(g,  tbl.reg$tf, tbl.reg$regionName, "edgeType") <- "bindsTo"
#
#    g <- graph::addEdge(tbl.reg$regionName, target.gene, g)
#    edgeData(g, tbl.reg$regionName, target.gene, "edgeType") <- "regulatorySiteFor"
#
    g

} # buildMultiModelGraph
#------------------------------------------------------------------------------------------------------------------------
test_multiple.geneRegulatoryModelsToGraph <- function()
{
   printf("--- test_multiple.geneRegulatoryModelToGraph")
   load(system.file(package="TReNA", "extdata", "twoAQP4modelsForTesting.RData"))
   models <- list(wt=list(tbl.model=x.wt$tbl.model, tbl.regulatoryRegions=x.wt$tbl.regulatoryRegions),
                  rs3875089=list(tbl.model=x.mut$tbl.model, tbl.regulatoryRegions=x.mut$tbl.regulatoryRegions))
   g <- buildMultiModelGraph("AQP4", models)

   noa.names <- sort(names(nodeDataDefaults(g)))
   checkEquals(length(noa.names), 33)
   checkEquals(length(grep("rs3875089.", noa.names, fixed=TRUE)), 11)
   checkEquals(length(grep("wt.", noa.names, fixed=TRUE)), 11)
   g.lo <- TReNA:::addGeneModelLayout(g)

   list(graph=g, modelNames=names(models))

} # test_multiple.geneRegulatoryModelToGraph
#------------------------------------------------------------------------------------------------------------------------
test.displayJSON <- function()
{
   # optional
   # try first with a known good example from the RCyjs test suite
   # standard.json.test.file <- system.file(package="RCyjs", "extdata", "g.json")
   # checkTrue(file.exists(standard.json.test.file))

   g <- test.tablesToFullGraph()
   g.lo <- addGeneModelLayout(g)

   g.json <- graphToJSON(g.lo)
   checkTrue(nchar(g.json) > 7000)
   g.json <- sprintf("network = %s", g.json)
   temp.filename <- tempfile(fileext=".json")
   write(g.json, file=temp.filename)

   PORTS=9047:9097
   rcy <- RCyjs(PORTS, graph=graphNEL())
   setBackgroundColor(rcy, "#FAFAFA")
   setDefaultEdgeColor(rcy, "blue")
   redraw(rcy)

   httpAddJsonGraphFromFile(rcy, temp.filename)
   fit(rcy)
   httpSetStyle(rcy, "style.js")
   #layout(rcy, "grid")
   checkEquals(nrow(getNodes(rcy)), length(nodes(g)))
   rcy

} # test.displayJSON
#------------------------------------------------------------------------------------------------------------------------
test.addGeneModelLayout <- function()
{
   printf("--- test.addGeneModelLayout")

   g <- test.tablesToFullGraph()
   checkEquals(as.integer(nodeData(g, attr="xPos")), rep(0, length(nodes(g))))
   checkEquals(as.integer(nodeData(g, attr="yPos")), rep(0, length(nodes(g))))

   g.lo <- addGeneModelLayout(g)
   checkEquals(as.integer(nodeData(g.lo, attr="xPos")),
               c(0, 1530, 1523, 1523, 1523, 1530, 1530, 1523, 1523, 1505, 1532, 1505, 1523, 1523, 1523, 1530, 1532))
   yPos <- as.integer(nodeData(g.lo, attr="yPos"))
     # some random placement employed for the TFs.  all footprints/motifs are at zero.  target.gene is at -200
   checkTrue(all(yPos) >= -200)
   checkTrue(all(yPos) <= 1200)
     # sample values:       c(-200, 698, 839, 1130, 651, 868, 1087, 604, 564, 541, 597, 0, 0, 0, 0, 0, 0))

} # test.addGeneModelLayout
#------------------------------------------------------------------------------------------------------------------------
test.geneModelLayoutNaNBug <- function()
{
   region <- 'chr17:50,201,552-50,201,727'
   target.gene <- "COL1A1"
   result <- createGeneModel(mtx=mtx.gtexPrimary, target.gene, region)
   tbl.model <- result$model
   tbl.reg   <- result$regulatoryRegions
   g <- tablesToFullGraph(tbl.model, tbl.reg)
   g.lo <- addGeneModelLayout(g)
   xPos <- unlist(nodeData(g.lo, attr="xPos"), use.names=FALSE)
   yPos <- unlist(nodeData(g.lo, attr="yPos"), use.names=FALSE)
   checkTrue(!any(is.nan(xPos)))
   checkTrue(!any(is.nan(yPos)))

} # test.geneModelLayoutNaNBug
#------------------------------------------------------------------------------------------------------------------------
test.g <- function(x)
{
   g <- x$graph
   modelNames <- x$modelNames
   tv <- TReNA.Viz(portRange=11011:11051)
   g.lo <- TReNA:::addGeneModelLayout(g)
   printf("addGraph")
   addGraph(tv, g.lo, modelNames)
   loadStyle(tv, "style.js")
   fit(tv)

   tv

} # test.g
#------------------------------------------------------------------------------------------------------------------------
