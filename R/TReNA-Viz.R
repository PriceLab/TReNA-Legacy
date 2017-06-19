#----------------------------------------------------------------------------------------------------
cyjsBrowserFile <- system.file(package="TReNA", "scripts", "viz", "trenaViz.html")
printf <- function(...) print(noquote(sprintf(...)))
#----------------------------------------------------------------------------------------------------
.TReNA.Viz <- setClass ("TReNA.Viz",
                        representation = representation(model="data.frame"),
                        contains = "BrowserVizClass",
                        prototype = prototype (uri="http://localhost", 9000)
                        )


#----------------------------------------------------------------------------------------------------
setGeneric('addGraph',         signature='obj', function(obj, graph, modelNames=list()) standardGeneric ('addGraph'))
setGeneric('httpAddGraph',     signature='obj', function(obj, graph) standardGeneric ('httpAddGraph'))
setGeneric('loadStyle',        signature='obj', function(obj, filename) standardGeneric ('loadStyle'))
setGeneric('fit',              signature='obj', function(obj, padding=30) standardGeneric('fit'))
setGeneric('fitSelected',      signature='obj', function(obj, padding=30) standardGeneric('fitSelectedContent'))
setGeneric('selectNodes',      signature='obj', function(obj, nodeIDs) standardGeneric('selectNodes'))
setGeneric('getSelectedNodes', signature='obj', function(obj) standardGeneric('getSelectedNodes'))
setGeneric('clearSelection',   signature='obj', function(obj) standardGeneric('clearSelection'))
setGeneric('sfn',              signature='obj', function(obj) standardGeneric('sfn'))
setGeneric('addBedTrackFromDataFrame',   signature='obj', function(obj, trackName, tbl.bed, displayMode="COLLAPSED", color)
                                  standardGeneric('addBedTrackFromDataFrame'))
setGeneric('addBedTrackFromHostedFile',   signature='obj',
                        function(obj, trackName, uri, index.uri=NA, displayMode="COLLAPSED", color)
                                  standardGeneric('addBedTrackFromHostedFile'))
setGeneric('showGenomicRegion',    signature='obj', function(obj, regionString) standardGeneric('showGenomicRegion'))
setGeneric('getGenomicRegion',     signature='obj', function(obj, regionString) standardGeneric('getGenomicRegion'))
setGeneric('layout',              signature='obj', function(obj, strategy) standardGeneric('layout'))
setGeneric('layoutStrategies',    signature='obj', function(obj) standardGeneric('layoutStrategies'))
setGeneric('geneRegulatoryModelToGraph',    signature='obj', function(obj, target.gene, tbl.gm, tbl.reg)
                standardGeneric('geneRegulatoryModelToGraph'))
#----------------------------------------------------------------------------------------------------
# constructor
TReNA.Viz = function(portRange=11000:11025, host="localhost", title="TReNA-Viz", quiet=TRUE)
{

   model <- data.frame()
   obj <- .TReNA.Viz(BrowserViz(portRange=portRange, host=host, title=title,
                                quiet=quiet, browserFile=cyjsBrowserFile,
                                httpQueryProcessingFunction=myQP),
                     model=model)

  while (!browserResponseReady(obj)){
      Sys.sleep(.1)
      }
   if(!quiet) {
      message(sprintf("BrowserViz ctor called from TReNA-Viz ctor got browser response"))
      print(getBrowserResponse(obj))
      }

   obj

} # TReNA.Viz constructor
#----------------------------------------------------------------------------------------------------
setMethod('addGraph', 'TReNA.Viz',

  function (obj, graph, modelNames=list()) {
     printf("TReNA.Viz::addGraph");
     print(graph)
     printf("--- calling .graphToJSON");
     g.json <- .graphToJSON(graph)
     payload <- list(graph=g.json, modelNames=modelNames)
     printf("about to send g.json: %d chars", nchar(g.json));
     send(obj, list(cmd="addGraph", callback="handleResponse", status="request", payload=payload))
     while (!browserResponseReady(obj)){
        Sys.sleep(.1)
        }
     printf("browserResponseReady")
     getBrowserResponse(obj);
     })

#----------------------------------------------------------------------------------------------------
setMethod('httpAddGraph', 'TReNA.Viz',

  function (obj, graph) {
     printf("TReNA.Viz::httpAddGraph");
     print(graph)
     g.json <- paste("network = ", as.character(biocGraphToCytoscapeJSON(graph)))
     filename <- "g.json"
     write(g.json, file=filename)
     send(obj, list(cmd="httpAddGraph", callback="handleResponse", status="request",
                    payload=filename))
     while (!browserResponseReady(obj)){
        Sys.sleep(.1)
        }
     printf("browserResponseReady")
     getBrowserResponse(obj);
     })

#----------------------------------------------------------------------------------------------------
setMethod('loadStyle', 'TReNA.Viz',

  function (obj, filename) {
     printf("TReNA.Viz::loadStyle");
     send(obj, list(cmd="httpSetStyle", callback="handleResponse", status="request", payload=filename))
     while (!browserResponseReady(obj)){
        Sys.sleep(.1)
        }
     printf("browserResponseReady")
     getBrowserResponse(obj);
     })

#----------------------------------------------------------------------------------------------------
setMethod('addBedTrackFromDataFrame', 'TReNA.Viz',

  function (obj, trackName, tbl.bed, displayMode="COLLAPSED", color) {
     printf("TReNA.Viz::addBedTrackFromDataFrame");
     temp.filename <- "tmp.bed"
     write.table(tbl.bed, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE, file=temp.filename)
     payload <- list(name=trackName, bedFileName=temp.filename, displayMode=displayMode, color=color)
     send(obj, list(cmd="addBedTrackFromDataFrame", callback="handleResponse", status="request", payload=payload))
     while (!browserResponseReady(obj)){
        Sys.sleep(.1)
        }
     printf("browserResponseReady")
     getBrowserResponse(obj);
     })

#----------------------------------------------------------------------------------------------------
setMethod('addBedTrackFromHostedFile', 'TReNA.Viz',

  function (obj, trackName, uri, index.uri, displayMode="COLLAPSED", color) {
     printf("TReNA.Viz::addBedTrackFromHostedFile");
     payload <- list(name=trackName, uri=uri, indexUri=index.uri, displayMode=displayMode, color=color)
     send(obj, list(cmd="addBedTrackFromHostedFile", callback="handleResponse",
                    status="request", payload=payload))
     while (!browserResponseReady(obj)){
        Sys.sleep(.1)
        }
     printf("browserResponseReady")
     getBrowserResponse(obj);
     })

#----------------------------------------------------------------------------------------------------
setMethod('showGenomicRegion', 'TReNA.Viz',

   function (obj, regionString) {
     payload <- list(regionString=regionString)
     send(obj, list(cmd="showGenomicRegion", callback="handleResponse", status="request", payload=payload))
     while (!browserResponseReady(obj)){
        Sys.sleep(.1)
        }
     printf("browserResponseReady")
     getBrowserResponse(obj);
     })

#----------------------------------------------------------------------------------------------------
setMethod('getGenomicRegion', 'TReNA.Viz',

   function (obj, regionString) {
     payload <- ""
     send(obj, list(cmd="getGenomicRegion", callback="handleResponse", status="request", payload=payload))
     while (!browserResponseReady(obj)){
        Sys.sleep(.1)
        }
     printf("browserResponseReady")
     getBrowserResponse(obj);
     })

#----------------------------------------------------------------------------------------------------
setMethod('fit', 'TReNA.Viz',

  function (obj, padding=30) {
     send(obj, list(cmd="fit", callback="handleResponse", status="request", payload=padding))
     while (!browserResponseReady(obj)){
        Sys.sleep(.1)
        }
     printf("browserResponseReady")
     getBrowserResponse(obj);
     })

#----------------------------------------------------------------------------------------------------
setMethod('fitSelected', 'TReNA.Viz',

  function (obj, padding=30) {
     send(obj, list(cmd="fitSelected", callback="handleResponse", status="request", payload=padding))
     while (!browserResponseReady(obj)){
        Sys.sleep(.1)
        }
     printf("browserResponseReady")
     getBrowserResponse(obj);
     })

#----------------------------------------------------------------------------------------------------
setMethod('selectNodes', 'TReNA.Viz',

  function (obj, nodeIDs) {
     payload <- list(nodeIDs=nodeIDs)
     send(obj, list(cmd="selectNodes", callback="handleResponse", status="request", payload=payload))
     while (!browserResponseReady(obj)){
        Sys.sleep(.1)
        }
     printf("browserResponseReady")
     getBrowserResponse(obj);
     })

#----------------------------------------------------------------------------------------------------
setMethod('getSelectedNodes', 'TReNA.Viz',

  function (obj) {
     payload <- ""
     send(obj, list(cmd="getSelectedNodes", callback="handleResponse", status="request", payload=payload))
     while (!browserResponseReady(obj)){
        Sys.sleep(.1)
        }
     printf("browserResponseReady")
     fromJSON(getBrowserResponse(obj))$id;
     })

#----------------------------------------------------------------------------------------------------
setMethod('sfn', 'TReNA.Viz',

  function (obj) {
     payload <- ""
     send(obj, list(cmd="sfn", callback="handleResponse", status="request", payload=payload))
     while (!browserResponseReady(obj)){
        Sys.sleep(.1)
        }
     printf("browserResponseReady")
     getBrowserResponse(obj);
     })

#----------------------------------------------------------------------------------------------------
setMethod('clearSelection', 'TReNA.Viz',

  function (obj) {
     payload <- ""
     send(obj, list(cmd="clearSelection", callback="handleResponse", status="request", payload=payload))
     while (!browserResponseReady(obj)){
        Sys.sleep(.1)
        }
     printf("browserResponseReady")
     getBrowserResponse(obj);
     })

#----------------------------------------------------------------------------------------------------
setMethod('layoutStrategies', 'TReNA.Viz',

  function (obj) {
     send(obj, list(cmd="layoutStrategies", callback="handleResponse", status="request",
                                  payload=""))
     while (!browserResponseReady(obj)){
        Sys.sleep(.1)
        }
     getBrowserResponse(obj)
     })

#----------------------------------------------------------------------------------------------------
setMethod('layout', 'TReNA.Viz',

  function (obj, strategy="random") {
     send(obj, list(cmd="doLayout", callback="handleResponse", status="request",
                                  payload=strategy))
     while (!browserResponseReady(obj)){
        Sys.sleep(.1)
        }
     getBrowserResponse(obj)
     })

#----------------------------------------------------------------------------------------------------
# {elements: [
#    {data: {id: 'a', score:5}, position: {x: 100, y: 200}},
#    {data: {id: 'b', score:100}, position: {x: 200, y: 200}},
#    {data: {id: 'e1', source: 'a', target: 'b'}}
#    ],  // elements array
# layout: { name: 'preset'},
# style: [{selector: 'node', style: {'content': 'data(id)'}}]
# }
.graphToJSON <- function(g)
{
    x <- '{"elements": [';
    nodes <- nodes(g)
    edgeNames <- edgeNames(g)
    edges <- strsplit(edgeNames, "~")  # a list of pairs
    edgeNames <- sub("~", "->", edgeNames)
    names(edges) <- edgeNames

    noa.names <- names(nodeDataDefaults(g))
    eda.names <- names(edgeDataDefaults(g))
    nodeCount <- length(nodes)
    edgeCount <- length(edgeNames)

    for(n in 1:nodeCount){
       node <- nodes[n]
       x <- sprintf('%s {"data": {"id": "%s"', x, node);
       nodeAttributeCount <- length(noa.names)
       for(i in seq_len(nodeAttributeCount)){
          noa.name <- noa.names[i];
          value <-  nodeData(g, node, noa.name)[[1]]
          if(is.numeric(value))
             x <- sprintf('%s, "%s": %s', x, noa.name, value)
          else
             x <- sprintf('%s, "%s": "%s"', x, noa.name, value)
          } # for i
       x <- sprintf('%s}', x)     # close off this node data element
       if(all(c("xPos", "yPos") %in% noa.names)){
           xPos <- as.integer(nodeData(g, node, "xPos"))
           yPos <- as.integer(nodeData(g, node, "yPos"))
           x <- sprintf('%s, "position": {"x": %d, "y": %d}', x, xPos, yPos)
           } # add position element
       x <- sprintf('%s}', x)     # close off this node data element
       if(n != nodeCount)
           x <- sprintf("%s,", x)  # another node coming, add a comma
       } # for n

    for(e in seq_len(edgeCount)) {
       edgeName <- edgeNames[e]
       edge <- edges[[e]]
       sourceNode <- edge[[1]]
       targetNode <- edge[[2]]
       x <- sprintf('%s, {"data": {"id": "%s", "source": "%s", "target": "%s"', x, edgeName, sourceNode, targetNode);
       edgeAttributeCount <- length(eda.names)
       for(i in seq_len(edgeAttributeCount)){
          eda.name <- eda.names[i];
          value <-  edgeData(g, sourceNode, targetNode, eda.name)[[1]]
          if(is.numeric(value))
             x <- sprintf('%s, "%s": %s', x, eda.name, value)
          else
             x <- sprintf('%s, "%s": "%s"', x, eda.name, value)
          } # for each edgeAttribute
       x <- sprintf('%s}}', x)     # close off this edge data element
       } # for e

    x <- sprintf("%s]}", x)

    x

} # .graphToJSON
#------------------------------------------------------------------------------------------------------------------------
myQP <- function(queryString)
{
   printf("=== TReNA-Viz::myQP");
   #print(queryString)
     # for reasons not quite clear, the query string comes in with extra characters
     # following the expected filename:
     #
     #  "?sampleStyle.js&_=1443650062946"
     #
     # check for that, cleanup the string, then see if the file can be found

   ampersand.loc <- as.integer(regexpr("&", queryString, fixed=TRUE))
   #printf("ampersand.loc: %d", ampersand.loc)

   if(ampersand.loc > 0){
      queryString <- substring(queryString, 1, ampersand.loc - 1);
      }

   questionMark.loc <- as.integer(regexpr("?", queryString, fixed=TRUE));
   #printf("questionMark.loc: %d", questionMark.loc)

   if(questionMark.loc == 1)
      queryString <- substring(queryString, 2, nchar(queryString))

   filename <- queryString;
   #printf("myQP filename: '%s'", filename)
   #printf("       exists?  %s", file.exists(filename));

   stopifnot(file.exists(filename))

   printf("--- about to scan %s", filename);
      # reconstitute linefeeds though collapsing file into one string, so json
      # structure is intact, and any "//" comment tokens only affect one line
   text <- paste(scan(filename, what=character(0), sep="\n", quiet=TRUE), collapse="\n")
   printf("%d chars read from %s", nchar(text), filename);

   return(text);

} # myQP
#----------------------------------------------------------------------------------------------------
setMethod('geneRegulatoryModelToGraph', 'TReNA.Viz',

  function (obj, target.gene, tbl.gm, tbl.reg) {

     required.geneModelColumnNames <- c("tf", "pearson", "spearman", "betaLasso", "randomForest", "pcaMax", "concordance")
     required.regulatoryRegionsColumnNames <- c("motifName", "chrom", "motifStart", "motifEnd", "strand",
                                                "motifScore", "motifRelativeScore", "match",
                                                "distance.from.tss",
                                                "chromStart", "chromEnd", "seq", "status", "tf")
     stopifnot(all(required.geneModelColumnNames %in% colnames(tbl.gm)))
     stopifnot(all(required.regulatoryRegionsColumnNames %in% colnames(tbl.reg)))

     printf("genes: %d, %d occurences of %d motifs", length(tbl.gm$tf), length(tbl.reg$motifName),
            length(unique(tbl.reg$motifName)))

     g <- graphNEL(edgemode = "directed")

     nodeDataDefaults(g, attr = "type") <- "undefined"             # targetGene, tf, footprint
     nodeDataDefaults(g, attr = "label") <- "default node label"
     nodeDataDefaults(g, attr = "distance") <- 0
     nodeDataDefaults(g, attr = "pearson") <- 0
     nodeDataDefaults(g, attr = "randomForest") <- 0
     nodeDataDefaults(g, attr = "pcaMax") <- 0
     nodeDataDefaults(g, attr = "concordance") <- 0
     nodeDataDefaults(g, attr = "betaLasso") <- 0
     nodeDataDefaults(g, attr = "motif") <- ""
     nodeDataDefaults(g, attr = "xPos") <- 0
     nodeDataDefaults(g, attr = "yPos") <- 0

     edgeDataDefaults(g, attr = "edgeType") <- "undefined"

     tfs <- tbl.gm$tf

     regRegions.names <- unlist(lapply(1:nrow(tbl.reg), function(i){
         distance.from.tss <- tbl.reg$distance.from.tss[i]
         region.size <- nchar(tbl.reg$match[i])
         motif.name <- tbl.reg$motifName[i]
         if(distance.from.tss < 0)
            sprintf("%s.fp.downstream.%05d.L%d.%s", target.gene, abs(distance.from.tss), region.size, motif.name)
          else
            sprintf("%s.fp.upstream.%05d.L%d.%s", target.gene, abs(distance.from.tss), region.size, motif.name)
          }))

   tbl.reg$regionName <- regRegions.names
   all.nodes <- unique(c(target.gene, tfs, regRegions.names))
   g <- addNode(all.nodes, g)

   nodeData(g, target.gene, "type") <- "targetGene"
   nodeData(g, tfs, "type")         <- "TF"
   nodeData(g, regRegions.names, "type")  <- "regulatoryRegion"
   nodeData(g, all.nodes, "label")  <- all.nodes
   nodeData(g, regRegions.names, "label") <- tbl.reg$motifName
   nodeData(g, regRegions.names, "distance") <- tbl.reg$distance
   nodeData(g, regRegions.names, "motif") <- tbl.reg$motifName

   nodeData(g, tfs, "pearson") <- tbl.gm$pearson
   nodeData(g, tfs, "betaLasso") <- tbl.gm$betaLasso
   nodeData(g, tfs, "randomForest") <- tbl.gm$randomForest
   nodeData(g, tfs, "pcaMax") <- tbl.gm$pcaMax
   nodeData(g, tfs, "concordance") <- tbl.gm$concordance

   #browser()
   g <- addEdge(tbl.reg$tf, tbl.reg$regionName, g)
   edgeData(g,  tbl.reg$tf, tbl.reg$regionName, "edgeType") <- "bindsTo"

   g <- graph::addEdge(tbl.reg$regionName, target.gene, g)
   edgeData(g, tbl.reg$regionName, target.gene, "edgeType") <- "regulatorySiteFor"

   g

   }) # geneRegulatoryModelToGraph

#------------------------------------------------------------------------------------------------------------------------
# {elements: [
#    {data: {id: 'a', score:5}, position: {x: 100, y: 200}},
#    {data: {id: 'b', score:100}, position: {x: 200, y: 200}},
#    {data: {id: 'e1', source: 'a', target: 'b'}}
#    ],  // elements array
# layout: { name: 'preset'},
# style: [{selector: 'node', style: {'content': 'data(id)'}}]
# }
graphToJSON <- function(g)
{
    x <- '{"elements": [';
    nodes <- nodes(g)
    edgeNames <- edgeNames(g)
    edges <- strsplit(edgeNames, "~")  # a list of pairs
    edgeNames <- sub("~", "->", edgeNames)
    names(edges) <- edgeNames

    noa.names <- names(nodeDataDefaults(g))
    eda.names <- names(edgeDataDefaults(g))
    nodeCount <- length(nodes)
    edgeCount <- length(edgeNames)

    for(n in 1:nodeCount){
       #printf("---- node %d", n)
       node <- nodes[n]
       x <- sprintf('%s {"data": {"id": "%s"', x, node);
       nodeAttributeCount <- length(noa.names)
       for(i in seq_len(nodeAttributeCount)){
          noa.name <- noa.names[i];
          #printf("node %s, noa.name: %s", node, noa.name)
          value <-  nodeData(g, node, noa.name)[[1]]
          if(is.numeric(value))
             x <- sprintf('%s, "%s": %s', x, noa.name, value)
          else
             x <- sprintf('%s, "%s": "%s"', x, noa.name, value)
          #browser();
          #xyz <- 99
          } # for i
       x <- sprintf('%s}', x)     # close off this node data element
       #printf("-- x partway: %s", x)
       if(all(c("xPos", "yPos") %in% noa.names)){
           xPos <- as.integer(nodeData(g, node, "xPos"))
           yPos <- as.integer(nodeData(g, node, "yPos"))
           x <- sprintf('%s, "position": {"x": %d, "y": %d}', x, xPos, yPos)
           #xyz <- 99
           } # add position element
       x <- sprintf('%s}', x)     # close off this node data element
       if(n != nodeCount)
           x <- sprintf("%s,", x)  # another node coming, add a comma
       } # for n

    #browser()
    #xyz <- 99

    for(e in seq_len(edgeCount)) {
       edgeName <- edgeNames[e]
       edge <- edges[[e]]
       sourceNode <- edge[[1]]
       targetNode <- edge[[2]]
       x <- sprintf('%s, {"data": {"id": "%s", "source": "%s", "target": "%s"', x, edgeName, sourceNode, targetNode);
       edgeAttributeCount <- length(eda.names)
       for(i in seq_len(edgeAttributeCount)){
          eda.name <- eda.names[i];
          value <-  edgeData(g, sourceNode, targetNode, eda.name)[[1]]
          if(is.numeric(value))
             x <- sprintf('%s, "%s": %s', x, eda.name, value)
          else
             x <- sprintf('%s, "%s": "%s"', x, eda.name, value)
          } # for each edgeAttribute
       x <- sprintf('%s}}', x)     # close off this edge data element
       } # for e

    x <- sprintf("%s]}", x)

    x

} # graphToJSON
#------------------------------------------------------------------------------------------------------------------------
addGeneModelLayout <- function(g)
{
   printf("--- addGeneModelLayout")
   all.distances <- sort(unique(unlist(nodeData(g, attr='distance'), use.names=FALSE)))
   print(all.distances)

   fp.nodes <- nodes(g)[which(unlist(nodeData(g, attr="type"), use.names=FALSE) == "regulatoryRegion")]
   tf.nodes <- nodes(g)[which(unlist(nodeData(g, attr="type"), use.names=FALSE) == "TF")]
   targetGene.nodes <- nodes(g)[which(unlist(nodeData(g, attr="type"), use.names=FALSE) == "targetGene")]

     # add in a zero in case all of the footprints are up or downstream of the 0 coordinate, the TSS
   span.endpoints <- range(c(0, as.numeric(nodeData(g, fp.nodes, attr="distance"))))
   span <- max(span.endpoints) - min(span.endpoints)
   footprintLayoutFactor <- 1
   if(span < 600)  #
       footprintLayoutFactor <- 600/span

   xPos <- as.numeric(nodeData(g, fp.nodes, attr="distance")) * footprintLayoutFactor
   yPos <- 0
   nodeData(g, fp.nodes, "xPos") <- xPos
   nodeData(g, fp.nodes, "yPos") <- yPos

   adjusted.span.endpoints <- range(c(0, as.numeric(nodeData(g, fp.nodes, attr="xPos"))))
   printf("raw span of footprints: %d   footprintLayoutFactor: %f  new span: %8.0f",
          span, footprintLayoutFactor, abs(max(adjusted.span.endpoints) - min(adjusted.span.endpoints)))

   tfs <- names(which(nodeData(g, attr="type") == "TF"))

   for(tf in tfs){
      footprint.neighbors <- edges(g)[[tf]]
      footprint.positions <- as.integer(nodeData(g, footprint.neighbors, attr="xPos"))
      new.xPos <- mean(footprint.positions)
      #printf("%8s: %5d", tf, new.xPos)
      nodeData(g, tf, "xPos") <- new.xPos
      nodeData(g, tf, "yPos") <- sample(300:1200, 1)
      } # for tf

   nodeData(g, targetGene.nodes, "xPos") <- 0
   nodeData(g, targetGene.nodes, "yPos") <- -200

   g

} # addGeneModelLayout
#------------------------------------------------------------------------------------------------------------------------

