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
setGeneric('addGraph',         signature='obj', function(obj, graph) standardGeneric ('addGraph'))
setGeneric('fit',              signature='obj', function(obj, padding=30) standardGeneric('fit'))
setGeneric('fitSelected',      signature='obj', function(obj, padding=30) standardGeneric('fitSelectedContent'))
setGeneric('selectNodes',      signature='obj', function(obj, nodeIDs) standardGeneric('selectNodes'))
setGeneric('sfn',              signature='obj', function(obj) standardGeneric('sfn'))
setGeneric('addBedTrack',      signature='obj', function(obj, trackName, tbl.bed) standardGeneric('addBedTrack'))


setGeneric('layout',              signature='obj', function(obj, strategy) standardGeneric('layout'))
setGeneric('layoutStrategies',    signature='obj', function(obj) standardGeneric('layoutStrategies'))
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

  function (obj, graph) {
     printf("TReNA.Viz::addGraph");
     print(graph)
     g.json <- .graphToJSON(graph)
     printf("about to send g.json: %d chars", nchar(g.json));
     send(obj, list(cmd="addGraph", callback="handleResponse", status="request",
                    payload=g.json))
     while (!browserResponseReady(obj)){
        Sys.sleep(.1)
        }
     printf("browserResponseReady")
     getBrowserResponse(obj);
     })

#----------------------------------------------------------------------------------------------------
setMethod('addBedTrack', 'TReNA.Viz',

  function (obj, trackName, tbl.bed) {
     printf("TReNA.Viz::addBedTrack");
     temp.filename <- "tmp.bed"
     write.table(tbl.bed, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE, file=temp.filename)
     payload <- list(name=trackName, bedFileName=temp.filename)
     send(obj, list(cmd="addBedTrack", callback="handleResponse", status="request", payload=payload))
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

