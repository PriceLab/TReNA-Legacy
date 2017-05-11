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
setGeneric('addGraph',     signature='obj', function(obj, graph) standardGeneric ('addGraph'))
setGeneric('fit',          signature='obj', function(obj, graph) standardGeneric ('fit'))
setGeneric('fitSelected',  signature='obj', function(obj, graph) standardGeneric ('fitSelected'))
setGeneric('layout',              signature='obj', function(obj, strategy) standardGeneric('layout'))
setGeneric('layoutStrategies',    signature='obj', function(obj) standardGeneric('layoutStrategies'))
#----------------------------------------------------------------------------------------------------
# constructor
TReNA.Viz = function(portRange=11000:11025, host="localhost", title="TReNA-Viz", quiet=TRUE)
{

   model <- data.frame()
   obj <- .TReNA.Viz(BrowserViz(portRange=portRange, host=host, title=title,
                                quiet=quiet, browserFile=cyjsBrowserFile,
                                httpQueryProcessingFunction=NULL),
                    model=model)

  while (!browserResponseReady(obj)){
      Sys.sleep(.1)
      }
   if(!quiet) {
      message(sprintf("BrowserViz ctor called from RCyjs ctor got browser response"))
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
setMethod('fit', 'TReNA.Viz',

  function (obj) {
     send(obj, list(cmd="fit", callback="handleResponse", status="request", payload=""))
     while (!browserResponseReady(obj)){
        Sys.sleep(.1)
        }
     printf("browserResponseReady")
     getBrowserResponse(obj);
     })

#----------------------------------------------------------------------------------------------------
setMethod('fitSelected', 'TReNA.Viz',

  function (obj) {
     send(obj, list(cmd="fitSelected", callback="handleResponse", status="request", payload=""))
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


