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

   tv <- TReNA.Viz(portRange=11011)
   TReNA:::addGraph(tv, g)
   layout(tv, "grid")

   tv

} # test_addDemoGraph
#----------------------------------------------------------------------------------------------------
test_.graphToJSON <- function()
{
   printf("--- test.graphToJSON")

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
