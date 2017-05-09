library(TReNA)
library(RUnit)
#----------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
   printf("--- test_Recipe.R")
   test_emptyConstructor()
   test_constructorWithParametersSpecified()
   test_cookRandomForestRecipe()

} # runTests
#----------------------------------------------------------------------------------------------------
test_emptyConstructor <- function()
{
   printf("--- test_emptyConstructor")
   recipe <- Recipe();
   checkEquals(is(recipe), "Recipe")

   checkTrue(is.na(getRecipeName(recipe)))
   checkTrue(is.na(getGenome(recipe)))
   checkEquals(getCandidateFilterSpec(recipe), list())
   checkEquals(getSolverSpec(recipe), list())
   checkTrue(is.na(getVariants(recipe)))
   checkTrue(is.na(getMinid(recipe)))

} # test_emptyConstructor
#----------------------------------------------------------------------------------------------------
test_constructorWithParametersSpecified <- function()
{
   printf("--- test_constructorWithParametersSpecified")

   name <- "mef2c.brain.hint"
   targetGene  <-  "MEF2C"
   genome <- "hg38"
      # build up the information needed for a footprint filter using TReNA package's built in chr5
      # brain hint footprint database
   db.address <- system.file(package="TReNA", "extdata")
   genome.db.uri    <- paste("sqlite:/", db.address, "genome.sub.db",  sep = "/")
   footprint.db.uri <- paste("sqlite:/", db.address, "project.sub.db", sep = "/")
   target.gene <- "MEF2C"
   promoter.length <- 1000

   candidateFilterSpec <- list(genomeDB=genome.db.uri,
                                footprintDB=footprint.db.uri,
                                geneCenteredSpec=list(
                                   targetGene=target.gene,
                                   tssUpstream=promoter.length,
                                   tssDownstream=promoter.length),
                                regionsSpec=c())
   candidateFilter.json <- toJSON(candidateFilterSpec, auto_unbox=TRUE)

   solverSpec <- list(solver="randomForest",
                       matrix="rosmap",
                       targetGene=target.gene,
                       candidateRegulators=NA_character_ # supplied by candidateFilter
                       )

   variants <- c("rs2710873", "rs7526076", "rs11584349", "rs4970401", "rs74048003")
   minid <- "minid.012345"

   recipe <- Recipe(name=name,
                    targetGene=targetGene,
                    genome=genome,
                    candidateFilter=candidateFilterSpec,
                    solver=solverSpec,
                    variants=variants,
                    minid=minid)

   checkEquals(is(recipe), "Recipe")

   checkEquals(getRecipeName(recipe), name)
   checkEquals(getTargetGene(recipe), targetGene)
   checkEquals(getGenome(recipe), genome)
   checkEquals(getCandidateFilterSpec(recipe), candidateFilterSpec)
   checkEquals(getSolverSpec(recipe), solverSpec)
   # checkEquals(getAssayMatrix(recipe), assayMatrix)
   checkEquals(getVariants(recipe), variants)
   checkEquals(getMinid(recipe), minid)
   # checkEquals(getRegions(recipe), regions)

} # test_constructorWithParametersSpecified
#----------------------------------------------------------------------------------------------------
test_cookRandomForestRecipe <- function()
{
   printf("--- test_cookRandomForestRecipe")

} # test_cookRandomForestRecipe
#----------------------------------------------------------------------------------------------------
