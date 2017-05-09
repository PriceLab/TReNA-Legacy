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
   test_cookTwoSimilarRandomForestRecipes()
   test_cookTwoDifferentRandomForestRecipes()

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
# reuse this in several tests
create.mef2c.sqlite.recipe <- function(geneCentered=TRUE, promoter.length=1000)
{
   recipe.name <- "mef2c.brain.hint.sqlite-geneCentered"
   targetGene  <-  "MEF2C"
   genome <- "hg38"
      # build up the information needed for a footprint filter using TReNA package's built in chr5
      # brain hint footprint database
   db.address <- system.file(package="TReNA", "extdata")
   genome.db.uri    <- paste("sqlite:/", db.address, "genome.sub.db",  sep = "/")
   linking.data.uri <- paste("sqlite:/", db.address, "project.sub.db", sep = "/")

   target.gene <- "MEF2C"

   chromosome <- "chr5"
   tss <- 88904257

   candidateFilterSpec <- list(filterType="DNaseFootprints",
                               genomeDB=genome.db.uri,
                               linkingDataURI=linking.data.uri,  # e.g., linking TFs via DNase footprints to genes
                               geneCenteredSpec=list(
                                  targetGene=target.gene,
                                  tssUpstream=promoter.length,
                                  tssDownstream=promoter.length),
                                regionsSpec=c())
   if(!geneCentered){
      recipe.name <- "mef2c.brain.hint.sqlite-regionSpecified"
      candidateFilterSpec$geneCenteredSpec=list()
      candidateFilterSpec$regionsSpec=c(sprintf("%s:%d-%d", chromosome, tss-1000, tss+1000))
      } # !geneCentered: specify region explicitly


   solverSpec <- list(solver="randomForest",
                       matrixName="rosmap",
                       targetGene=target.gene,
                       candidateRegulators=NA_character_ # supplied by candidateFilter
                       )

   variants <- c("rs2710873", "rs7526076", "rs11584349", "rs4970401", "rs74048003")
   minid <- "minid.012345"

   recipe <- Recipe(name=recipe.name,
                    targetGene=targetGene,
                    genome=genome,
                    candidateFilter=candidateFilterSpec,
                    solver=solverSpec,
                    variants=variants,
                    minid=minid)
   recipe

} # create.mef2c.sqlite.recipe
#----------------------------------------------------------------------------------------------------
test_constructorWithParametersSpecified <- function()
{
   printf("--- test_constructorWithParametersSpecified")

   recipe <- create.mef2c.sqlite.recipe()
   checkEquals(is(recipe), "Recipe")

   checkEquals(getRecipeName(recipe), "mef2c.brain.hint.sqlite")
   checkEquals(getTargetGene(recipe), "MEF2C")
   checkEquals(getGenome(recipe), "hg38")

   cfSpec <- getCandidateFilterSpec(recipe)
   checkEquals(cfSpec$filterType, "DNaseFootprints")
   checkEquals(sort(names(cfSpec)),
               c("filterType", "geneCenteredSpec", "genomeDB", "linkingDataURI", "regionsSpec"))

   expected.uri <- sprintf("sqlite://%s/%s/%s", system.file(package="TReNA"), "extdata", "genome.sub.db")
   checkEquals(cfSpec$genomeDB, expected.uri)

   expected.uri <- sprintf("sqlite://%s/%s/%s", system.file(package="TReNA"), "extdata", "project.sub.db")
   checkEquals(cfSpec$linkingDataURI, expected.uri)


   checkTrue(is.null(cfSpec$regionsSpec))  # we used a gene-centered

   solverSpec <- getSolverSpec(recipe)
   checkEquals(solverSpec$solver, "randomForest")
   checkEquals(solverSpec$matrixName, "rosmap")
   checkEquals(solverSpec$targetGene, "MEF2C")
   checkTrue(is.na(solverSpec$candidateRegulators))

   checkEquals(getSolverSpec(recipe), solverSpec)
   checkEquals(getVariants(recipe), c("rs2710873", "rs7526076", "rs11584349", "rs4970401", "rs74048003"))
   checkEquals(getMinid(recipe), "minid.012345")

} # test_constructorWithParametersSpecified
#----------------------------------------------------------------------------------------------------
test_cookRandomForestRecipe <- function()
{
   printf("--- test_cookRandomForestRecipe")

   recipe <- create.mef2c.sqlite.recipe()
   recipe.name <- getRecipeName(recipe)

   cook <- Cook(recipes=list(recipe))
   model <- createModelsFromRecipes(cook)

   checkEquals(names(model), recipe.name)
   model.1 <- model[[recipe.name]]

   checkEquals(sort(names(model.1)), c("edges", "r2"))
   tbl.tfs <- subset(model.1$edges, IncNodePurity > 1)

       # a quick sanity check
       # FOXO4 comes out on top in this recipe.  GeneCards
       # reports it as an inhibitor of MEF2C
       # http://www.genecards.org/cgi-bin/carddisp.pl?gene=FOXO4

   checkEquals(rownames(tbl.tfs)[1], "FOXO4")

   checkEqualsNumeric(tbl.tfs$IncNodePurity[1], 10.52238, tol=1)
   checkEqualsNumeric(tbl.tfs$gene.cor[1], -0.7013734, tol=1e-1)

} # test_cookRandomForestRecipe
#----------------------------------------------------------------------------------------------------
# two very similar recipes, differing only in the argument to the included CandidateFilterSpec
# field: the region of interest for the first - the area in which to look for footprints -
# is specified implicitly as a gene and its up- and downstream promoter region from the tss,
# which is looked up.  the second recipe provides explicit chrN:start-end region (as a string)
# we expect both calculated models to be identical modulo small random numerical differences
test_cookTwoSimilarRandomForestRecipes <- function()
{
   printf("--- test_cookTwoSimilarRandomForestRecipes")

   recipe.1 <- create.mef2c.sqlite.recipe(geneCentered=TRUE)
   recipe.1.name <- getRecipeName(recipe.1)

   recipe.2 <- create.mef2c.sqlite.recipe(geneCentered=FALSE)
   recipe.2.name <- getRecipeName(recipe.2)

   cook <- Cook(recipes=list(recipe.1, recipe.2))
   model <- createModelsFromRecipes(cook)

   checkEquals(names(model),
               c("mef2c.brain.hint.sqlite-geneCentered", "mef2c.brain.hint.sqlite-regionSpecified"))

   model.1 <- model[[recipe.1.name]]

   checkEquals(sort(names(model.1)), c("edges", "r2"))
   tbl.tfs <- subset(model.1$edges, IncNodePurity > 1)

       # a quick sanity check
       # FOXO4 comes out on top in this recipe.  GeneCards
       # reports it as an inhibitor of MEF2C
       # http://www.genecards.org/cgi-bin/carddisp.pl?gene=FOXO4

   checkEquals(rownames(tbl.tfs)[1], "FOXO4")
   checkEqualsNumeric(tbl.tfs$IncNodePurity[1], 10.52238, tol=1)
   checkEqualsNumeric(tbl.tfs$gene.cor[1], -0.7013734, tol=1e-1)

   model.2 <- model[[recipe.2.name]]

   checkEquals(sort(names(model.2)), c("edges", "r2"))
   tbl.tfs <- subset(model.2$edges, IncNodePurity > 1)

       # a quick sanity check
       # FOXO4 comes out on top in this recipe.  GeneCards
       # reports it as an inhibitor of MEF2C
       # http://www.genecards.org/cgi-bin/carddisp.pl?gene=FOXO4

   checkEquals(rownames(tbl.tfs)[1], "FOXO4")
   checkEqualsNumeric(tbl.tfs$IncNodePurity[1], 10.52238, tol=1)
   checkEqualsNumeric(tbl.tfs$gene.cor[1], -0.7013734, tol=1e-1)

} # test_cookTwoSimilarRandomForestRecipes
#----------------------------------------------------------------------------------------------------
# two very similar recipes, differing only in the argument to the included CandidateFilterSpec
# field: the region of interest for the first - the area in which to look for footprints -
# is specified implicitly as a gene and its up- and downstream promoter region from the tss,
# which is looked up.  the second recipe provides explicit chrN:start-end region (as a string)
# we expect both calculated models to be identical modulo small random numerical differences
test_cookTwoDifferentRandomForestRecipes <- function()
{
   printf("--- test_cookTwoDifferentRandomForestRecipes")

   recipe.1 <- create.mef2c.sqlite.recipe(geneCentered=TRUE, promoter.length=2000)
   recipe.1.name <- getRecipeName(recipe.1)

   recipe.2 <- create.mef2c.sqlite.recipe(geneCentered=FALSE)
   recipe.2.name <- getRecipeName(recipe.2)

   cook <- Cook(recipes=list(recipe.1, recipe.2))
   models <- createModelsFromRecipes(cook)

   checkEquals(names(models),
               c("mef2c.brain.hint.sqlite-geneCentered", "mef2c.brain.hint.sqlite-regionSpecified"))

   model.1 <- models[[recipe.1.name]]

   checkEquals(sort(names(model.1)), c("edges", "r2"))
   tbl.tfs <- subset(model.1$edges, IncNodePurity > 1)

       # a quick sanity check
       # FOXO4 comes out on top in this recipe.  GeneCards
       # reports it as an inhibitor of MEF2C
       # http://www.genecards.org/cgi-bin/carddisp.pl?gene=FOXO4

   checkEquals(rownames(tbl.tfs)[1], "HLF")
   checkEqualsNumeric(tbl.tfs$IncNodePurity[1], 11.9, tol=1)
   checkEqualsNumeric(tbl.tfs$gene.cor[1], 0.7459, tol=1e-1)

   model.2 <- models[[recipe.2.name]]

   checkEquals(sort(names(model.2)), c("edges", "r2"))
   tbl.tfs <- subset(model.2$edges, IncNodePurity > 1)

       # a quick sanity check
       # FOXO4 comes out on top in this recipe.  GeneCards
       # reports it as an inhibitor of MEF2C
       # http://www.genecards.org/cgi-bin/carddisp.pl?gene=FOXO4

   checkEquals(rownames(tbl.tfs)[1], "FOXO4")
   checkEqualsNumeric(tbl.tfs$IncNodePurity[1], 10.52238, tol=1)
   checkEqualsNumeric(tbl.tfs$gene.cor[1], -0.7013734, tol=1e-1)

} # test_cookTwoDifferentRandomForestRecipes
#----------------------------------------------------------------------------------------------------
