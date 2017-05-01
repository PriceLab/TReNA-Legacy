.Recipe <- setClass('Recipe',
                       representation(name="character",
                                      targetGene="character",
                                      genome="character",
                                      candidateFilterSpec="list",
                                      solverSpec="list",
                                      variants="character",
                                      minid="character"
                                      )
                        )
#------------------------------------------------------------------------------------------------------------------------
setGeneric("getRecipeName",              signature="obj", function(obj) standardGeneric ("getRecipeName"))
setGeneric("getTargetGene",              signature="obj", function(obj) standardGeneric ("getTargetGene"))
setGeneric("getGenome",                  signature="obj", function(obj) standardGeneric ("getGenome"))
setGeneric("getCandidateFilterSpec",     signature="obj", function(obj) standardGeneric ("getCandidateFilterSpec"))
setGeneric("getSolverSpec",              signature="obj", function(obj) standardGeneric ("getSolverSpec"))
setGeneric("getVariants",                signature="obj", function(obj) standardGeneric ("getVariants"))
setGeneric("getMinid",                   signature="obj", function(obj) standardGeneric ("getMinid"))
#------------------------------------------------------------------------------------------------------------------------
Recipe <- function(name=NA_character_,
                   targetGene=NA_character_,
                   genome=NA_character_,
                   candidateFilterSpec=list(),
                   solverSpec=list(),
                   variants=NA_character_,
                   minid=NA_character_)
{
   .Recipe(name=name, targetGene=targetGene, genome=genome,
           candidateFilterSpec=candidateFilterSpec,
           solverSpec=solverSpec, variants=variants, minid=minid)

} # Recipe constructor
#------------------------------------------------------------------------------------------------------------------------
setMethod("getRecipeName", "Recipe",
          function(obj){
             return(obj@name)
          })
#------------------------------------------------------------------------------------------------------------------------
setMethod("getTargetGene", "Recipe",
          function(obj){
             return(obj@targetGene)
          })
#------------------------------------------------------------------------------------------------------------------------
setMethod("getGenome", "Recipe",
          function(obj){
             return(obj@genome)
          })
#------------------------------------------------------------------------------------------------------------------------
setMethod("getCandidateFilterSpec", "Recipe",
          function(obj){
             return(obj@candidateFilterSpec)
          })
#------------------------------------------------------------------------------------------------------------------------
setMethod("getSolverSpec", "Recipe",
          function(obj){
             return(obj@solverSpec)
          })
#------------------------------------------------------------------------------------------------------------------------
setMethod("getVariants", "Recipe",
          function(obj){
             return(obj@variants)
          })
#------------------------------------------------------------------------------------------------------------------------
setMethod("getMinid", "Recipe",
          function(obj){
             return(obj@minid)
          })
#------------------------------------------------------------------------------------------------------------------------
