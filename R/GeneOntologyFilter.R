#----------------------------------------------------------------------------------------------------
.GeneOntologyFilter <- setClass("GeneOntologyFilter",
                            contains="CandidateFilter",
                            representation (organismDatabase="OrgDb",
                                            regulatory.genes="character",
                                            quiet="logical"))

#----------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#----------------------------------------------------------------------------------------------------
GeneOntologyFilter <- function(organismDatabase,  quiet=TRUE)
{
    data.file <- system.file(package="TReNA", "extdata", "human.regulatory.genes.RData")
    load(data.file)
    regulatory.genes <- sort(unique(unlist(all.regulatory.genes, use.names=FALSE)))
    .GeneOntologyFilter(organismDatabase=organismDatabase,
                        regulatory.genes=regulatory.genes, quiet=quiet)

} # GeneOntologyFilter, the constructor
#----------------------------------------------------------------------------------------------------
setMethod("getCandidates", "GeneOntologyFilter",

    function(obj, argsList){

          # Collect arguments from argsList
        expected.args <- c("mode")
        missing.args <- setdiff(names(argsList), expected.args)
        if(length(missing.args) > 0){
            printf("GeneOntologyFilter::getCandidates, missing fields in argsList: %s",
                   paste(missing.args, collapse=","))
            stop()
            }

       mode <- argsList[["mode"]]
       stopifnot(mode %in% c("direct")) # , "calculate"))
       list(tbl=data.frame(), tfs=obj@regulatory.genes)
    }) # getCandidates

#----------------------------------------------------------------------------------------------------
