#------------------------------------------------------------------------------------------------------------------------
.Solver <- setClass ("Solver",
                     representation = representation(mtx.assay="matrix",
                                                     quiet="logical",
                                                     state="environment")
                     )

#------------------------------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#------------------------------------------------------------------------------------------------------------------------
setGeneric("getSolverName",   signature="obj", function(obj, target.gene, tfs) standardGeneric ("getSolverName"))
setGeneric("run",             signature="obj", function(obj, target.gene, tfs,
                                                        splitDataPercentage=100.0) standardGeneric ("run"))
setGeneric("getAssayData",    signature="obj", function(obj) standardGeneric ("getAssayData"))
#------------------------------------------------------------------------------------------------------------------------
Solver <- function(mtx.assay=matrix(), quiet=TRUE)
{
   env <- new.env(parent=emptyenv())
  .Solver(mtx.assay=mtx.assay, quiet=quiet, state=env)

} # TReNA, the constructor
#------------------------------------------------------------------------------------------------------------------------
setMethod("getAssayData", "Solver",

   function (obj){
      obj@mtx.assay
      })

#------------------------------------------------------------------------------------------------------------------------
