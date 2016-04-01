#------------------------------------------------------------------------------------------------------------------------
.TReNA <- setClass ("TReNA",
                    representation = representation(solver="Solver",
                                                    quiet="logical")
                    )

#------------------------------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#------------------------------------------------------------------------------------------------------------------------
setGeneric("getModel",                 signature="obj", function(obj) standardGeneric ("getModel"))
setGeneric("solve",                    signature="obj", function(obj, target.gene, tfs) standardGeneric ("solve"))
#------------------------------------------------------------------------------------------------------------------------
TReNA <- function(mtx.assay=matrix(), solver="lasso", quiet=TRUE)
{
  stopifnot(solver %in% c("lasso", "randomForest"))

  if(solver == "lasso")
      solver <- LassoSolver(mtx.assay)
  else if(solver == "randomForest")
      solver <- RandomForestSolver(mtx.assay)

  .TReNA(solver=solver, quiet=quiet)

} # TReNA, the constructor
#------------------------------------------------------------------------------------------------------------------------
setMethod("solve", "TReNA",

   function (obj, target.gene, tfs){
       run(obj@solver, target.gene, tfs)
   })

#------------------------------------------------------------------------------------------------------------------------
