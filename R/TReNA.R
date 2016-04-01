#------------------------------------------------------------------------------------------------------------------------
.TReNA <- setClass ("TReNA",
                    representation = representation(solver="Solver",
                                                    quiet="logical")
                    )

#------------------------------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#------------------------------------------------------------------------------------------------------------------------
# setGeneric("getModel",                 signature="obj", function(obj) standardGeneric ("getModel"))
setGeneric("solve",                    signature="obj", function(obj, target.gene, tfs) standardGeneric ("solve"))
setGeneric("trainModel",               signature="obj", function(obj, target.gene, tfs, training.samples) standardGeneric ("trainModel"))
setGeneric("predictFromModel",         signature="obj", function(obj, model, tfs, test.samples) standardGeneric ("predictFromModel"))
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
setMethod("trainModel", "TReNA",

   function (obj, target.gene, tfs, training.samples){
       fit <- trainModel(obj@solver, target.gene, tfs, training.samples)
       return(fit)
   })

#------------------------------------------------------------------------------------------------------------------------
setMethod("predictFromModel", "TReNA",

   function (obj, model, tfs, test.samples){
       prediction <- predictFromModel(obj@solver, model, tfs, test.samples)
       return(prediction)
       })

#------------------------------------------------------------------------------------------------------------------------
