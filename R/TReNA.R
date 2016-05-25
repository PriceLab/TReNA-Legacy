#------------------------------------------------------------------------------------------------------------------------
.TReNA <- setClass ("TReNA",
                    representation = representation(solver="Solver",
                                                    quiet="logical")
                    )

#------------------------------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#------------------------------------------------------------------------------------------------------------------------
setGeneric("solve",                    signature="obj", function(obj, target.gene, tfs,
                                                                 tf.weights=rep(1, length(tfs)) , ... )
                                                           standardGeneric ("solve"))
setGeneric("trainModel",               signature="obj", function(obj, target.gene, tfs, training.samples,
                                                                 tf.weights=rep(1, length(tfs)))
                                                           standardGeneric ("trainModel"))
setGeneric("predictFromModel",         signature="obj", function(obj, model, tfs, test.samples)
                                                           standardGeneric ("predictFromModel"))
#------------------------------------------------------------------------------------------------------------------------
TReNA <- function(mtx.assay=matrix(), solver="lasso", quiet=TRUE)
{
  stopifnot(solver %in% c("lasso", "randomForest", "bayesSpike"))

  if(solver == "lasso")
      solver <- LassoSolver(mtx.assay)
  else if(solver == "randomForest")
      solver <- RandomForestSolver(mtx.assay)
  else if(solver == "bayesSpike")
      solver <- BayesSpikeSolver(mtx.assay)

  .TReNA(solver=solver, quiet=quiet)

} # TReNA, the constructor
#------------------------------------------------------------------------------------------------------------------------
setMethod("solve", "TReNA",

   function (obj, target.gene, tfs, tf.weights=rep(1, length(tfs))){
      run(obj@solver, target.gene, tfs, tf.weights , ... )
      })

#------------------------------------------------------------------------------------------------------------------------
setMethod("trainModel", "TReNA",

   function (obj, target.gene, tfs, training.samples, tf.weights=rep(1, length(tfs))){
      fit <- trainModel(obj@solver, target.gene, tfs, training.samples, tf.weights)
      return(fit)
      })

#------------------------------------------------------------------------------------------------------------------------
setMethod("predictFromModel", "TReNA",

   function (obj, model, tfs, test.samples){
      prediction <- predictFromModel(obj@solver, model, tfs, test.samples)
      return(prediction)
      })

#------------------------------------------------------------------------------------------------------------------------
