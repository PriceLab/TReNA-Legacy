#' The central class of the TReNA package
#'
#' @param mtx.assay An assay matrix of gene expression data
#' @param solver A string matching the designated solver for relating a target gene to transcription factors. (default = "lasso")
#'
#' @return An object of the TReNA class

#------------------------------------------------------------------------------------------------------------------------
.TReNA <- setClass ("TReNA",
                    representation = representation(solver="Solver",
                                                    quiet="logical")
                    )

#------------------------------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#------------------------------------------------------------------------------------------------------------------------
setGeneric("solve",                    signature="obj", function(obj, target.gene, tfs,
                                                                 tf.weights=rep(1, length(tfs)), extraArgs=list())
                                                           standardGeneric ("solve"))
setGeneric("trainModel",               signature="obj", function(obj, target.gene, tfs, training.samples,
                                                                 tf.weights=rep(1, length(tfs)))
                                                           standardGeneric ("trainModel"))
setGeneric("predictFromModel",         signature="obj", function(obj, model, tfs, test.samples)
                                                           standardGeneric ("predictFromModel"))
#------------------------------------------------------------------------------------------------------------------------
TReNA <- function(mtx.assay=matrix(), solver="lasso", quiet=TRUE)
{
  stopifnot(solver %in% c("lasso", "randomForest", "bayesSpike", "pearson", "spearman"))

  if(solver == "lasso")
      solver <- LassoSolver(mtx.assay)
  else if(solver == "randomForest")
      solver <- RandomForestSolver(mtx.assay)
  else if(solver == "bayesSpike")
      solver <- BayesSpikeSolver(mtx.assay)
  else if(solver == "pearson")
      solver <- PearsonSolver(mtx.assay)
  else if(solver == "spearman")
      solver <- SpearmanSolver(mtx.assay)

  .TReNA(solver=solver, quiet=quiet)

} # TReNA, the constructor
#------------------------------------------------------------------------------------------------------------------------
setMethod("solve", "TReNA",

   function (obj, target.gene, tfs, tf.weights=rep(1, length(tfs)), extraArgs=list()){
      # printf("entering TReNA::solve")
      run(obj@solver, target.gene, tfs, tf.weights, extraArgs)
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
