#' An S4 class to represent a TReNA object
#'
#' @name TReNA-class
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
#------------------------------------------------------------------------------------------------------------------------
TReNA <- function(mtx.assay=matrix(), solver="lasso", quiet=TRUE)
{
  stopifnot(solver %in% c("lasso", "randomForest", "bayesSpike", "pearson", "spearman","sqrtlasso","lassopv","ridge"))

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
  else if(solver == "sqrtlasso")
      solver <- SqrtLassoSolver(mtx.assay)
  else if(solver == "lassopv")
      solver <- LassoPVSolver(mtx.assay)
  else if(solver == "ridge")
      solver <- RidgeSolver(mtx.assay)

  .TReNA(solver=solver, quiet=quiet)

} # TReNA, the constructor
#------------------------------------------------------------------------------------------------------------------------
#' Solve the TReNA object
#' @name solve
#'
#' @param obj An object of class TReNA
#' @param target.gene A designated target gene that should be part of the mtx.assay data
#' @param tfs The designated set of transcription factors that could be associated with the target gene.
#' @param tf.weights A set of weights on the transcription factors (default = rep(1, length(tfs)))
#' @param extraArgs Modifiers to the solver
#'
#' @return A data frame containing coefficients relating the target gene to each transcription factor

setMethod("solve", "TReNA",

   function (obj, target.gene, tfs, tf.weights=rep(1, length(tfs)), extraArgs=list()){
      # printf("entering TReNA::solve")
      run(obj@solver, target.gene, tfs, tf.weights, extraArgs)
      })
#------------------------------------------------------------------------------------------------------------------------
