#' @title Class TReNA
#'
#' @description
#' Class \code{TReNA} defines a TReNA object and contains an assay matrix, which contains expression data over a set of
#' samples for a group of genes, and a string representing the name of a chosen solver. 
#'
#' @name TReNA-class
#' @rdname TReNA-class
#' @aliases TReNA
#' 
#' @slot mtx.assay An assay matrix of gene expression data
#' @slot solver A string matching the designated solver for relating a target gene to transcription factors. 
#' TReNA currently supports 9 solver strings (default = "lasso"):
#' \itemize{
#' \item{\link[=solve.Lasso]{"lasso"}}
#' \item{\link[=solve.Ridge]{"ridge"}}
#' \item{\link[=solve.RandomForest]{"randomForest"}}
#' \item{\link[=solve.BayesSpike]{"bayesSpike"}}
#' \item{\link[=solve.SqrtLasso]{"sqrtlasso"}}
#' \item{\link[=solve.LassoPV]{"lassopv"}}
#' \item{\link[=solve.Pearson]{"pearson"}}
#' \item{\link[=solve.Spearman]{"spearman"}}
#' \item{\link[=solve.Ensemble]{"ensemble"}}
#' }
#' 
#' @slot quiet A logical denoting whether or not the TReNA object should print output
#' 
#' @return An object of the TReNA class
#'
#' @seealso \code{\link{solve}}

#----------------------------------------------------------------------------------------------------
.TReNA <- setClass ("TReNA",
                    representation = representation(solver="Solver",
                                                    quiet="logical")
                    )

#----------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#----------------------------------------------------------------------------------------------------

setGeneric("solve",                    signature="obj", function(obj, target.gene, tfs,
                                                                 tf.weights=rep(1, length(tfs)), extraArgs=list())
                                                           standardGeneric ("solve"))
#----------------------------------------------------------------------------------------------------
#'
#' @name TReNA-class
#' @rdname TReNA-class

TReNA <- function(mtx.assay=matrix(), solver="lasso", quiet=TRUE)
{
    stopifnot(solver %in% c("lasso", "randomForest", "bayesSpike", "pearson",
                            "spearman","sqrtlasso","lassopv","ridge", "ensemble"))

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
    else if(solver == "ensemble")
        solver <- EnsembleSolver(mtx.assay)

    .TReNA(solver=solver, quiet=quiet)    

} # TReNA, the constructor
#----------------------------------------------------------------------------------------------------
#' Solve a TReNA Object
#'
#' A TReNA object contains an assay matrix with expression data for genes of interest and a string
#' representing the chosen solver. The \code{solve} method runs the specified solver given a target
#' gene and a designated set of transcription factors, returning a list of parameters that quantify
#' the relationship between the transcription factors and the target gene. 
#'
#' @name solve-methods
#' @rdname solve
#' @aliases solve solve.TReNA
#' 
#' @param obj An object of class TReNA
#' @param target.gene A designated target gene that should be part of the mtx.assay data
#' @param tfs The designated set of transcription factors that could be associated with the target gene.
#' @param tf.weights A set of weights on the transcription factors (default = rep(1, length(tfs)))
#' @param extraArgs Modifiers to the solver
#'
#' @return A data frame containing coefficients relating the target gene to each transcription factor
#'
#' @seealso \code{\link{TReNA}}
#'
#' @examples
#' # Load included Alzheimer's data, create a TReNA object with LASSO as the solver, and solve
#' load(system.file(package="TReNA", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
#' trena <- TReNA(mtx.assay = mtx.sub, solver = "lasso")
#' target.gene <- "APOE"
#' tfs <- setdiff(rownames(mtx.sub), target.gene)
#' tbl <- solve(trena, target.gene, tfs)

setMethod("solve", "TReNA",

   function (obj, target.gene, tfs, tf.weights=rep(1, length(tfs)), extraArgs=list()){
      # printf("entering TReNA::solve")
      run(obj@solver, target.gene, tfs, tf.weights, extraArgs)
      })
#----------------------------------------------------------------------------------------------------
