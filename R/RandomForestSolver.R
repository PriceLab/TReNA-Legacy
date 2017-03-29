#' Class RandomForestSolver
#'
#' @import randomForest
#' @include Solver.R
#' @import methods
#' 
#' @name RandomForestSolver-class
#' @rdname RandomForestSolver-class

.RandomForestSolver <- setClass ("RandomForestSolver", contains="Solver")
#------------------------------------------------------------------------------------------------------------------------
#' Create a Solver class object using the Random Forest solver
#' 
#' @param mtx.assay An assay matrix of gene expression data
#' @param quiet A logical denoting whether or not the solver should print output
#' @return A Solver class object with Random Forest as the solver
#'
#' @export
#' 
#' @examples
#' solver <- RandomForestSolver()


RandomForestSolver <- function(mtx.assay=matrix(), quiet=TRUE)
{
    obj <- .RandomForestSolver(Solver(mtx.assay=mtx.assay, quiet=quiet))

    # Send a warning if there's a row of zeros
    if(!is.na(max(mtx.assay)) & any(rowSums(mtx.assay) == 0))
       warning("One or more gene has zero expression; this may yield warnings when using Random Forest.")

    obj

} # RandomForestSolver, the constructor
#------------------------------------------------------------------------------------------------------------------------
#' Get RandomForest Solver Name
#'
#' @param obj An object of class RandomForestSolver
#' 
#' @return "RandomForestSolver"
#'
#' @export
#' 
#' @examples
#' solver <- RandomForestSolver()
#' getSolverName(solver)

setMethod("getSolverName", "RandomForestSolver",

  function (obj){
     return("RandomForestSolver")
     })

#----------------------------------------------------------------------------------------------------
#' Run the Random Forest Solver
#'
#' @rdname solve.RandomForest
#' @aliases run.RandomForestSolver solve.RandomForest
#'
#' @param obj An object of class TReNA with "randomForest" as the solver string
#' @param target.gene A designated target gene that should be part of the mtx.assay data
#' @param tfs The designated set of transcription factors that could be associated with the target gene.
#' @param tf.weights A set of weights on the transcription factors (default = rep(1, length(tfs)))
#' @param extraArgs Modifiers to the Random Forest solver
#'
#' @return A list containing various parameters of the Random Forest fit.
#'
#' @seealso \code{\link{randomForest}}
#' 
#' @family solver methods
#' 
#' @examples
#' # Load included Alzheimer's data, create a TReNA object with Random Forest as solver, and solve
#' load(system.file(package="TReNA", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
#' trena <- TReNA(mtx.assay = mtx.sub, solver = "randomForest")
#' target.gene <- "MEF2C"
#' tfs <- setdiff(rownames(mtx.sub), target.gene)
#' tbl <- solve(trena, target.gene, tfs)


setMethod("run", "RandomForestSolver",

  function (obj, target.gene, tfs, tf.weights=rep(1,length(tfs), extraArgs=list())){

      # Check if target.gene is in the bottom 10% in mean expression; if so, send a warning      
      if(rowMeans(obj@mtx.assay)[target.gene] < stats::quantile(rowMeans(obj@mtx.assay), probs = 0.1)){          
          warning("Target gene mean expression is in the bottom 10% of all genes in the assay matrix")          
      }      
      
     mtx <- obj@mtx.assay
     stopifnot(target.gene %in% rownames(mtx))
     stopifnot(all(tfs %in% rownames(mtx)))
     if(length(tfs)==0) return(NULL)

        # we don't try to handle tf self-regulation
     deleters <- grep(target.gene, tfs)
     if(length(deleters) > 0){
       tfs <- tfs[-deleters]
       tf.weights = tf.weights[-deleters]
     }
     if(length(tfs)==0) return(NULL)

     x = t(mtx[tfs,,drop=F])
     y = as.vector(t(mtx[target.gene,])) # Change y to a vector to avoid RF warning

     fit <- randomForest( x = x, y = y )
     edges = as.data.frame(fit$importance)
     pred.values = stats::predict(fit)
     r2 = stats::cor(pred.values , mtx[target.gene,])^2
     gene.cor <- sapply(rownames(edges), function(tf) stats::cor(mtx[tf,], mtx[target.gene,]))
     edges$gene.cor <- gene.cor
     edges <- edges[order(edges$IncNodePurity, decreasing=TRUE),]
     return(list(edges = edges , r2 = r2))
     })

#----------------------------------------------------------------------------------------------------
