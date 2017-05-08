#' Class RandomForestSolver
#'
#' @import randomForest
#' @include Solver.R
#' @import methods
#'
#' @name RandomForestSolver-class
#' @rdname RandomForestSolver-class

.RandomForestSolver <- setClass ("RandomForestSolver",
                                 contains="Solver",
                                 slots=c(targetGene="character",
                                         candidateRegulators="character",
                                         regulatorWeights="numeric")
                                 )
#----------------------------------------------------------------------------------------------------
#' Create a Solver class object using the Random Forest solver
#'
#' @param mtx.assay An assay matrix of gene expression data
#' @param quiet A logical denoting whether or not the solver should print output
#'
#' @return A Solver class object with Random Forest as the solver
#'
#' @seealso  \code{\link{solve.RandomForest}}, \code{\link{getAssayData}}
#'
#' @family Solver class objects
#' 
#' @export
#'
#' @examples
#' solver <- RandomForestSolver()

RandomForestSolver <- function(mtx.assay=matrix(), targetGene, candidateRegulators,
                               regulatorWeights=rep(1, length(candidateRegulators)),
                               quiet=TRUE)
{
   if(any(grepl(targetGene, candidateRegulators)))
      candidateRegulators <- candidateRegulators[-grep(targetGene, candidateRegulators)]

   candidateRegulators <- intersect(candidateRegulators, rownames(mtx.assay))
   stopifnot(length(candidateRegulators) > 0)

   obj <- .RandomForestSolver(Solver(mtx.assay=mtx.assay, quiet=quiet),
                              targetGene=targetGene,
                              candidateRegulators=candidateRegulators,
                              regulatorWeights=regulatorWeights)

    # Send a warning if there's a row of zeros
    if(!is.na(max(mtx.assay)) & any(rowSums(mtx.assay) == 0))
       warning("One or more gene has zero expression; this may yield warnings when using Random Forest.")

    obj

} # RandomForestSolver, the constructor
#----------------------------------------------------------------------------------------------------
setMethod('show', 'RandomForestSolver',

    function(object) {
       regulator.count <- length(object@candidateRegulators)
       if(regulator.count > 10){
          regulatorString <- paste(object@candidateRegulators[1:10], collapse=",")
          regulatorString <- sprintf("%s...", regulatorString);
          }
       else
          regulatorString <- paste(object@candidateRegulators, collapse=",")

       msg = sprintf("RandomForestSolver with mtx.assay (%d, %d), targetGene %s, %d candidate regulators %s",
                     nrow(object@mtx.assay), ncol(object@mtx.assay),
                     object@targetGene, regulator.count, regulatorString)
       cat (msg, '\n', sep='')
       })

#----------------------------------------------------------------------------------------------------
#' Run the Random Forest Solver
#'
#' @rdname solve.RandomForest
#' @aliases run.RandomForestSolver solve.RandomForest
#'
#' @description
#' Given a TReNA object with RandomForest as the solver, use the \code{\link{randomForest}} function
#' to estimate coefficients for each transcription factor as a predictor of the target gene's
#' expression level.
#' This method should be called using the \code{\link{solve}} method on an appropriate TReNA object.

#' @param obj An object of class TReNA with "randomForest" as the solver string
#' @param targetGene A designated target gene that should be part of the mtx.assay data
#' @param candidateRegulators The designated set of transcription factors that could be associated with the target gene.
#' @param tf.weights A set of weights on the transcription factors (default = rep(1, length(candidateRegulators)))
#' @param extraArgs Modifiers to the Random Forest solver
#'
#' @return A list containing various parameters of the Random Forest fit.
#'
#' @seealso \code{\link{randomForest}}, \code{\link{RandomForestSolver}}
#' 
#' @family solver methods
#'
#' @examples
#' # Load included Alzheimer's data, create a TReNA object with Random Forest as solver, and solve
#' load(system.file(package="TReNA", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
#' trena <- TReNA(mtx.assay = mtx.sub, solver = "randomForest")
#' targetGene <- "MEF2C"
#' candidateRegulators <- setdiff(rownames(mtx.sub), targetGene)
#' tbl <- solve(trena, targetGene, candidateRegulators)


setMethod("run", "RandomForestSolver",

  function (obj){

      # Check if target.gene is in the bottom 10% in mean expression; if so, send a warning      
      if(rowMeans(getAssayData(obj))[target.gene] < stats::quantile(rowMeans(getAssayData(obj)), probs = 0.1)){          
          warning("Target gene mean expression is in the bottom 10% of all genes in the assay matrix")          
      }      
      
     mtx <- getAssayData(obj)
     stopifnot(target.gene %in% rownames(mtx))
     stopifnot(all(tfs %in% rownames(mtx)))
     if(length(tfs)==0) return(NULL)

     mtx <- obj@mtx.assay
     if(length(obj@candidateRegulators)==0) return(NULL)

     x = t(mtx[obj@candidateRegulators,,drop=F])
     y = as.vector(t(mtx[obj@targetGene,])) # Change y to a vector to avoid RF warning

     fit <- randomForest( x = x, y = y )
     edges = as.data.frame(fit$importance)
     pred.values = stats::predict(fit)
     r2 = stats::cor(pred.values , mtx[obj@targetGene,])^2
     gene.cor <- sapply(rownames(edges), function(tf) stats::cor(mtx[tf,], mtx[obj@targetGene,]))
     edges$gene.cor <- gene.cor
     edges <- edges[order(edges$IncNodePurity, decreasing=TRUE),]
     return(list(edges = edges , r2 = r2))
     })
#----------------------------------------------------------------------------------------------------
