#------------------------------------------------------------------------------------------------------------------------
#' An S4 class to represent a LASSO P-Value solver
#'
#' @include Solver.R
#' @import lassopv
#' @import methods
#' 
#' @name LassoPVSolver-class

.LassoPVSolver <- setClass ("LassoPVSolver", contains="Solver")
#------------------------------------------------------------------------------------------------------------------------
#' Create a Solver class object using the LASSO P-Value solver
#'
#' @param mtx.assay An assay matrix of gene expression data
#' @param target.gene A designated target gene that should be part of the mtx.assay data
#' @param candidateRegulators The designated set of transcription factors that could be associated
#' with the target gene
#' @param quiet A logical denoting whether or not the solver should print output
#' 
#' @return A Solver class object with LASSO P-Value as the solver
#'
#' @seealso  \code{\link{solve.LassoPV}}, \code{\link{getAssayData}}
#'
#' @family Solver class objects
#' 
#' @export
#' 
#' @examples
#' solver <- LassoPVSolver()

LassoPVSolver <- function(mtx.assay=matrix(), targetGene, candidateRegulators, quiet=TRUE)
{
    # Remove the targetGene from candidateRegulators
    if(any(grepl(targetGene, candidateRegulators)))        
        candidateRegulators <- candidateRegulators[-grep(targetGene, candidateRegulators)]    

    # Check to make sure the matrix contains some of the candidates
    candidateRegulators <- intersect(candidateRegulators, rownames(mtx.assay))    
    stopifnot(length(candidateRegulators) > 0)
    
    obj <- .LassoPVSolver(Solver(mtx.assay=mtx.assay,
                                 quiet=quiet,
                                 targetGene = targetGene,
                                 candidateRegulators = candidateRegulators))

    # Send a warning if there's a row of zeros
    if(!is.na(max(mtx.assay)) & any(rowSums(mtx.assay) == 0))
       warning("One or more gene has zero expression; this may cause problems when using P-Value LASSO. You may want to try 'lasso' or 'ridge' instead.")

    obj

} # LassoPVSolver, the constructor
#----------------------------------------------------------------------------------------------------
setMethod('show', 'LassoPVSolver',

    function(obj) {
       regulator.count <- length(getRegulators(obj))
       if(regulator.count > 10){
          regulatorString <- paste(getRegulators(obj)[1:10], collapse=",")
          regulatorString <- sprintf("%s...", regulatorString);
          }
       else
          regulatorString <- paste(getRegulators(obj), collapse=",")

       msg = sprintf("LassoPVSolver with mtx.assay (%d, %d), targetGene %s, %d candidate regulators %s",
                     nrow(getAssayData(obj)), ncol(getAssayData(obj)),
                     getTarget(obj), regulator.count, regulatorString)
       cat (msg, '\n', sep='')
    })
#----------------------------------------------------------------------------------------------------
#' Run the LASSO P-Value Solver
#'
#' @rdname solve.LassoPV
#' @aliases run.LassoPVSolver solve.LassoPV
#' 
#' @description Given a TReNA object with LASSO P-Value as the solver, use the \code{\link{lassopv}}
#' function to estimate coefficients for each transcription factor as a predictor of the target
#' gene's expression level.
#' 
#' @param obj An object of class LassoPVSolver
#'
#' @return A data frame containing the p-values for each transcription factor pertaining to the target gene
#' plus the Pearson correlations between each transcription factor and the target gene.
#'
#' @seealso \code{\link{lassopv}}, , \code{\link{LassoPVSolver}}
#'
#' @family solver methods
#' 
#' @examples
#' # Load included Alzheimer's data, create a TReNA object with Bayes Spike as solver, and solve
#' load(system.file(package="TReNA", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
#' target.gene <- "MEF2C"
#' tfs <- setdiff(rownames(mtx.sub), target.gene)
#' lassopv.solver <- LassoPVSolver(mtx.sub, target.gene, tfs)
#' tbl <- run(lassopv.solver)

setMethod("run", "LassoPVSolver",

          function (obj){
              
              mtx <- getAssayData(obj)
              target.gene <- getTarget(obj)
              tfs <- getRegulators(obj)             

              # Check if target.gene is in the bottom 10% in mean expression; if so, send a warning         
              if(rowMeans(mtx)[target.gene] < stats::quantile(rowMeans(mtx), probs = 0.1)){                 
                  warning("Target gene mean expression is in the bottom 10% of all genes in the assay matrix")            
              }
              
              # we don't try to handle tf self-regulation           
              deleters <- grep(target.gene, tfs)
              if(length(deleters) > 0){                  
                  tfs <- tfs[-deleters]                  
                  if(!obj@quiet)                   
                      message(sprintf("LassoPVSolver removing target.gene from candidate regulators: %s", target.gene))                  
              }                          
              
              if(length(tfs) == 0)                  
                  return(data.frame()) 
              stopifnot(target.gene %in% rownames(mtx))             
              stopifnot(all(tfs %in% rownames(mtx)))              
              features <- t(mtx[tfs,,drop=FALSE ])              
              target <- as.numeric(mtx[target.gene,])

              # Run LASSO P-Value and return the P-values, ordered by increasing value
              fit <- lassopv(features, target)
              fit <- fit[order(fit, decreasing=FALSE)]

              # Add pearson correlations and make a data frame              
              correlations.of.betas.to.targetGene <- unlist(lapply(names(fit),
                                                                   function(x) stats::cor(mtx[x,], mtx[target.gene,])))
              tbl <- data.frame(row.names = names(fit),
                                p.values = fit,
                                gene.cor=correlations.of.betas.to.targetGene)
              return(tbl)
})
#----------------------------------------------------------------------------------------------------
