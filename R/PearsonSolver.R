#----------------------------------------------------------------------------------------------------
#' An S4 class to represent a Pearson solver
#'
#' @include Solver.R
#' @import methods
#' 
#' @name PearsonSolver-class
#' 

.PearsonSolver <- setClass ("PearsonSolver",contains = "Solver")
#----------------------------------------------------------------------------------------------------
#' Create a Solver class object using  Pearson correlation coefficients as the solver
#'
#' @param mtx.assay An assay matrix of gene expression data
#' @param target.gene A designated target gene that should be part of the mtx.assay data
#' @param candidateRegulators The designated set of transcription factors that could be associated
#' with the target gene
#' @param quiet A logical denoting whether or not the solver should print output
#' 
#' @return A Solver class object with Pearson correlation coefficients as the solver
#'
#' @seealso  \code{\link{solve.Pearson}}, \code{\link{getAssayData}}
#'
#' @family Solver class objects
#' 
#' @export
#' 
#' @examples
#' solver <- PearsonSolver()

PearsonSolver <- function(mtx.assay = matrix(), targetGene, candidateRegulators, quiet=TRUE)
{
    # Remove the targetGene from candidateRegulators
    if(any(grepl(targetGene, candidateRegulators)))        
        candidateRegulators <- candidateRegulators[-grep(targetGene, candidateRegulators)]    

    # Check to make sure the matrix contains some of the candidates
    candidateRegulators <- intersect(candidateRegulators, rownames(mtx.assay))    
    stopifnot(length(candidateRegulators) > 0)
    
    obj <- .PearsonSolver(Solver(mtx.assay=mtx.assay,
                                 quiet=quiet,
                                 targetGene=targetGene,                                 
                                 candidateRegulators=candidateRegulators))                          

    # Send a warning if there's a row of zeros
    if(!is.na(max(mtx.assay)) & any(rowSums(mtx.assay) == 0))
       warning("One or more gene has zero expression; this may yield 'NA' results and warnings when using Pearson correlations")

    obj

} #PearsonSolver, the constructor
#----------------------------------------------------------------------------------------------------
setMethod('show', 'PearsonSolver',

    function(obj) {
       regulator.count <- length(getRegulators(obj))
       if(regulator.count > 10){
          regulatorString <- paste(getRegulators(obj)[1:10], collapse=",")
          regulatorString <- sprintf("%s...", regulatorString);
          }
       else
          regulatorString <- paste(getRegulators(obj), collapse=",")

       msg = sprintf("PearsonSolver with mtx.assay (%d, %d), targetGene %s, %d candidate regulators %s",
                     nrow(getAssayData(obj)), ncol(getAssayData(obj)),
                     getTarget(obj), regulator.count, regulatorString)
       cat (msg, '\n', sep='')
    })
#----------------------------------------------------------------------------------------------------
#' Run the Pearson Solver
#'
#' @rdname solve.Pearson
#' @aliases run.PearsonSolver solve.Pearson
#'
#' @description Given a TReNA object with Pearson as the solver, use the \code{\link{cor}} function
#' to estimate coefficients for each transcription factor as a perdictor of the target gene's
#' expression level. 
#'
#' @param obj An object of class PearsonSolver
#' 
#' @return The set of Pearson Correlation Coefficients between each transcription factor and the target gene.
#'
#' @seealso \code{\link{cor}}, \code{\link{PearsonSolver}}
#'
#' @family solver methods
#' 
#' @examples
#' # Load included Alzheimer's data, create a TReNA object with Bayes Spike as solver, and solve
#' load(system.file(package="TReNA", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
#' target.gene <- "MEF2C"
#' tfs <- setdiff(rownames(mtx.sub), target.gene)
#' pearson.solver <- PearsonSolver(mtx.sub, target.gene, tfs)
#' tbl <- run(pearson.solver)

setMethod("run", "PearsonSolver",

          function (obj){

              # Check if target.gene is in the bottom 10% in mean expression; if so, send a warning              
              if(rowMeans(getAssayData(obj))[getTarget(obj)] < stats::quantile(rowMeans(getAssayData(obj)), probs = 0.1)){                  
                  warning("Target gene mean expression is in the bottom 10% of all genes in the assay matrix")                  
              }
              
              mtx <- getAssayData(obj)
              target.gene <- getTarget(obj)
              tfs <- getRegulators(obj)
              
              # Check that target gene and tfs are all part of the matrix
              stopifnot(target.gene %in% rownames(mtx))
              stopifnot(all(tfs %in% rownames(mtx)))
              # If given no tfs, return nothing
              if (length(tfs)==0) return(NULL)

              # Don't handle tf self-regulation, so take target gene out of tfs
              deleters <- grep(target.gene, tfs)
              if(length(deleters) > 0){
                  tfs <- tfs[-deleters]
              }
              # If target gene was the only tf, then return nothing
              if(length(tfs)==0) return(NULL)

              x = t(mtx[tfs,,drop=FALSE])
              y = as.vector(t(mtx[target.gene,])) # Make target gene levels into a vector

              # Calculate Pearson correlation coefficients
              fit <- stats::cor( x = x, y = y)

              # Return the coefficients as a data frame 
              tbl <- data.frame(row.names = rownames(fit)[order(abs(fit), decreasing = TRUE)],
                                coefficient = fit[order(abs(fit), decreasing = TRUE)])

              return(tbl)
          })
#----------------------------------------------------------------------------------------------------
