#------------------------------------------------------------------------------------------------------------------------
#' An S4 class to represent a Square Root LASSO solver
#'
#' @import flare
#' @import doParallel
#' @import foreach
#' @import methods
#' 
#' @include Solver.R
#' 
#' @name SqrtLassoSolver-class

.SqrtLassoSolver <- setClass ("SqrtLassoSolver",
                              contains="Solver",
                              slots = c(regulatorWeights="numeric",
                                        lambda = "numeric",
                                        nCores = "numeric")
                              )
#----------------------------------------------------------------------------------------------------
#' Create a Solver class object using the Square Root LASSO solver
#'
#' @param mtx.assay An assay matrix of gene expression data
#' @param target.gene A designated target gene that should be part of the mtx.assay data
#' @param candidateRegulators The designated set of transcription factors that could be associated
#' with the target gene
#' @param regulatorWeights A set of weights on the transcription factors
#' (default = rep(1, length(tfs)))
#' @param lambda A tuning parameter that determines the severity of the penalty function imposed
#' on the elastic net regression. If unspecified, lambda will be determined via
#' permutation testing (default = numeric(0)).
#' @param nCores An integer specifying the number of computational cores to devote to this
#' square root LASSO solver. This solver is generally quite slow and is greatly sped up when using
#' multiple cores (default = 4)
#' @param quiet A logical denoting whether or not the solver should print output
#' 
#' @return A Solver class object with Square Root LASSO as the solver
#'
#' @seealso  \code{\link{solve.SqrtLasso}}, \code{\link{getAssayData}}
#'
#' @family Solver class objects
#' 
#' @export
#' 
#' @examples
#' solver <- SqrtLassoSolver()

SqrtLassoSolver <- function(mtx.assay=matrix(), targetGene, candidateRegulators,
                            regulatorWeights = rep(1, length(candidateRegulators)),
                            lambda = numeric(0), nCores = 4, quiet=TRUE)
{
    if(any(grepl(targetGene, candidateRegulators)))
      candidateRegulators <- candidateRegulators[-grep(targetGene, candidateRegulators)]

    candidateRegulators <- intersect(candidateRegulators, rownames(mtx.assay))
    
    stopifnot(length(candidateRegulators) > 0)    

    obj <- .SqrtLassoSolver(mtx.assay=mtx.assay,                        
                        quiet=quiet,                        
                        targetGene=targetGene,                        
                        candidateRegulators=candidateRegulators,                        
                        regulatorWeights=regulatorWeights,                        
                        lambda = lambda,                        
                        nCores = nCores                        
                       )

    # Send a warning if there's a row of zeros
    if(!is.na(max(mtx.assay)) & any(rowSums(mtx.assay) == 0))
       warning("One or more gene has zero expression; this may cause problems when using Square Root LASSO. You may want to try 'lasso' or 'ridge' instead.")

    obj

} # SqrtLassoSolver, the constructor
#----------------------------------------------------------------------------------------------------
setMethod('show', 'SqrtLassoSolver',

    function(obj) {
       regulator.count <- length(getRegulators(obj))
       if(regulator.count > 10){
          regulatorString <- paste(getRegulators(obj)[1:10], collapse=",")
          regulatorString <- sprintf("%s...", regulatorString);
          }
       else
          regulatorString <- paste(getRegulators(obj), collapse=",")

       msg = sprintf("SqrtLassoSolver with mtx.assay (%d, %d), targetGene %s, %d candidate regulators %s,  with %d cores",
                     nrow(getAssayData(obj)), ncol(getAssayData(obj)),
                     getTarget(obj), regulator.count, regulatorString, obj@nCores)
       cat (msg, '\n', sep='')
    })
#----------------------------------------------------------------------------------------------------
#' Run the Square Root LASSO Solver
#'
#' @rdname solve.SqrtLasso
#' @aliases run.SqrtLassoSolver solve.SqrtLasso
#' 
#' @description Given a TReNA object with Square Root LASSO as the solver,
#' use the \code{\link{slim}} function to estimate coefficients
#' for each transcription factor as a predictor of the target gene's expression level.
#' This method should be called using the \code{\link{solve}} method on an appropriate TReNA object.
#' 
#' @param obj An object of class Solver with "sqrtlasso" as the solver string
#'
#' @return A data frame containing the coefficients relating the target gene to
#' each transcription factor, plus other fit parameters.
#'
#' @seealso \code{\link{slim}}, \code{\link{SqrtLassoSolver}}
#'
#' @family solver methods
#'
#' @examples
#' # Load included Alzheimer's data, create a TReNA object with Square Root LASSO as solver, and solve
#' load(system.file(package="TReNA", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
#' target.gene <- "MEF2C"
#' tfs <- setdiff(rownames(mtx.sub), target.gene)
#' sqrt.solver <- SqrtLassoSolver(mtx.sub, target.gene, tfs)
#' tbl <- run(sqrt.solver)
#'
#' # Solve the same problem but use 8 cores
#' sqrt.solver <- SqrtLassoSolver(mtx.sub, target.gene, tfs, nCores = 8)
#' tbl <- run(sqrt.solver)

setMethod("run", "SqrtLassoSolver",

          function (obj){

              mtx <- getAssayData(obj)
              target.gene <- getTarget(obj)
              tfs <- getRegulators(obj)
              lambda <- obj@lambda
              nCores <- obj@nCores
              
              if(length(tfs) == 0)                  
                  return(data.frame())              

              # Check if target.gene is in the bottom 10% in mean expression; if so, send a warning
              if(rowMeans(mtx)[target.gene] < stats::quantile(rowMeans(mtx), probs = 0.1)){
                  warning("Target gene mean expression is in the bottom 10% of all genes in the assay matrix")
                  }

              # we don't try to handle tf self-regulation              
              deleters <- grep(target.gene, tfs)              
              if(length(deleters) > 0){                  
                  tfs <- tfs[-deleters]                  
                  if(!obj@quiet)
                      message(sprintf("SqrtLassoSolver removing target.gene from candidate regulators: %s", target.gene))
              }
              
              if( length(tfs) == 0 ) return( data.frame() )
                           
              stopifnot(target.gene %in% rownames(mtx))             
              stopifnot(all(tfs %in% rownames(mtx)))              
              stopifnot(class(lambda) %in% c("NULL","numeric"))              
              features <- t(mtx[tfs,,drop=FALSE ])              
              target <- as.numeric(mtx[target.gene,])
              
              if( length(tfs) == 1 ) {                  
                  fit = stats::lm( target ~ features )                  
                  mtx.beta = stats::coef(fit)                  
                  cor.target.feature = stats::cor( target , features )[1,1]                 
                  mtx.beta = data.frame( beta = mtx.beta[2] , intercept = mtx.beta[1] , gene.cor = cor.target.feature )                  
                  rownames(mtx.beta) = tfs                  
                  return( mtx.beta )                  
              }

              # If no lambda, run a binary search for the best lasso using permutation of the data set
              if(length(lambda) == 0){

                  target.mixed <- sample(target)
                  threshold <- 1E-15
                  lambda.change <- 10^(-4)
                  lambda <- 1
                  
                  cl <- parallel::makeCluster(nCores)
                  doParallel::registerDoParallel(cl)
                  
                  lambda.list <- foreach::foreach(i = 1:30) %dopar% {

                      # Do a binary search
                      step.size <- lambda/2 # Start at 0.5
                      while(step.size > lambda.change){
                          # Get the fit
                          fit <- slim(features, target.mixed, method = "lq", verbose = FALSE, lambda = lambda)
                          # Case 1: nonsense, need to lower lambda
                          if(max(fit$beta) < threshold){
                              lambda <- lambda - step.size
                          }
                          # Case 2: sense, need to raise lambda
                          else{
                              lambda <- lambda + step.size
                          }
                          # Halve the step size and re-scramble the target
                          step.size <- step.size/2
                          target.mixed <- sample(target)
                      }
                      lambda
                  }
                  # Grab the lambdas and average them
                  lambda.list <- unlist(lambda.list)                  
                  parallel::stopCluster(cl)                  
                  lambda <- mean(lambda.list) + (stats::sd(lambda.list)/sqrt(length(lambda.list)))                                 
              }              

              # Run square root lasso and return an object of class "slim"              
              fit <- slim(features, target, method = "lq", lambda = lambda, verbose=FALSE)
              
              # Pull out the coefficients        
              mtx.beta <- as.matrix(fit$beta)
              colnames(mtx.beta) <- "beta"
              rownames(mtx.beta) <- colnames(features)
              deleters <- as.integer(which(mtx.beta[,1] == 0))
              if( all( mtx.beta[,1] == 0 ) ) return( data.frame() )
              if(length(deleters) > 0)
                  mtx.beta <- mtx.beta[-deleters, , drop=FALSE]

              # put the intercept, admittedly with much redundancy, into its own column
              mtx.beta <- cbind(mtx.beta, intercept=rep(fit$intercept, nrow(mtx.beta)))
              
              #browser()              
              correlations.of.betas.to.targetGene <- unlist(lapply(rownames(mtx.beta),
                                                                   function(x) stats::cor(mtx[x,], mtx[target.gene,])))

              mtx.beta <- as.matrix(cbind( mtx.beta, gene.cor=correlations.of.betas.to.targetGene))
#              if(!obj@quiet)
#                  plot(fit$nlambda, label=TRUE)

              if( nrow(mtx.beta) > 1 ) {
                  ordered.indices <- order(abs(mtx.beta[, "beta"]), decreasing=TRUE)
                  mtx.beta <- mtx.beta[ordered.indices,]
              }

              return(as.data.frame(mtx.beta))
})
#----------------------------------------------------------------------------------------------------
