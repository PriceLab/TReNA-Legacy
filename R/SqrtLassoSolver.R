#------------------------------------------------------------------------------------------------------------------------
#' An S4 class to represent a Square Root LASSO solver
#'
#' @include Solver.R
#' @name SqrtLassoSolver-class

.SqrtLassoSolver <- setClass ("SqrtLassoSolver", contains="Solver")
#------------------------------------------------------------------------------------------------------------------------
#' Create a Solver class object using the Square Root LASSO solver
#'
#' @param mtx.assay An assay matrix of gene expression data
#'
#' @return A Solver class object with Square Root LASSO as the solver
#'
#' @examples
#' solver <- SqrtLassoSolver()

SqrtLassoSolver <- function(mtx.assay=matrix(), quiet=TRUE)
{
    obj <- .SqrtLassoSolver(Solver(mtx.assay=mtx.assay, quiet=quiet))

    obj

} # SqrtLassoSolver, the constructor
#------------------------------------------------------------------------------------------------------------------------
#' Get SqrtLasso Solver name
#'
#' @param obj An object of class SqrtLassoSolver
#' 
#' @return "SqrtLassoSolver"
#'
#' @examples
#' solver <- SqrtLassoSolver()
#' getSolverName(solver)

setMethod("getSolverName", "SqrtLassoSolver",

  function (obj){
     return("SqrtLassoSolver")
     })

#----------------------------------------------------------------------------------------------------
#' Run the Square Root LASSO Solver
#'
#' @aliases run.SqrtLassoSolver
#' @description Given a TReNA object with Square Root LASSO as the solver, use the \code{\link{slim}} function to estimate coefficients
#' for each transcription factor as a predictor of the target gene's expression level. 
#'
#' @param obj An object of class SqrtLassoSolver
#' @param target.gene A designated target gene that should be part of the mtx.assay data
#' @param tfs The designated set of transcription factors that could be associated with the target gene.
#' @param extraArgs Modifiers to the Square RootLASSO solver
#'
#' @return A data frame containing the coefficients relating the target gene to each transcription factor, plus other fit parameters.
#'
#' @seealso \code{\link{slim}}
#'
#' @examples
#' 

### Note: I've removed all references to alpha as I don't believe slim uses it
### Similar note for tf.weights and keep.metrics
setMethod("run", "SqrtLassoSolver",

          function (obj, target.gene, tfs, extraArgs=list()){
              
              if(length(tfs) == 0)                  
                  return(data.frame())              

              # Set default lambda and number of lambdas
              lambda <- NULL              
              nlambda <- 5
                    
              if("lambda" %in% names(extraArgs))                  
                  lambda <- extraArgs[["lambda"]]

              if("nlambda" %in% names(extraArgs))
                  nlambda <- 5
        
        # we don't try to handle tf self-regulation
                  deleters <- grep(target.gene, tfs)              
              if(length(deleters) > 0){                  
                  tfs <- tfs[-deleters]                  
                  tf.weights <- tf.weights[-deleters]                  
                  if(!obj@quiet)                   
                      message(sprintf("SqrtLassoSolver removing target.gene from candidate regulators: %s", target.gene))                  
              }
              
              if( length(tfs) == 0 ) return( data.frame() )
              
              mtx <- obj@mtx.assay              
              stopifnot(target.gene %in% rownames(mtx))             
              stopifnot(all(tfs %in% rownames(mtx)))              
              stopifnot(class(lambda) %in% c("NULL","numeric"))              
              features <- t(mtx[tfs,,drop=F ])              
              target <- as.numeric(mtx[target.gene,])
              
              if( length(tfs) == 1 ) {                  
                  fit = lm( target ~ features )                  
                  mtx.beta = coef(fit)                  
                  cor.target.feature = cor( target , features )[1,1]                 
                  mtx.beta = data.frame( beta = mtx.beta[2] , intercept = mtx.beta[1] , gene.cor = cor.target.feature )                  
                  rownames(mtx.beta) = tfs                  
                  return( mtx.beta )                  
              }
              
              # Run square root lasso and return an object of class "slim"              
              fit <- slim(features, target, method = "lq", lambda = lambda, verbose=FALSE)
              
              # extract the exponents of the fit             
              #tbl.out <- as.matrix(fit$beta)              
              #deleters <- as.integer(which(tbl.out[,1] == 0))              
              #if(length(deleters) > 0)              
              #   tbl.out <- tbl.out[-deleters, , drop=FALSE]              
              #colnames(tbl.out) <- "beta"
              
              # Pull out the coefficients from the final lambda        
              mtx.beta <- as.matrix(fit$beta[,nlambda])
              colnames(mtx.beta) <- "beta"
              rownames(mtx.beta) <- colnames(features)
              deleters <- as.integer(which(mtx.beta[,1] == 0))
              if( all( mtx.beta[,1] == 0 ) ) return( data.frame() )
              if(length(deleters) > 0)
                  mtx.beta <- mtx.beta[-deleters, , drop=FALSE]

              # put the intercept, admittedly with much redundancy, into its own column
              mtx.beta <- cbind(mtx.beta, intercept=rep(fit$intercept[[nlambda]], nrow(mtx.beta)))
              
              #browser()              
              correlations.of.betas.to.targetGene <- unlist(lapply(rownames(mtx.beta),
                                                                   function(x) cor(mtx[x,], mtx[target.gene,])))

              #browser()
              mtx.beta <- as.matrix(cbind( mtx.beta, gene.cor=correlations.of.betas.to.targetGene))
              if(!obj@quiet)
                  plot(fit$nlambda, label=TRUE)

              if( nrow(mtx.beta) > 1 ) {
                  ordered.indices <- order(abs(mtx.beta[, "beta"]), decreasing=TRUE)
                  mtx.beta <- mtx.beta[ordered.indices,]
              }

              return(as.data.frame(mtx.beta))
})
#----------------------------------------------------------------------------------------------------
#' Rescale Square Root LASSO Predictor Weights
#'
#' @aliases rescalePredictorWeights.SqrtLassoSolver
#'
#' @param obj An object of class SqrtLassoSolver
#' @param rawValue.min The minimum value of the raw expression values
#' @param rawValue.max The maximum value of the raw expression values
#' @param rawValues A matrix of raw expression values
#'
#' @return A matrix of the raw values re-scaled using the minimum and maximum values

# lasso penalizes predictors on a scale of 1 (full weight) to infinity (zero weight)
# here we wish to support incoming rawValues scaled between a possibly theoretical
# rawValue.min and rawValue.max
# we have empirical evidence that <large but non-infinite number> functions as a full penalty
# without distorting the scale so much that even good rawValues get reduced to nothing
# which is what .Machine$double.xmax would do

setMethod("rescalePredictorWeights", "SqrtLassoSolver",

   function (obj, rawValue.min, rawValue.max, rawValues){
      1 - ((rawValues-rawValue.min)/(rawValue.max-rawValue.min))
      })

#----------------------------------------------------------------------------------------------------

