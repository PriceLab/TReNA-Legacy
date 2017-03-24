#----------------------------------------------------------------------------------------------------
#' Run the Elastic Net Solvers
#' @description Given a TReNA object with either LASSO or Ridge Regression as the solver, use the \code{\link{glmnet}} function to estimate coefficients for each transcription factor as a predictor of the target gene's expression level.
#'
#' @param obj An object of class Solver
#' @param target.gene A designated target gene that should be part of the mtx.assay data
#' @param tfs The designated set of transcription factors that could be associated with the target gene.
#' @param tf.weights A set of weights on the transcription factors (default = rep(1, length(tfs)))
#' @param alpha The LASSO/Ridge tuning parameter
#' @param lambda The penalty tuning parameter for elastic net
#' @param keep.metrics A binary variable indicating whether or not to keep metrics
#'
#' @return A data frame containing the coefficients relating the target gene to each transcription factor, plus other fit parameters
#'
#' @seealso \code{\link{glmnet}}
#'
#' @examples
#' 

.elasticNetSolver <-  function (obj, target.gene, tfs, tf.weights, alpha, lambda, keep.metrics){

   if(length(tfs) == 0)
       return(data.frame())
    
        # we don't try to handle tf self-regulation
    deleters <- grep(target.gene, tfs)
    if(length(deleters) > 0){
       tfs <- tfs[-deleters]
       tf.weights <- tf.weights[-deleters]
     if(!obj@quiet)
	message(sprintf("Removing target.gene from candidate regulators: %s", target.gene))
       }

     if( length(tfs) == 0 ) return( data.frame() )

     mtx <- obj@mtx.assay
     stopifnot(target.gene %in% rownames(mtx))
     stopifnot(all(tfs %in% rownames(mtx)))
     stopifnot(class(lambda) %in% c("NULL","numeric"))
     features <- t(mtx[tfs,,drop=F ])
     target <- as.numeric(mtx[target.gene,])

     if( length(tfs) == 1 ) {
       fit = stats::lm( target ~ features )
       mtx.beta = stats::coef(fit)
       cor.target.feature = stats::cor( target , features )[1,1]
       mtx.beta = data.frame( beta = mtx.beta[2] , intercept = mtx.beta[1] , gene.cor = cor.target.feature )
       rownames(mtx.beta) = tfs
       if( keep.metrics == F ) return( mtx.beta )
       if( keep.metrics == T ) return( list( mtx.beta = mtx.beta , lambda = NA , r2 = cor.target.feature^2 ) )
     }

     if( is.null(lambda) ) {

         # Run Permutation testing to find lambda
         if( alpha != 0 )
             alpha.perm = alpha
         else(alpha.perm = 0.1)
         target.mixed <- sample(target)
         threshold <- 1E-15
         lambda.change <- 10^(-4)
         lambda <- 1
         lambda.list <- numeric(length=50)
         
         for(i in 1:length(lambda.list)){             
             # Do a binary search             
             step.size <- lambda/2 # Start at 0.5             
             while(step.size > lambda.change){
                 
                 # Get the fit
                 fit <- glmnet(features, target.mixed, penalty.factor = tf.weights, alpha=alpha.perm, lambda=lambda)
                 
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
             lambda.list[[i]] <- lambda
         }
         # Give lambda as 1 + 1se
         lambda <- mean(lambda.list) + (stats::sd(lambda.list)/sqrt(length(lambda.list)))
         
         fit <- glmnet(features, target, penalty.factor=tf.weights, alpha=alpha, lambda=lambda)
         
         }

         # For non-LASSO
 #        else{
 #            fit <- cv.glmnet(features, target, penalty.factor=tf.weights, grouped=FALSE , alpha = alpha )         
 #            lambda.min <- fit$lambda.min         
 #            lambda <-fit$lambda.1se          
 #        }

     else if(is.numeric(lambda)){
         fit = glmnet(features, target, penalty.factor=tf.weights, alpha=alpha, lambda=lambda)
         }

    # extract the exponents of the fit
    mtx.beta <- as.matrix( stats::predict( fit , newx = features , type = "coef" , s = lambda ) )
    colnames(mtx.beta) <- "beta"
    deleters <- as.integer(which(mtx.beta[,1] == 0))
    if( all( mtx.beta[,1] == 0 ) ) return( data.frame() )
    if(length(deleters) > 0)
        mtx.beta <- mtx.beta[-deleters, , drop=FALSE]

    # put the intercept, admittedly with much redundancy, into its own column    
    intercept <- mtx.beta[1,1]    
    mtx.beta <- mtx.beta[-1, , drop=FALSE]    
    mtx.beta <- cbind(mtx.beta, intercept=rep(intercept, nrow(mtx.beta)))    
    correlations.of.betas.to.targetGene <- unlist(lapply(rownames(mtx.beta), function(x) stats::cor(mtx[x,], mtx[target.gene,])))
    

     mtx.beta <- as.matrix(cbind( mtx.beta, gene.cor=correlations.of.betas.to.targetGene))
     #if(!obj@quiet)
     #   graphics::plot(fit.nolambda, xvar='lambda', label=TRUE)

     if( nrow(mtx.beta) > 1 ) {
        ordered.indices <- order(abs(mtx.beta[, "beta"]), decreasing=TRUE)
        mtx.beta <- mtx.beta[ordered.indices,]
     }

     mtx.beta = as.data.frame(mtx.beta)

     if( keep.metrics == TRUE ) {
        pred.values = predict( fit , newx = features , s = lambda , type = "link" )
        r2 = (cor( target , pred.values )[1,1])^2
        return( list( mtx.beta = mtx.beta , lambda = lambda , r2 = r2 ) )
     }

     if( keep.metrics == FALSE )
        return(as.data.frame(mtx.beta))
     }
# ElasticNetSolver 
#----------------------------------------------------------------------------------------------------
