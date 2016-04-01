#------------------------------------------------------------------------------------------------------------------------
.LassoSolver <- setClass ("LassoSolver", contains="Solver")
#------------------------------------------------------------------------------------------------------------------------
LassoSolver <- function(mtx.assay=matrix(), quiet=TRUE)
{
    obj <- .LassoSolver(Solver(mtx.assay=mtx.assay, quiet=quiet))

    obj

} # TReNA, the constructor
#------------------------------------------------------------------------------------------------------------------------
setMethod("getSolverName", "LassoSolver",

  function (obj){
     return("LassoSolver")
     })

#----------------------------------------------------------------------------------------------------
setMethod("run", "LassoSolver",

  function (obj, target.gene, tfs){

     mtx <- obj@mtx.assay
     stopifnot(target.gene %in% rownames(mtx))
     stopifnot(all(tfs %in% rownames(mtx)))
     features <- t(mtx[tfs, ])
     target <- t(mtx[target.gene, , drop=FALSE])

     if(!obj@quiet)
         printf("begining cross-validation for glmnet, using %d tfs, target %s", length(tfs), target.gene)

     cv.out <- cv.glmnet(features, target, grouped=FALSE)
     lambda.min <- cv.out$lambda.min

     if(!obj@quiet) printf("glmnet fit with lambda %f", lambda.min)
     fit = glmnet(features, target, lambda=lambda.min)

       # extract the exponents of the fit
     mtx.beta <- as.matrix(coef(fit, s="lambda.min"))
     deleters <- as.integer(which(mtx.beta[,1] == 0))
     if(length(deleters) > 0)
        mtx.beta <- mtx.beta[-deleters, , drop=FALSE]

     return(mtx.beta)
     # for sake of speed, no longer calculatining r-squared here
     # the test it provides can now be accomplished via trainModel and predictFromModel
     #predicted <- predict.glmnet(fit, features)
     #SSE <- (predicted-target)^2;
     #SST <- (target-mean(target))^2
     #R.squared <- 1-(SSE/SST)
     #return(list(scores=mtx.beta, rSquared=R.squared))
     })


#----------------------------------------------------------------------------------------------------
setMethod("trainModel", "LassoSolver",

   function (obj, target.gene, tfs, training.samples){

      mtx <- obj@mtx.assay
      stopifnot(target.gene %in% rownames(mtx))
      stopifnot(all(tfs %in% rownames(mtx)))
      features <- t(mtx[tfs, ])
      target <- t(mtx[target.gene, , drop=FALSE])

      if(!obj@quiet)
          printf("begining cross-validation for glmnet, using %d tfs, target %s", length(tfs), target.gene)

      cv.out <- cv.glmnet(features, target, grouped=FALSE)
      lambda.min <- cv.out$lambda.min
      glmnet(features, target, lambda=lambda.min)
      })

#----------------------------------------------------------------------------------------------------
setMethod("predictFromModel", "LassoSolver",

   function (obj, model, tfs, test.samples){
      mtx <- obj@mtx.assay
      predict.glmnet(model, t(mtx[tfs, test.samples]))
      })

#----------------------------------------------------------------------------------------------------

