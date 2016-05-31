#------------------------------------------------------------------------------------------------------------------------
.RandomForestSolver <- setClass ("RandomForestSolver", contains="Solver")
#------------------------------------------------------------------------------------------------------------------------
RandomForestSolver <- function(mtx.assay=matrix(), quiet=TRUE)
{
    obj <- .RandomForestSolver(Solver(mtx.assay=mtx.assay, quiet=quiet))

    obj

} # TReNA, the constructor
#------------------------------------------------------------------------------------------------------------------------
setMethod("getSolverName", "RandomForestSolver",

  function (obj){
     return("RandomForestSolver")
     })

#----------------------------------------------------------------------------------------------------
setMethod("run", "RandomForestSolver",

  function (obj, target.gene, tfs, tf.weights=rep(1,length(tfs), extraArgs=list())){

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
     y = t(mtx[target.gene,])

     fit <- randomForest( x = x, y = y )
     edges = as.data.frame(fit$importance)
     pred.values = predict(fit)
     r2 = cor( pred.values , mtx[target.gene,])^2
     return( list( edges = edges , r2 = r2 ) )
     })


#----------------------------------------------------------------------------------------------------
