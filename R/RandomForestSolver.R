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

  function (obj, target.gene, tfs, tf.weights=rep(1,length(tfs))){

     mtx <- obj@mtx.assay
     stopifnot(target.gene %in% rownames(mtx))
     stopifnot(all(tfs %in% rownames(mtx)))
     result <- randomForest(t(mtx[tfs,]), t(mtx[target.gene,]))
     return(list(scores=result$importance, rSquared=NA))
     })


#----------------------------------------------------------------------------------------------------
