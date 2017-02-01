# This is a simple Spearman solver
.SpearmanSolver <- setClass ("SpearmanSolver", contains = "Solver")
#----------------------------------------------------------------------------------------------------
SpearmanSolver <- function(mtx.assay = matrix(), quiet=TRUE)
{
    obj <- .SpearmanSolver(Solver(mtx.assay=mtx.assay, quiet=quiet))

    obj

}
#----------------------------------------------------------------------------------------------------
setMethod("getSolverName", "SpearmanSolver",

          function (obj){
              return("SpearmanSolver")
          })

#----------------------------------------------------------------------------------------------------
setMethod("run", "SpearmanSolver",

          function (obj, target.gene, tfs, tf.weights = rep(1,length(tfs), extraArgs=list())){

              mtx <- obj@mtx.assay
              # Check that target gene and tfs are all part of the matrix
              stopifnot(target.gene %in% rownames(mtx))
              stopifnot(all(tfs %in% rownames(mtx)))
              # If given no tfs, return nothing
              if (length(tfs)==0) return(NULL)

              # Don't handle tf self-regulation, so take target gene out of tfs
              deleters <- grep(target.gene, tfs)
              if(length(deleters) > 0){
                  tfs <- tfs[-deleters]
                  tf.weights = tf.weights[-deleters]
              }
              # If target gene was the only tf, then return nothing
              if(length(tfs)==0) return(NULL)

              x = t(mtx[tfs,,drop=F])
              y = as.vector(t(mtx[target.gene,])) # Make target gene levels into a vector

              # Calculate Spearman correlation coefficients
              fit <- cor( x = x, y = y, method = "spearman")

                                        # For now, just return the fit
              return(fit)
})