#------------------------------------------------------------------------------------------------------------------------
.TReNA <- setClass ("TReNA",
                    representation = representation(mtx.assay="matrix",
                                                    quiet="logical")
                    )

#------------------------------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#------------------------------------------------------------------------------------------------------------------------
setGeneric("getAssayData",             signature="obj", function(obj) standardGeneric ("getAssayData"))
setGeneric("fit",                      signature="obj", function(obj, target.gene, tfs) standardGeneric ("fit"))
setGeneric("getModel",                 signature="obj", function(obj) standardGeneric ("getModel"))
#------------------------------------------------------------------------------------------------------------------------
TReNA <- function(mtx.assay=matrix(), quiet=TRUE)
{
  .TReNA(mtx.assay=mtx.assay, quiet=quiet)

} # TReNA, the constructor
#------------------------------------------------------------------------------------------------------------------------
setMethod("getAssayData", "TReNA",

  function (obj){
     invisible(obj@mtx.assay)
     })

#------------------------------------------------------------------------------------------------------------------------
setMethod("fit", "TReNA",

  function (obj, target.gene, tfs){

     mtx <- getAssayData(obj)
     tfs.of.interest <- intersect(tfs, colnames(mtx))
     features <- mtx[, tfs.of.interest]
     stopifnot(target.gene %in% colnames(mtx))

     target <- mtx[, target.gene]

     if(!obj@quiet)
         printf("begining cross-validation for glmnet, using %d tfs, target %s", length(tfs), target.gene)

     cv.out <- cv.glmnet(features, target, grouped=FALSE)
     lambda.min <- cv.out$lambda.min

     if(!obj@quiet) printf("glmnet fit with labmda %f", lambda.min)
     fit = glmnet(features, target, lambda=lambda.min)

       # extract the exponents of the fit
     mtx.beta <- as.matrix(coef(cv.out, s="lambda.min"))
     deleters <- as.integer(which(mtx.beta[,1] == 0))
     if(length(deleters) > 0)
        mtx.beta <- mtx.beta[-deleters, , drop=FALSE]
     return(mtx.beta)
     })

#------------------------------------------------------------------------------------------------------------------------
