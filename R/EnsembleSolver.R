#----------------------------------------------------------------------------------------------------
#' An S4 class to represent an Ensemble Solver
#'
#' @include Solver.R
#' @name EnsembleSolver-class
#'
.EnsembleSolver <- setClass("EnsembleSolver", contains="Solver")
#----------------------------------------------------------------------------------------------------
#' Create a Solver class object using an ensemble approach
#'
#' @param mtx.assay An assay matrix of gene expression data
#' @param quiet A logical denoting whether or not the solver should print output
#' @return A Solver class object with Ensemble as the solver
#'
#' @examples
#' solver <- EnsembleSolver()
#'
#' @family solver methods

EnsembleSolver <- function(mtx.assay=matrix(), quiet=TRUE)
{
    obj <- .EnsembleSolver(Solver(mtx.assay=mtx.assay, quiet=quiet))

    obj

} # EnsembleSolver, the constructor
#----------------------------------------------------------------------------------------------------
#' Get Ensemble Solver name
#'
#' @param obj An object of class EnsembleSolver
#'
#' @return "EnsembleSolver"
#'
#' @examples
#' solver <- EnsembleSolver()
#' getSolverName(solver)

setMethod("getSolverName", "EnsembleSolver",

          function(obj){
              return("EnsembleSolver")
          })

#----------------------------------------------------------------------------------------------------
#' Run the Ensemble Solver
#'
#' @rdname EnsembleSolver
#' @aliases run.EnsembleSolver
#' @description Given a TReNA object with Ensemble as the solver and a list of solvers (default = all solvers), estimate coefficients for each transcription factor as a predictor of the target gene's expression level. The final scores for the ensemble method combine all specified solvers to create a composite score for each transcription factor. 
#' @param obj An object of class EnsembleSolver
#' @param target.gene A designated target gene that should be part of the mtx.assay data
#' @param tfs The designated set of transcription factors that could be associated with the target gene
#' @param extraArgs Modifiers to the Ensemble solver
#'
#' @return A data frame containing the scores for all solvers and a composite score relating the target gene to each transcription factor
#'
#' @examples
#' # Load included Alzheimer's data, create a TReNA object with LASSO as solver, and solve
#' load(system.file(package="TReNA", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
#' trena <- TReNA(mtx.assay = mtx.sub, solver = "ensemble")
#' target.gene <- "APOE"
#' tfs <- setdiff(rownames(mtx.sub), target.gene)
#' tbl <- solve(trena, target.gene, tfs)

setMethods("run", "EnsembleSolver",

           function(obj, target.gene, tfs, extraArgs = list()){
               # Specify the % of total genes you want
               gene.cutoff <- 0.1

               # Check for the gene cutoff and set it if it's there
               if("gene.cutoff" %in% extraArgs)
                   gene.cutoff <- extraArgs[["gene.cutoff"]]
               
               # For now, run the Lasso and RF solvers, and that's it
               solver.list <- c("lasso","randomForest")
               out.list <- list(length = length(solver.list))

               browser()
               for(i in 1:length(solver.list)){

                   # Create and solve the TReNA object
                   trena <- TReNA(obj@mtx.assay, solver = solver.list[[i]] )
                   out.list[[i]] <- solve(trena, target.gene, tfs)
                   names(out.list)[i] <- paste("out",solver.list[[i]],sep=".")
                   
               }

               # Drop the gene.cor from lasso
               out.lasso <- out.list$out.lasso
               out.lasso$gene <- rownames(out.lasso)
               out.lasso <- out.lasso[, c("beta","gene")]
               rownames(out.lasso) <- NULL

               out.randomForest <- out.list$out.randomForest$edges
               out.randomForest$gene <- rownames(out.randomForest)
               rownames(out.randomForest) <- NULL


               # Grab the top "how.many" genes for each solver
               how.many <- round(length(tfs)*gene.cutoff)
               all.genes <- character(length = 2*how.many)
               while(length(all.genes) > gene.cutoff * length(tfs)){ # Needs some work

                   out.lasso.top <- head(out.lasso$gene, how.many)
                   out.randomForest.top <- head(out.randomForest$gene, how.many)

                   all.genes <- unique(c(out.lasso.top,
                                         out.randomForest.top))

                   how.many <- how.many - 1        
               }

               # Pull out the specified genes
               out.lasso.sub <- subset(out.lasso, out.lasso$gene %in% all.genes)
               out.randomForest.sub <- subset(out.randomForest, out.randomForest$gene %in% all.genes)
               
               # Merge the tables
               tbl.all <- merge(out.lasso.sub, out.randomForest.sub, by = "gene", all = TRUE)

               # Replace missing values and scale the data
               tbl.all[is.na(tbl.all)] <- 0
               tbl.scale <- scale(tbl.all[,-1])
               rownames(tbl.scale) <- tbl.all$gene

               # Transform via PCA and compute the ensemble score
               pca <- prcomp(tbl.scale, center=FALSE, scale.=FALSE)
               extr <- apply(pca$x[,pca$sdev > 0.1],1, function(x) {sqrt(sum(x*x))})
               extr <- as.data.frame(extr)
               extr$gene <- rownames(extr)
               rownames(extr) <- NULL
               tbl.all <- merge(tbl.all, extr, by = "gene", all = TRUE)
               tbl.all <- tbl.all[order(tbl.all$extr, decreasing = TRUE),]

               return(tbl.all)

               
               
               
           })
#----------------------------------------------------------------------------------------------------
