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
#' @param extraArgs Modifiers to the Ensemble solver, including "solver.list", and "gene.cutoff"
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

               # Check if target.gene is in the bottom 10% in mean expression; if so, send a warning               
               if(rowMeans(obj@mtx.assay)[target.gene] < quantile(rowMeans(obj@mtx.assay), probs = 0.1)){                   
                   warning("Target gene mean expression is in the bottom 10% of all genes in the assay matrix")                   
               }               
               
               # Specify defaults for gene cutoff and solver list
               gene.cutoff <- 0.1
               solver.list <- "all.solvers"
               
               # Check for the gene cutoff and solvers, then and set them if they're not there
               if("gene.cutoff" %in% names(extraArgs))
                   gene.cutoff <- extraArgs[["gene.cutoff"]]

               if("solver.list" %in% names(extraArgs))
                   solver.list <- extraArgs[["solver.list"]]

               # Convert the "all" solvers argument
               if(solver.list[1] == "all.solvers"){
                   solver.list <- c("lasso",                                    
                                    "randomForest",                                    
#                                    "bayesSpike",                                    
                                    "pearson",                                    
                                    "spearman",                                    
#                                    "sqrtlasso",                                    
                                    "lassopv",
                                    "ridge")
               }
               out.list <- list()

               for(i in 1:length(solver.list)){
                   # Create and solve the TReNA object for each solver
                   trena <- TReNA(obj@mtx.assay, solver = solver.list[[i]] )

                   # if there's extraArgs, pass them
                   if(solver.list[[i]] %in% names(extraArgs)){                       
                       extraParams <- extraArgs[[solver.list[[i]]]]}                   
                   else{ extraParams <- list()}
                   
                      
                   out.list[[i]] <- solve(trena, target.gene, tfs, extraArgs = extraParams)
                   names(out.list)[i] <- paste("out",solver.list[[i]],sep=".")
               }

               # Output lasso with beta
               if("lasso" %in% solver.list){
                   out.list$out.lasso$gene <- rownames(out.list$out.lasso)                   
                   out.list$out.lasso <- out.list$out.lasso[, c("beta","gene")]                   
                   rownames(out.list$out.lasso) <- NULL
                   names(out.list$out.lasso) <- c("beta.lasso", "gene")
                   lasso.med <- median(out.list$out.lasso$beta.lasso)
                   lasso.scale <- mad(out.list$out.lasso$beta.lasso)
               }               

               # Output randomForest IncNodePurity
               if("randomForest" %in% solver.list){
                   out.list$out.randomForest <- out.list$out.randomForest$edges                   
                   out.list$out.randomForest$gene <- rownames(out.list$out.randomForest)
                   out.list$out.randomForest <- out.list$out.randomForest[, c("IncNodePurity","gene")]
                   rownames(out.list$out.randomForest) <- NULL
                   names(out.list$out.randomForest) <- c("rf.score", "gene")
                   randomForest.med <- median(out.list$out.randomForest$rf.score)
                   randomForest.scale <- sqrt(mean(
                       out.list$out.randomForest$rf.score*out.list$out.randomForest$rf.score))
               }

               # Output the z-score from bayesSpike
               if("bayesSpike" %in% solver.list){
                   out.list$out.bayesSpike$gene <- rownames(out.list$out.bayesSpike)
                   rownames(out.list$out.bayesSpike) <- NULL
                   out.list$out.bayesSpike <- out.list$out.bayesSpike[, c("z", "gene")]
                   names(out.list$out.bayesSpike) <- c("bayes.z", "gene")
                   bayesSpike.med <- median(out.list$out.bayesSpike$bayes.z)
                   bayesSpike.scale <- mad(out.list$out.bayesSpike$bayes.z)
               }

               # Pearson
               if("pearson" %in% solver.list){
                   out.list$out.pearson$gene <- rownames(out.list$out.pearson)                   
                   rownames(out.list$out.pearson) <- NULL
                   names(out.list$out.pearson) <- c("pearson.coeff","gene")
                   pearson.med <- median(out.list$out.pearson$pearson.coeff)
                   pearson.scale <- mad(out.list$out.pearson$pearson.coeff)
               }
               
               #Spearman
               if("spearman" %in% solver.list){
                   out.list$out.spearman$gene <- rownames(out.list$out.spearman)
                   rownames(out.list$out.spearman) <- NULL
                   names(out.list$out.spearman) <- c("spearman.coeff", "gene")
                   spearman.med <- median(out.list$out.spearman$spearman.coeff)
                   spearman.scale <- mad(out.list$out.spearman$spearman.coeff)
               }
               
               #LassoPV
               if("lassopv" %in% solver.list){
                   out.list$out.lassopv$gene <- rownames(out.list$out.lassopv)
                   rownames(out.list$out.lassopv) <- NULL
                   out.list$out.lassopv <- out.list$out.lassopv[, c("p.values","gene")]
                   names(out.list$out.lassopv) <- c("lasso.p.value", "gene")
                   p.log10 <- -log10(out.list$out.lassopv$lasso.p.value)
                   lassopv.med <- median(p.log10)
                   lassopv.scale <- sqrt(mean(p.log10*p.log10))
               }

               #SqrtLasso
               if("sqrtlasso" %in% solver.list){
                   out.list$out.sqrtlasso$gene <- rownames(out.list$out.sqrtlasso)
                   rownames(out.list$out.sqrtlasso) <- NULL
                   out.list$out.sqrtlasso <- out.list$out.sqrtlasso[, c("beta", "gene")]
                   names(out.list$out.sqrtlasso) <- c("beta.sqrtlasso", "gene")
                   sqrtlasso.med <- median(out.list$out.sqrtlasso$beta.sqrtlasso)
                   sqrtlasso.scale <- mad(out.list$out.sqrtlasso$beta.sqrtlasso)
               }
               
               #Ridge
               if("ridge" %in% solver.list){
                   out.list$out.ridge$gene <- rownames(out.list$out.ridge)
                   out.list$out.ridge <- out.list$out.ridge[, c("beta","gene")]
                   rownames(out.list$out.ridge) <- NULL
                   names(out.list$out.ridge) <- c("beta.ridge", "gene")
                   ridge.med <- median(out.list$out.ridge$beta.ridge)
                   ridge.scale <- mad(out.list$out.ridge$beta.ridge)
               }               

               # Grab the top "how.many" genes for each solver
               how.many <- round(length(tfs)*gene.cutoff)
               all.genes <- character(length = 2*how.many)
               top.list <- list()
                   
               while(length(all.genes) > gene.cutoff * length(tfs)){

                   # Grab the top "how.many" of each result
                   for(i in 1:length(out.list)){
                       top.list[[i]] <- head(out.list[[i]]$gene, how.many)
                   }
                   
                   all.genes <- unique(unlist(top.list))

                   how.many <- how.many - 1        
               }

               # Pull out the specified genes
               sub.list <- list()
               for(i in 1:length(out.list)){

                   sub.list[[i]] <- subset(out.list[[i]], out.list[[i]]$gene %in% all.genes)
               }
               
               # Merge the tables
               tbl.all <- merge(sub.list[[1]], sub.list[[2]], by = "gene", all = TRUE)
               if(length(sub.list) > 2){
                   for(i in 3:(length(sub.list))){                  
                       tbl.all <- merge(tbl.all, sub.list[[i]], by = "gene", all = TRUE)
                   }
               }
               
               # Replace missing values and scale the data
               # Use the *.med and *.scale values to center/scale everything
               tbl.all[is.na(tbl.all)] <- 0
               tbl.scale <- tbl.all[,-1]

               if("lasso.p.value" %in% names(tbl.scale)){                  
                   tbl.scale$lasso.p.value <- -log10(tbl.scale$lasso.p.value)
                   tbl.scale$lasso.p.value <- scale(tbl.scale$lasso.p.value,
                                                    center = lassopv.med,
                                                    scale = lassopv.scale)
               }

               if("beta.lasso" %in% names(tbl.scale)){
                   tbl.scale$beta.lasso <- scale(tbl.scale$beta.lasso,
                                                 center = lasso.med,
                                                 scale = lasso.scale)
               }
               
               if("beta.ridge" %in% names(tbl.scale)){
                   tbl.scale$beta.ridge <- scale(tbl.scale$beta.ridge,
                                                 center = ridge.med,
                                                 scale = ridge.scale)
               }
               
               if("pearson.coeff" %in% names(tbl.scale)){
                   tbl.scale$pearson.coeff <- scale(tbl.scale$pearson.coeff,
                                                    center = pearson.med,
                                                    scale = pearson.scale)
                   }

               if("spearman.coeff" %in% names(tbl.scale)){
                   tbl.scale$spearman.coeff <- scale(tbl.scale$spearman.coeff,
                                                     center = spearman.med,
                                                     scale = spearman.scale)
                   }

               if("beta.sqrtlasso" %in% names(tbl.scale)){
                   tbl.scale$beta.sqrtlasso <- scale(tbl.scale$beta.sqrtlasso,
                                                     center = sqrtlasso.med,
                                                     scale = sqrtlasso.scale)
                   }

               if("bayes.z" %in% names(tbl.scale)){
                   tbl.scale$bayes.z <- scale(tbl.scale$bayes.z,
                                              center = bayesSpike.med,
                                              scale = bayesSpike.scale)
                   }

               if("rf.score" %in% names(tbl.scale)){
                   tbl.scale$rf.score <- scale(tbl.scale$rf.score,
                                               center = randomForest.med,
                                               scale = randomForest.scale)
               }
               
#               browser()
#               tbl.scale <- scale(tbl.scale)
               rownames(tbl.scale) <- tbl.all$gene

               # Transform via PCA and compute the ensemble score
               pca <- prcomp(tbl.scale, center=FALSE, scale.=FALSE)
               extr <- apply(pca$x[,pca$sdev > 0.1],1, function(x) {sqrt(mean(x*x))})
               extr <- as.data.frame(extr)
               extr$gene <- rownames(extr)
               rownames(extr) <- NULL
               tbl.all <- merge(tbl.all, extr, by = "gene", all = TRUE)

               # Compute the scaled "composite score"
               comp <- apply(pca$x[, pca$sdev > 0.1], 1,
                             function(x) {sqrt(mean((2*atan(x)/pi)^2))})
               comp <- as.data.frame(comp)
               comp$gene <- rownames(comp)
               rownames(comp) <- NULL
               tbl.all <- merge(tbl.all, comp, by = "gene", all = TRUE)               
               tbl.all <- tbl.all[order(tbl.all$comp, decreasing = TRUE),]

               return(tbl.all)
           })
#----------------------------------------------------------------------------------------------------
