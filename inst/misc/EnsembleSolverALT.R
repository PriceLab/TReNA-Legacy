#----------------------------------------------------------------------------------------------------
#' Class EnsembleSolver
#'
#' @include Solver.R
#' @import methods
#'
#' @name EnsembleSolver-class
#'
.EnsembleSolver <- setClass("EnsembleSolver",
                            slots = c(solverList = "list",
                                      geneCutoff = "numeric",
                                      quiet = "logical")
                            )
#----------------------------------------------------------------------------------------------------
#' Create a Solver class object using an ensemble of solvers
#'
#' @param solverList A list object containing objects of the Solver class
#' @param geneCutoff A decimal indicating the fraction of genes to use from each solver object
#' @param quiet A logical denoting whether or not the solver should print output
#'
#' @return A Solver class object with Ensemble as the solver
#'
#' @seealso  \code{\link{solve.Ensemble}}, \code{\link{getAssayData}}
#'
#' @family Solver class objects
#'
#' @export
#'
#' @examples
#' solver <- EnsembleSolver()

EnsembleSolver <- function(solverList,
                           geneCutoff = 0.1,
                           quiet=TRUE)
{ 
    obj <- .EnsembleSolver(solverList = solverList
                           geneCutoff = geneCutoff,               
                           quiet=quiet)
    obj

} # EnsembleSolver, the constructor
#----------------------------------------------------------------------------------------------------
#' Run the Ensemble Solver
#'
#' @name run,EnsembleSolver-method
#' @rdname solve.Ensemble
#' @aliases run.EnsembleSolver solve.Ensemble
#'
#' @description Given a TReNA object with Ensemble as the solver and a list of solvers
#' (default = "default.solvers"), estimate coefficients for each transcription factor
#' as a predictor of the target gene's expression level. The final scores for the ensemble
#' method combine all specified solvers to create a composite score for each transcription factor.
#' This method should be called using the \code{\link{solve}} method on an appropriate TReNA object.
#'
#' @param obj An object of class EnsembleSolver
#'
#' @return A data frame containing the scores for all solvers and two composite scores
#' relating the target gene to each transcription factor. The two new scores are:
#' \itemize{
#' \item{"concordance": a composite score created similarly to "extreme_score", but with each solver's
#' score scaled using *atan(x)*. This score scales from 0-1}
#' \item{"pcaMax": a composite score created using the root mean square of the principal
#' components of the individual solver scores}
#' }
#'
#' @seealso \code{\link{EnsembleSolver}}
#' 
#' @family solver methods
#'
#' @examples
#' # Load included Alzheimer's data, create 3 solver objects, and solve
#' load(system.file(package="TReNA", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
#' target.gene <- "MEF2C"
#' tfs <- setdiff(rownames(mtx.sub), target.gene)
#' lasso.solver <- LassoSolver(mtx.sub, target.gene, tfs)
#' rf.solver <- RandomForestSolver(mtx.sub, target.gene, tfs)
#' bayes.solver <- BayesSpikeSolver(mtx.sub, target.gene, tfs)
#' ensemble.solver <- EnsembleSolver(list(lasso.solver, rf.solver, bayes.solver))
#' tbl <- run(ensemble.solver)
#'

setMethod("run", "EnsembleSolver",

          function(obj){

              solver.list <- obj@solverList
              mtx <- getAssayData(solver.list[[1]])
              target.gene <- getTarget(solver.list[[1]])
              tfs <- getRegulators(solver.list[[1]])
          
               # Check if target.gene is in the bottom 10% in mean expression; if so, send a warning
               if(rowMeans(mtx)[target.gene] < stats::quantile(rowMeans(mtx), probs = 0.1)){
                   warning("Target gene mean expression is in the bottom 10% of all genes in the assay matrix")
               }

              # Check to make sure matrices are the same
              

              
              # Check to make sure targets are all the same; send a warning if not
              all.targets <- vapply(solver.list, getTarget, FUN.VALUE = "character")
              if(!all(all.targets == target.gene)){
                  warning("One or more of the supplied target genes does not match")
              }
              
              # Check to make sure regulators are all the same; send a warning if not
              all.tfs <- sapply(solver.list, getRegulators)
              tf.compare <- sapply(all.tfs, setequal, tfs)
              if(!all(tf.compare)){
                  warning("One or more of the supplied candidate regulator sets do not match")
              }

              # Pick out the solvers present
              all.solvers <- sapply(solver.list, class)

              # Weed out any duplicate solvers

              
              out.list <- list()
              
              solver <- switch(requested.solver.root,
                     "lasso" = LassoSolver(mtx.assay, solverName),
                     "randomForest" = RandomForestSolver(mtx.assay, solverName),
                     "bayesSpike" = BayesSpikeSolver(mtx.assay, solverName),
                     "pearson" = PearsonSolver(mtx.assay, solverName),
                     "spearman" = SpearmanSolver(mtx.assay, solverName),
                     "sqrtlasso" = SqrtLassoSolver(mtx.assay, solverName),
                     "lassopv" = LassoPVSolver(mtx.assay, solverName),
                     "ridge" = RidgeSolver(mtx.assay, solverName),
                     "ensemble" = EnsembleSolver(mtx.assay, solverName))

              for(i in 1:length(solver.list)){                  
                  # Create and solve the TReNA object for each solver                  
                  trena <- TReNA(getAssayData(obj), solver = solver.list[[i]] )                  
                  out.list[[i]] <- solve(trena, target.gene, tfs, tf.weights, extraArgs = extraParams)                  
                  names(out.list)[i] <- paste("out",solver.list[[i]],sep=".")                  
              }
              
               # Output lasso with beta
               if("lasso" %in% solver.list){
                   out.list$out.lasso$gene <- rownames(out.list$out.lasso)
                   out.list$out.lasso <- out.list$out.lasso[, c("beta","gene")]
                   rownames(out.list$out.lasso) <- NULL
                   names(out.list$out.lasso) <- c("beta.lasso", "gene")
                   lasso.med <- stats::median(out.list$out.lasso$beta.lasso)
                   lasso.scale <- stats::mad(out.list$out.lasso$beta.lasso)
               }

               # Output randomForest IncNodePurity
               if("randomForest" %in% solver.list){
                   out.list$out.randomForest <- out.list$out.randomForest$edges
                   out.list$out.randomForest$gene <- rownames(out.list$out.randomForest)
                   out.list$out.randomForest <- out.list$out.randomForest[, c("IncNodePurity","gene")]
                   rownames(out.list$out.randomForest) <- NULL
                   names(out.list$out.randomForest) <- c("rf.score", "gene")
                   randomForest.med <- stats::median(out.list$out.randomForest$rf.score)
                   randomForest.scale <- sqrt(mean(
                       out.list$out.randomForest$rf.score*out.list$out.randomForest$rf.score))
               }

               # Output the z-score from bayesSpike
               if("bayesSpike" %in% solver.list){
                   out.list$out.bayesSpike$gene <- rownames(out.list$out.bayesSpike)
                   rownames(out.list$out.bayesSpike) <- NULL
                   out.list$out.bayesSpike <- out.list$out.bayesSpike[, c("z", "gene")]
                   names(out.list$out.bayesSpike) <- c("bayes.z", "gene")
                   bayesSpike.med <- stats::median(out.list$out.bayesSpike$bayes.z)
                   bayesSpike.scale <- stats::mad(out.list$out.bayesSpike$bayes.z)
               }

               # Pearson
               if("pearson" %in% solver.list){
                   out.list$out.pearson$gene <- rownames(out.list$out.pearson)
                   rownames(out.list$out.pearson) <- NULL
                   names(out.list$out.pearson) <- c("pearson.coeff","gene")
                   pearson.med <- stats::median(out.list$out.pearson$pearson.coeff)
                   pearson.scale <- stats::mad(out.list$out.pearson$pearson.coeff)
               }

               #Spearman
               if("spearman" %in% solver.list){
                   out.list$out.spearman$gene <- rownames(out.list$out.spearman)
                   rownames(out.list$out.spearman) <- NULL
                   names(out.list$out.spearman) <- c("spearman.coeff", "gene")
                   spearman.med <- stats::median(out.list$out.spearman$spearman.coeff)
                   spearman.scale <- stats::mad(out.list$out.spearman$spearman.coeff)
               }

               #LassoPV
               if("lassopv" %in% solver.list){
                   out.list$out.lassopv$gene <- rownames(out.list$out.lassopv)
                   rownames(out.list$out.lassopv) <- NULL
                   out.list$out.lassopv <- out.list$out.lassopv[, c("p.values","gene")]
                   names(out.list$out.lassopv) <- c("lasso.p.value", "gene")
                   p.log10 <- -log10(out.list$out.lassopv$lasso.p.value)
                   lassopv.med <- stats::median(p.log10)
                   lassopv.scale <- sqrt(mean(p.log10*p.log10))
               }

               #SqrtLasso
               if("sqrtlasso" %in% solver.list){
                   out.list$out.sqrtlasso$gene <- rownames(out.list$out.sqrtlasso)
                   rownames(out.list$out.sqrtlasso) <- NULL
                   out.list$out.sqrtlasso <- out.list$out.sqrtlasso[, c("beta", "gene")]
                   names(out.list$out.sqrtlasso) <- c("beta.sqrtlasso", "gene")
                   sqrtlasso.med <- stats::median(out.list$out.sqrtlasso$beta.sqrtlasso)
                   sqrtlasso.scale <- stats::mad(out.list$out.sqrtlasso$beta.sqrtlasso)
               }

               #Ridge
               if("ridge" %in% solver.list){
                   out.list$out.ridge$gene <- rownames(out.list$out.ridge)
                   out.list$out.ridge <- out.list$out.ridge[, c("beta","gene")]
                   rownames(out.list$out.ridge) <- NULL
                   names(out.list$out.ridge) <- c("beta.ridge", "gene")
                   ridge.med <- stats::median(out.list$out.ridge$beta.ridge)
                   ridge.scale <- stats::mad(out.list$out.ridge$beta.ridge)
               }

               # Grab all genes for each solver to start with
               how.many <- length(tfs)
               top.list <- lapply(out.list, function(x) head(x$gene, how.many))
               all.genes <- unique(as.character(unlist(top.list)))

               # Run same thing in a while loop until the cutoff or 10 is reached
               while(length(all.genes) > gene.cutoff * length(tfs) & length(all.genes) > 10){

                   how.many <- round(0.75*how.many)
                   top.list <- lapply(out.list, function(x) utils::head(x$gene, how.many))
                   all.genes <- unique(as.character(unlist(top.list)))
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

               rownames(tbl.scale) <- tbl.all$gene

               # Compute the scaled "concordance score"
               pca <- stats::prcomp(tbl.scale, center=FALSE, scale.=FALSE)
               pca$x <- pca$x / sqrt(length(which(pca$sdev > 0.1)))
               concordance <- apply(pca$x[, pca$sdev > 0.1], 1,
                             function(x) {sqrt(mean((2*atan(x)/pi)^2))})
               concordance <- as.data.frame(concordance)
               concordance$gene <- rownames(concordance)
               rownames(concordance) <- NULL
               tbl.all <- merge(tbl.all, concordance, by = "gene", all = TRUE)

               # Transform via PCA and compute the pcaMax score
               pcaMax <- apply(pca$x[,pca$sdev > 0.1],1, function(x) {sqrt(mean(x*x))})
               pcaMax <- as.data.frame(pcaMax)
               pcaMax$gene <- rownames(pcaMax)
               rownames(pcaMax) <- NULL
               tbl.all <- merge(tbl.all, pcaMax, by = "gene", all = TRUE)

               # Sort by pcaMax
               tbl.all <- tbl.all[order(tbl.all$pcaMax, decreasing = TRUE),]

               return(tbl.all)
           })
#----------------------------------------------------------------------------------------------------
