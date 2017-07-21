#----------------------------------------------------------------------------------------------------
#' Class EnsembleSolver
#'
#' @include Solver.R
#' @import methods
#'
#' @name EnsembleSolver-class
#'
.EnsembleSolver <- setClass("EnsembleSolver",
                            contains="Solver",
                            slots = c(solverNames = "character",
                                      geneCutoff = "numeric",
                                      alpha.lasso = "numeric",
                                      alpha.ridge = "numeric",
                                      lambda.lasso = "numeric",
                                      lambda.ridge = "numeric",
                                      lambda.sqrt = "numeric",
                                      nCores.sqrt = "numeric",
                                      nOrderings.bayes = "numeric"
                                      )
                            )

#----------------------------------------------------------------------------------------------------
#' Create a Solver class object using an ensemble of solvers
#'
#' @param mtx.assay An assay matrix of gene expression data
#' @param targetGene A designated target gene that should be part of the mtx.assay data
#' @param candidateRegulators The designated set of transcription factors that could be associated
#' with the target gene
#' @param solverNames A character vector of strings denoting 
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

EnsembleSolver <- function(mtx.assay=matrix(), targetGene, candidateRegulators,
                           solverNames = c("lasso",
                                           "lassopv",
                                           "pearson",
                                           "randomForest",
                                           "ridge",
                                           "spearman"),
                           geneCutoff = 0.1,
                           alpha.lasso = 0.9,
                           alpha.ridge = 0.0,
                           lambda.lasso = numeric(0),
                           lambda.ridge = numeric(0),
                           lambda.sqrt = numeric(0),
                           nCores.sqrt = 4,
                           nOrderings.bayes = 10,
                           quiet=TRUE)
{
    if(any(grepl(targetGene, candidateRegulators)))        
        candidateRegulators <- candidateRegulators[-grep(targetGene, candidateRegulators)]    

    candidateRegulators <- intersect(candidateRegulators, rownames(mtx.assay))    
    stopifnot(length(candidateRegulators) > 0)    

    # Send a warning if there's a row of zeros
    if(!is.na(max(mtx.assay)) & any(rowSums(mtx.assay) == 0))
        warning("One or more gene has zero expression; this may cause difficulty when using Bayes Spike. You may want to try 'lasso' or 'ridge' instead.")

    obj <- .EnsembleSolver(mtx.assay = mtx.assay,
                           targetGene = targetGene,
                           candidateRegulators = candidateRegulators,
                           solverNames = solverNames,
                           geneCutoff = geneCutoff,
                           alpha.lasso = alpha.lasso,
                           alpha.ridge = alpha.ridge,
                           lambda.lasso = lambda.lasso,
                           lambda.ridge = lambda.ridge,
                           lambda.sqrt = lambda.sqrt,
                           nCores.sqrt = nCores.sqrt,
                           nOrderings.bayes = nOrderings.bayes,
                           quiet=quiet)
    obj

} # EnsembleSolver, the constructor
#----------------------------------------------------------------------------------------------------
setMethod('show', 'EnsembleSolver',

    function(obj) {
       regulator.count <- length(getRegulators(obj))
       if(regulator.count > 10){
          regulatorString <- paste(getRegulators(obj)[1:10], collapse=",")
          regulatorString <- sprintf("%s...", regulatorString);
          }
       else
          regulatorString <- paste(getRegulators(obj), collapse=",")

       msg = sprintf("EnsembleSolver with mtx.assay (%d, %d), targetGene %s, %d candidate regulators %s,  solvers: %s",
                     nrow(getAssayData(obj)), ncol(getAssayData(obj)),
                     getTarget(obj), regulator.count, regulatorString,
                     paste(obj@solverNames,collapse = ", "))
       cat (msg, '\n', sep='')
    })
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
#' @param obj An object of class Solver with "ensemble" as the solver string
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
#' # Load included Alzheimer's data, create an Ensemble object with default solvers, and solve
#' load(system.file(package="TReNA", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
#' target.gene <- "MEF2C"
#' tfs <- setdiff(rownames(mtx.sub), target.gene)
#' ensemble.solver <- EnsembleSolver(mtx.sub, target.gene, tfs)
#' tbl <- run(ensemble.solver)
#'
#' # Solve the same problem, but supply extra arguments that change alpha for LASSO to 0.8 and also
#' # Change the gene cutoff from 10% to 20%
#' ensemble.solver <- EnsembleSolve(mtx.sub, target.gene, tfs, gene.cutoff = 0.2, alpha.lasso = 0.8)))
#' tbl <- run(ensemble.solver)
#'
#' # Solve the original problem with default cutoff and solver parameters, but use only 4 solvers
#' ensemble.solver <- EnsembleSolver(mtx.sub, target.gene, tfs,
#' solverNames = c("lasso", "randomForest", "pearson", "bayesSpike"))
#' tbl <- run(ensemble.solver)

setMethod("run", "EnsembleSolver",

          function(obj){

              mtx <- getAssayData(obj)
              target.gene <- getTarget(obj)
              tfs <- getRegulators(obj)
              gene.cutoff <- obj@geneCutoff

               # Check if target.gene is in the bottom 10% in mean expression; if so, send a warning
               if(rowMeans(mtx)[target.gene] < stats::quantile(rowMeans(mtx), probs = 0.1)){
                   warning("Target gene mean expression is in the bottom 10% of all genes in the assay matrix")
               }

              # Create a list of solvers and a list for solutions
              out.list <- list()              
              solver.list <- obj@solverNames

              for(i in 1:length(solver.list)){
                  # Find the correct solver
                  solver <- switch(tolower(solver.list[i]),
                                   "lasso" = LassoSolver(mtx, target.gene, tfs,
                                                         alpha = obj@alpha.lasso, lambda = obj@lambda.lasso),
                                   "randomforest" = RandomForestSolver(mtx, target.gene, tfs),
                                   "bayesspike" = BayesSpikeSolver(mtx, target.gene, tfs,
                                                                   nOrderings = obj@nOrderings.bayes),           
                                   "pearson" = PearsonSolver(mtx, target.gene, tfs),                                   
                                   "spearman" = SpearmanSolver(mtx, target.gene, tfs),                                   
                                   "sqrtlasso" = SqrtLassoSolver(mtx, target.gene, tfs,
                                                                 lambda = obj@lambda.sqrt, nCores = obj@nCores.sqrt),             
                                   "lassopv" = LassoPVSolver(mtx, target.gene, tfs),                           
                                   "ridge" = RidgeSolver(mtx, target.gene, tfs,
                                                         alpha = obj@alpha.ridge, lambda = obj@lambda.ridge))
                  
                  # Solve each Solver object and save it to the output list                  
                  out.list[[i]] <- run(solver)
                  names(out.list)[i] <- paste("out",tolower(solver.list[[i]]),sep=".")                  
              }

               # Output lasso with beta
               if("lasso" %in% tolower(solver.list)){
                   out.list$out.lasso$gene <- rownames(out.list$out.lasso)
                   out.list$out.lasso <- out.list$out.lasso[, c("beta","gene")]
                   rownames(out.list$out.lasso) <- NULL
                   names(out.list$out.lasso) <- c("beta.lasso", "gene")
                   lasso.med <- stats::median(out.list$out.lasso$beta.lasso)
                   lasso.scale <- stats::mad(out.list$out.lasso$beta.lasso)
               }

               # Output randomforest IncNodePurity
               if("randomforest" %in% tolower(solver.list)){
                   out.list$out.randomforest <- out.list$out.randomforest$edges
                   out.list$out.randomforest$gene <- rownames(out.list$out.randomforest)
                   out.list$out.randomforest <- out.list$out.randomforest[, c("IncNodePurity","gene")]
                   rownames(out.list$out.randomforest) <- NULL
                   names(out.list$out.randomforest) <- c("rf.score", "gene")
                   randomforest.med <- stats::median(out.list$out.randomforest$rf.score)
                   randomforest.scale <- sqrt(mean(
                       out.list$out.randomforest$rf.score*out.list$out.randomforest$rf.score))
               }

               # Output the z-score from bayesspike
               if("bayesspike" %in% tolower(solver.list)){
                   out.list$out.bayesspike$gene <- rownames(out.list$out.bayesspike)
                   rownames(out.list$out.bayesspike) <- NULL
                   out.list$out.bayesspike <- out.list$out.bayesspike[, c("z", "gene")]
                   names(out.list$out.bayesspike) <- c("bayes.z", "gene")
                   bayesspike.med <- stats::median(out.list$out.bayesspike$bayes.z)
                   bayesspike.scale <- stats::mad(out.list$out.bayesspike$bayes.z)
               }

               # Pearson
               if("pearson" %in% tolower(solver.list)){
                   out.list$out.pearson$gene <- rownames(out.list$out.pearson)
                   rownames(out.list$out.pearson) <- NULL
                   names(out.list$out.pearson) <- c("pearson.coeff","gene")
                   pearson.med <- stats::median(out.list$out.pearson$pearson.coeff)
                   pearson.scale <- stats::mad(out.list$out.pearson$pearson.coeff)
               }

               #Spearman
               if("spearman" %in% tolower(solver.list)){
                   out.list$out.spearman$gene <- rownames(out.list$out.spearman)
                   rownames(out.list$out.spearman) <- NULL
                   names(out.list$out.spearman) <- c("spearman.coeff", "gene")
                   spearman.med <- stats::median(out.list$out.spearman$spearman.coeff)
                   spearman.scale <- stats::mad(out.list$out.spearman$spearman.coeff)
               }

               #LassoPV
               if("lassopv" %in% tolower(solver.list)){
                   out.list$out.lassopv$gene <- rownames(out.list$out.lassopv)
                   rownames(out.list$out.lassopv) <- NULL
                   out.list$out.lassopv <- out.list$out.lassopv[, c("p.values","gene")]
                   names(out.list$out.lassopv) <- c("lasso.p.value", "gene")
                   p.log10 <- -log10(out.list$out.lassopv$lasso.p.value)
                   lassopv.med <- stats::median(p.log10)
                   lassopv.scale <- sqrt(mean(p.log10*p.log10))
               }

               #SqrtLasso
               if("sqrtlasso" %in% tolower(solver.list)){
                   out.list$out.sqrtlasso$gene <- rownames(out.list$out.sqrtlasso)
                   rownames(out.list$out.sqrtlasso) <- NULL
                   out.list$out.sqrtlasso <- out.list$out.sqrtlasso[, c("beta", "gene")]
                   names(out.list$out.sqrtlasso) <- c("beta.sqrtlasso", "gene")
                   sqrtlasso.med <- stats::median(out.list$out.sqrtlasso$beta.sqrtlasso)
                   sqrtlasso.scale <- stats::mad(out.list$out.sqrtlasso$beta.sqrtlasso)
               }

               #Ridge
               if("ridge" %in% tolower(solver.list)){
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
                                              center = bayesspike.med,
                                              scale = bayesspike.scale)
                   }

               if("rf.score" %in% names(tbl.scale)){
                   tbl.scale$rf.score <- scale(tbl.scale$rf.score,
                                               center = randomforest.med,
                                               scale = randomforest.scale)
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
