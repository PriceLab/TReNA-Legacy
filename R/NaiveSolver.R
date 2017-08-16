# NaiveSolver Prototype
# Utilizing lasso to eliminate variables and allowing lm to be more accurate
# Be warned! This is known to function much better on "larger" data
#----------------------------------------------------------------------------------------------------
# solver class defintion
.NaiveSolver <- setClass("NaiveSolver",
                         contains = "Solver",
                         slots = c(alpha = "numeric",
                                   lambda = "numeric",
                                   rfSelector = "logical"))
#----------------------------------------------------------------------------------------------------
# solver constructor function
NaiveSolver <- function(mtx.assay = matrix(), targetGene, candidateRegulators,
                        alpha = 1, lambda = numeric(), rfSelector = FALSE, quiet = TRUE)
{
  # Remove the targetGene from candidateRegulators
  if(any(grepl(targetGene, candidateRegulators)))      
    candidateRegulators <- candidateRegulators[-grep(targetGene, candidateRegulators)]    
  
  # Check to make sure the matrix contains some of the candidates
  candidateRegulators <- intersect(candidateRegulators, rownames(mtx.assay))    
  stopifnot(length(candidateRegulators) > 0)
  
  obj <- .NaiveSolver(mtx.assay = mtx.assay,targetGene = targetGene,
                      candidateRegulators = candidateRegulators, alpha = alpha,
                      lambda = lambda, rfSelector = rfSelector, quiet = quiet)
  
  obj
}
#----------------------------------------------------------------------------------------------------
setMethod('show', 'NaiveSolver',
          
          function(obj) {
            regulator.count <- length(getRegulators(obj))
            if(regulator.count > 10){
              regulatorString <- paste(getRegulators(obj)[1:10], collapse=",")
              regulatorString <- sprintf("%s...", regulatorString);
            }
            else
              regulatorString <- paste(getRegulators(obj), collapse=",")
            
            msg = sprintf("NaiveSolver with mtx.assay (%d, %d), targetGene %s, %d candidate regulators %s, alpha = %f",
                          nrow(getAssayData(obj)), ncol(getAssayData(obj)),
                          getTarget(obj), regulator.count, regulatorString, obj@alpha)
            cat (msg, '\n', sep='')
          })
#----------------------------------------------------------------------------------------------------
# run function
setMethod("run", "NaiveSolver",
          function (obj){
            
            # begin with a TReNA object and extract the slots/data
            mtx <- getAssayData(obj)
            target.gene <- getTarget(obj)
            tfs <- getRegulators(obj)
            alpha <- obj@alpha
            lambda <- obj@lambda
            rfSelector <- obj@rfSelector
            
            # Check if target.gene is in bottom 10% mean expressions
            if(rowMeans(mtx)[target.gene] < stats::quantile(rowMeans(mtx), probs = 0.1)){
              warning("Target gene mean expression is in the bottom 10% of all genes in the assay matrix")
            }
            
            # Check if data is optimal
            if(dim(mtx)[1] > dim(mtx)[2]){
              warning("There are fewer variables than entries (this solver may not be optimal for this data)")
            }
            
            # We also have to make a few checks to determine the data is formatted/submitted properly
            stopifnot(target.gene %in% rownames(mtx))
            stopifnot(all(tfs %in% rownames(mtx)))
            if (length(tfs) == 0)
              return(NULL)
            
            deleters <- grep(target.gene, tfs)
            if (length(deleters) > 0){
              tfs <- tfs[-deleters]
            }
            if (length(tfs) == 0)
              return(NULL)
            
            features <- t(mtx[tfs,,drop=FALSE])
            target <- as.numeric(mtx[target.gene,])
            tf.weights <- rep(1, length(tfs))
            
            if(rfSelector){
            # RandomForest Filtering Method
              # Solve with RF
              solver <- RandomForestSolver(mtx,target.gene,tfs)
              tbl <- run(solver)$edges
              # Pull out threshold %
              # Say threshold for now is 10%
              threshold <- nrow(tbl) * 0.1
              selected.names <- rownames(tbl[1:threshold,]) %>%
                append(target.gene)
              selected.matrix <- as.data.frame(t(mtx[which(rownames(mtx) %in% selected.names),]))
              # run lm on that
              lin.mod <- lm(paste0(target.gene,"~."),data = selected.matrix)
              # return formatted betas/p.vals
              coef.summary <- summary(lin.mod)$coefficients
              tbl <- data.frame(row.names = rownames(coef.summary),
                                beta = coef.summary[,1],
                                p.value = coef.summary[,4])
              tbl <- tbl[-1,,drop=FALSE]
              tbl <- tbl[order(tbl$p.value),]
              return(tbl)
            }else{
            # Lasso Filtering Method
            
            # glmnet prep
            if( length(lambda) == 0 ) {
              
              # Run Permutation testing to find lambda
              if( alpha != 0 )
                alpha.perm = alpha
              else
                (alpha.perm = 0.1)
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
            
            else if(is.numeric(lambda)){
              fit = glmnet(features, target, penalty.factor=tf.weights, alpha=alpha, lambda=lambda)
            }
            
            # Pull out the non-zero coefficients and those matching variables
            fit.coefs <- as.matrix(coef(fit))
            nz.indices <- which(!(fit.coefs == 0))
            non.zeroes <- rownames(as.matrix(fit.coefs[nz.indices,,drop=FALSE]))[-1]
            non.zeroes <- append(non.zeroes, target.gene)
            selected.matrix <- mtx[which(rownames(mtx) %in% non.zeroes),]
            
            # lm the new matrix
            selected.matrix <- as.data.frame(t(selected.matrix))
            lin.mod <- lm(paste0(target.gene,"~."),data = selected.matrix)
            
            # Create, format, and return the data frame from the lm
            coef.summary <- summary(lin.mod)$coefficients
            tbl <- data.frame(row.names = rownames(coef.summary),
                              beta = coef.summary[,1],
                              p.value = coef.summary[,4])
            tbl <- tbl[-1,,drop=FALSE]
            tbl <- tbl[order(tbl$p.value),]
            return(tbl)
          }
            })

#----------------------------------------------------------------------------------------------------