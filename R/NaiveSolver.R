# NaiveSolver Prototype
# Utilizing lasso to eliminate variables and allowing lm to be more accurate
# Be warned! This is known to function much better on "larger" data
#----------------------------------------------------------------------------------------------------
# solver class defintion
.NaiveSolver <- setClass("NaiveSolver", contains = "Solver")
#----------------------------------------------------------------------------------------------------
# solver constructor function
NaiveSolver <- function(mtx.assay = matrix(), quiet = TRUE)
{
  obj <- .NaiveSolver(Solver(mtx.assay = mtx.assay, quiet = quiet))
  
  obj
}
#----------------------------------------------------------------------------------------------------
# run function
setMethod("run", "NaiveSolver",
          function (obj, target.gene, tfs, tf.weights=rep(1,length(tfs)), extraArgs=list()){
            
            # begin with a TReNA object and extract the slots/data
            mtx <- getAssayData(obj)
            
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
            if (length(tfs) == 0) return(NULL)
            
            deleters <- grep(target.gene, tfs)
            if (length(deleters) > 0){
              tfs <- tfs[-deleters]
            }
            if (length(tfs) == 0) return(NULL)
            
            # Default setting
            alpha <- 1
            lambda <- NULL
            
            if("alpha" %in% names(extraArgs))
              alpha <- extraArgs[["alpha"]]
            if("lambda" %in% names(extraArgs))
              lambda <- extraArgs[["lambda"]]
            
            features <- t(mtx[tfs,,drop=FALSE])
            target <- as.numeric(mtx[target.gene,])
            tf.weights <- rep(1, length(tfs))
            
            # glmnet prep
            if( is.null(lambda) ) {
              
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
            selected.matrix <- mtx.asinh[which(rownames(mtx.asinh) %in% non.zeroes),]
            
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
            invisible(tbl)
          })

#----------------------------------------------------------------------------------------------------