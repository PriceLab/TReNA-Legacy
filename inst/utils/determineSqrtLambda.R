library(flare)

#----------------------------------------------------------------------------------------------------
# Determine lambda for square root lasso via permutation testing and a brute force method

bruteForceLambda <- function(nlambda = 1000){

    # 0) Load and segregate the data
    load(system.file(package="TReNA", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
    target.gene <- "MEF2C"
    mtx.asinh <- asinh(mtx.sub)
    tfs <- setdiff(rownames(mtx.asinh), "MEF2C")
    features <- t(mtx.asinh[tfs,,drop=F])
    target <- as.numeric(mtx.asinh[target.gene,])

    # 1) Permute the response vector
    set.seed(101010)
    target <- sample(target)

    # 2) Create lambda values and find the fit for each one
    fit <- slim(features, target, method = "lq", verbose = FALSE, nlambda = nlambda)

    # 3) Find the lambda threshold below which the algorithm returns something other than nonsense
    threshold <- 10^(-15)
    below <- NULL

    for(i in 1:ncol(fit$beta)){
        if (max(fit$beta[,i] > threshold)){
            below <- i
            break
        }}

    return(fit$lambda[[below-1]])
}

#----------------------------------------------------------------------------------------------------
# Determine lambda for square root lasso via permutation testing and a binary method

findBinaryLambda <- function(){

    # 0) Load and segregate the data
    load(system.file(package="TReNA", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
    target.gene <- "MEF2C"
    mtx.asinh <- asinh(mtx.sub)
    tfs <- setdiff(rownames(mtx.asinh), "MEF2C")
    features <- t(mtx.asinh[tfs,,drop=F])
    target <- as.numeric(mtx.asinh[target.gene,])

    # 1) Permute the response vector
    set.seed(101010)
    target <- sample(target)

    # 2) Create a threshold, a measure for the change in lambda and a starting value

    threshold <- 10^(-15)
    lambda.change <- 10^(-7)
    lambda <- 1

    # 3) Write a while loop that does a binary search

    step.size <- lambda/2 # Start at 0.5
    while(step.size > lambda.change){

        # Get the fit
        fit <- slim(features, target, method = "lq", verbose = FALSE, lambda = lambda)

        # Evaluate the betas and change lambda
        # Case 1: nonsense, need to lower lambda
        if(max(fit$beta) < threshold){
            lambda <- lambda - step.size
        }
        # Case 2: sense, need to raise lambda
        else{
            lambda <- lambda + step.size
        }

        # Halve the step size
        step.size <- step.size/2
    }

    # 4) Return the final lambda
    return(lambda)
}
