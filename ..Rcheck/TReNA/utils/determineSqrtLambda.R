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

#----------------------------------------------------------------------------------------------------
# Create dummy data and verify the LASSO and SqrtLasso both work

testLassos <- function(){

    # Create the dummy data
    n = 100
    p = 200
    sigma = 0.1

    # ten nonzero elements in beta
    beta_nz_data <- 1:10
    beta_nz_data <- sapply(beta_nz_data,function(x){x*(-1)^x})    
    beta_data <- c(beta_nz_data,numeric(p-length(beta_nz_data)))

    # X mean centered and normalized so sum of sqs in each col = N
    X_data <- matrix(rnorm(n*p),nrow=n,ncol=p)
    X_data <- scale(X_data)

    # y mean centered and normalized so sum of sqs = N
    y_data <- X_data %*% beta_data + sigma*rnorm(n);
    y_data <- scale(y_data)

    # Transpose matrix to work with solvers and add rownames
    mtx.assay <- t(cbind(X_data,y_data))
    rownames(mtx.assay) <- paste("Gene", 1:(p+1), sep="_")
    
    # Try fitting it using the LASSO
    trena <- TReNA(mtx.assay=mtx.assay,solver="lasso",quiet=FALSE)
    tfs <- setdiff(rownames(mtx.assay), "Gene_201")
    tbl.lasso <- solve(trena, "Gene_201", tfs)

    
    # Try fitting it using the Sqrt LASSO
    trena <- TReNA(mtx.assay=mtx.assay,solver="sqrtlasso",quiet=FALSE)
    tbl.sqrtlasso <- solve(trena, "Gene_201", tfs)
    
}
