mcTRN<-function(expr, physnet, N=20, ratio=0.8,lambda = 0.1, cores = 4, alpha = 1){
  # Monte Carlo version of TRN: building TRN using resampled sub-datasets;
  # N: number of resampling times
  # ratio: percent of samples used to build TRN, 0.8, defaults
  # Hongdong@ISB, Sep. 23, 2015
  

  # basic parameters
  if (!is.matrix(expr)) {expr=as.matrix(expr)}  
  ngene=nrow(expr)
  nsample=ncol(expr)
  nsubset=floor(nsample*ratio)
  tf=intersect(rownames(expr),colnames(physnet))
  ntf=length(tf)
  allgenes=rownames(expr)
  
  
  
  # create variables to store information needed
  freq=matrix(0,nrow=ngene,ncol=ntf)
  rownames(freq)=allgenes;
  colnames(freq)=tf
  
  
  beta_mean=matrix(0,nrow=ngene,ncol=ntf)
  rownames(beta_mean)=allgenes
  colnames(beta_mean)=tf
  
  
  beta_sd=matrix(0,nrow=ngene,ncol=ntf)
  rownames(beta_sd)=allgenes
  colnames(beta_sd)=tf
  
  
  tmplate=matrix(0,nrow=ngene,ncol=ntf)
  rownames(tmplate)=allgenes
  colnames(tmplate)=tf
  
  # loop for building TRN using resampled sub-datasets
  for (i in 1:N){
    # subdata
    subindex=sample.int(nsample,nsubset)
    subsample=expr[,subindex]

    
    # build TRN
    trn=trimTRN(fitTRN(subsample, physnet,lambda,cores,alpha))
    
    
    # tmplate
    beta=tmplate # recover to original full size
    beta[rownames(trn),colnames(trn)]=trn
    

    # record frequency
    knonzero=which(abs(beta)>0.000000001)
    freq[knonzero]=freq[knonzero]+1
   
    
    # update mean and sd  
    for (j in 1:length(knonzero)){
      address=knonzero[j]
      
      m0=beta_mean[address]
      s0=beta_sd[address]
      xnew=beta[address]
      nnew=freq[address]
      
      
      m1=(m0*(nnew-1)+xnew)/nnew
      s1=s0+(xnew-m0)*(xnew-m1)
      
      beta_mean[address]=m1
      beta_sd[address]=s1
      
    }
  
    rm(beta)
  } # end of each subset
  
 
  # calculate SD
  knonzero=which(freq>1)
  for (j in 1:length(knonzero)){
    address=knonzero[j]
    beta_sd[address]=sqrt(beta_sd[address]/(freq[address]-1))
  }
  
  
  ## clear unused
  rm(tmplate)
  
  # output beta and frequency
  output=list(freq=freq,beta_mean=beta_mean,beta_sd=beta_sd,N=N,ratio=ratio)
  
  return(output)
  
}


