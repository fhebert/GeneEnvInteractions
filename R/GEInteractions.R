PopulationPhenotypeGEInter = function(X,E,beta0,betaX,I.betaX,betaE,gamma,I.gamma,U=NULL,betaU=NULL){
  S = matrix(rep(E,ncol(X)),ncol=ncol(X))*X
  Y = beta0+E*betaE+X[,I.betaX,drop=FALSE]%*%betaX+S[,I.gamma,drop=FALSE]%*%gamma
  if(!is.null(U)){
    Y = Y+U%*%betaU
  }
  Y = 1/(1+exp(-Y))
  Y = rbinom(length(Y),1,Y)
  return(Y)
}

ParametricResamplingGEInter = function(probs,N=1000){
  Y0 = sapply(1:N,function(i){rbinom(length(probs),1,probs)})
  return(Y0)
}

PCsel = function(X,thresh){
  X = scale(X)
  S = cor(X)
  ev = eigen(S)
  ind = which(ev$values>=thresh)
  XX = X%*%ev$vectors[,ind]
  A = matrix(1,ncol=1,nrow=nrow(X))%*%matrix(sqrt(ev$values[ind]),nrow=1)
  XX = XX/A
  return(XX)
}