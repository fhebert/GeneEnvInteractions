SampleGEInter = function(X,E,beta0,betaX,I.betaX,betaE,gamma,I.gamma,n0,n1,U=NULL,betaU=NULL){
  Y = PopulationPhenotypeGEInter(X,E,beta0,betaX,I.betaX,betaE,gamma,I.gamma,U=NULL,betaU=NULL)
  i = which(Y==0)
  ind = sample(c(sample(i,n0,FALSE),sample((1:length(Y))[-i],n1,FALSE)))
  GG = X[ind,]
  EE = matrix(E[ind],ncol=1)
  Y = matrix(Y[ind],ncol=1)
  U = U[ind,]
  return(list(SNP=GG,Env=EE,Phenotype=Y,Covariates=U))
}
