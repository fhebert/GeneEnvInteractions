GEInterScoreTest = function(X,E,Y,perm.method="parametric",N=1000,U=NULL){
  GG = PCsel(X,0.001*ncol(X))
  W = cbind(1,E,GG,U)
  S = matrix(rep(E,ncol(X)),ncol=ncol(X))*X
  hatY0 = LogitProbs(W,Y)
  Z = as.vector(t(Y-hatY0)%*%S)
  s2Y = sum((Y-hatY0)^2)/length(Y)
  WtS = t(W)%*%(S)
  A = (crossprod(S)-t(WtS)%*%solve(crossprod(W))%*%WtS)
  Gamma = s2Y*A
  d = sqrt(diag(Gamma))
  Z = Z/d
  D = matrix(rep(d,ncol(S)),ncol=ncol(S))
  Sigma = Gamma/(D*t(D))
  if(perm.method=="parametric"){
    Y0 = ParametricResamplingGEInter(hatY0,N)
    hatY0 = PermLogitProbs(W,Y0)
  }else{
    Y0 = sapply(1:N,function(i){sample(Y)})
    hatY0 = PermLogitProbs(W,Y0)
  }
  d0 = Y0-hatY0
  Z0 = t(d0)%*%S
  sY0 = sqrt((rep(1,length(Y))%*%(d0^2))[1,]/length(Y))
  dA = sqrt(diag(A))
  Z0 = Z0/(matrix(sY0,ncol=1)%*%matrix(dA,nrow=1))
  Z0 = as.matrix(Z0)
  diag(Sigma) = 1
  return(list(Z=Z,Z0=Z0,Sigma=Sigma,Y0=Y0,S=S,W=W))
}
