
#' @importFrom stats  model.matrix p.adjust pchisq quantile rnorm

lh1_match<-function(y,x,others_abun,others_pres,rand_sub,rand_set,beta1,beta01,beta2,beta02,lamda1,lamda2,gamma1,gamma2,phi){
  quad.n<-3
  Q<-quad.n^3
  
  y=as.numeric(y)
  beta1m=matrix(rep(beta1,Q),ncol = Q)
  if(is.null(others_abun)){
    mu=1/(1+exp(-(tcrossprod(x,t(beta1m))+lamda1*rand_sub+gamma1*rand_set)))
  }else {
    beta01m=matrix(rep(beta01,Q),ncol = Q)
    mu=1/(1+exp(-(tcrossprod(others_abun,t(beta01m))+tcrossprod(x,t(beta1m))+lamda1*rand_sub+gamma1*rand_set)))
  }
  mu[mu<1e-10]=1e-10
  mu[mu>(1-1e-10)]=1-1e-10
  positive.ind=ifelse(y>0,1,0)
  
  y_m=matrix(rep(y,Q),ncol=Q)
  nonzero.density=ifelse(y_m>0,(y_m^(phi*mu-1))*((1-y_m)^(phi*(1-mu)-1))*gamma(phi)/(gamma(phi*mu)*gamma(phi*(1-mu))),0)
  
  
  beta2m=matrix(rep(beta2,Q),ncol = Q)
  if(is.null(others_pres)){
    presence=1/(1+exp(-(tcrossprod(x, t(beta2m))+lamda2*rand_sub+gamma2*rand_set)))
  } else {
    beta02m=matrix(rep(beta02,Q),ncol = Q)
    presence=1/(1+exp(-(tcrossprod(others_pres,t(beta02m))+tcrossprod(x,t(beta2m))+lamda2*rand_sub+gamma2*rand_set)))
  }
  
  presence[presence<1e-4]=1e-4
  presence[presence>(1-1e-4)]=1-1e-4
  
  
  l1=(1-positive.ind)*(1-presence)+positive.ind*presence*nonzero.density
  return(l1)
}
