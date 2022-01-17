#' @export


lh1_match<-function(y,phi,x,beta1,beta2,lamda1,lamda2,gamma1,gamma2,rand_sub,rand_set){
  quad.n=3
  Q=quad.n^3
  beta1=matrix(rep(beta1,Q),ncol = Q)
  mu=1/(1+exp(-(x%*%beta1+lamda1*rand_sub+gamma1*rand_set)))
  # mu[mu<1e-12]=1e-12
  # mu[mu>(1-1e-12)]=1-1e-12
  positive.ind=ifelse(y>0,1,0)


  #nonzero.density=nonzero_density(y,mu,phi)

  y_m=matrix(rep(y,Q),ncol=Q)
  nonzero.density=ifelse(y_m>0,(y_m^(phi*mu-1))*((1-y_m)^(phi*(1-mu)-1))*gamma(phi)/(gamma(phi*mu)*gamma(phi*(1-mu))),0)


  beta2=matrix(rep(beta2,Q),ncol = Q)
  presence=1/(1+exp(-(x%*%beta2+lamda2*rand_sub+gamma2*rand_set)))
  # presence[presence<1e-4]=1e-4
  # presence[presence>(1-1e-4)]=1-1e-4


  l1=(1-positive.ind)*(1-presence)+positive.ind*presence*nonzero.density
  return(l1)
}

