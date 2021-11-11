#' @export

lh1_single<-function(y,phi,x,beta1,beta2,lamda1,lamda2,rand_sub){
  quad.n<-30
  beta1=matrix(rep(beta1,quad.n),ncol = quad.n)
  mu=1/(1+exp(-(x%*%beta1+lamda1*rand_sub)))
  mu[mu<1e-10]=1e-10
  mu[mu>(1-1e-10)]=1-1e-10
  positive.ind=ifelse(y>0,1,0)
  
  
  nonzero.density=apply(mu,2,function(s){
    d=ifelse(y>0,(y^(phi*s-1))*((1-y)^(phi*(1-s)-1))*gamma(phi)/(gamma(phi*s)*gamma(phi*(1-s))),0)
    return(d)}
  )
  
  beta2=matrix(rep(beta2,quad.n),ncol = quad.n)
  presence=1/(1+exp(-(x%*%beta2+lamda2*rand_sub)))
  presence[presence<1e-4]=1e-4
  presence[presence>(1-1e-4)]=1-1e-4
  
  
  l1=(1-positive.ind)*(1-presence)+positive.ind*presence*nonzero.density
  return(l1)
}

