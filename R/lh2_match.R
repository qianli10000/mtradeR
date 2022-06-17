
#' @importFrom stats  model.matrix p.adjust pchisq quantile rnorm

lh2_match<-function(o,z,rand_sub,rand_set,alph){
  quad.n<-3
  Q<-quad.n^3
  alpham=matrix(rep(alph,Q),ncol = Q)
  prob=1/(1+exp(-(crossprod(t(z),alpham)+rand_sub+rand_set)))
  l2=o*prob+(1-o)*(1-prob)
  return(l2)
}
