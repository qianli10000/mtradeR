#' @export

lh2_match<-function(o,z,rand_sub,rand_set){
  
  alpham=matrix(rep(alph,Q),ncol = Q)
  prob=1/(1+exp(-(crossprod(t(z),alpham)+rand_sub+rand_set)))
  l2=o*prob+(1-o)*(1-prob)
  return(l2)
}
