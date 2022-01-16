#' @export

lh2_match<-function(o,alpha,z,rand_sub,rand_set){
  quad.n<-4
  alpha=matrix(rep(alpha,quad.n^3),ncol = quad.n^3)
  prob=1/(1+exp(-(z%*%alpha+rand_sub+rand_set)))
  l2=o*prob+(1-o)*(1-prob)
  return(l2)
}
