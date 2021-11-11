#' @export

lh2_single<-function(o,alpha,z,rand_sub){
  quad.n<-30
  alpha=matrix(rep(alpha,quad.n),ncol = quad.n)
  prob=1/(1+exp(-(z%*%alpha+rand_sub)))
  l2=o*prob+(1-o)*(1-prob)
  return(l2)
}
