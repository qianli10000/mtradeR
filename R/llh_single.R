#' @export

llh_single<-function(par,taxa,longi_design_all,logistic_design_all,outcome,longi_idset,logistic_idset,rand.var,shrinkage){

  quad.n<-30
  longi_dim=ncol(longi_design_all)
  logistic_dim=ncol(logistic_design_all)

  beta1<-par[1:longi_dim]
  lamda1=par[(longi_dim+1)]

  beta2<-par[(longi_dim+2):(2*longi_dim+1)]
  lamda2=par[(2*longi_dim+2)]

  alpha<-par[(2*longi_dim+3):(2*longi_dim+2+logistic_dim)]
  phi<-par[length(par)-1]
  sig_sub<-par[length(par)]

  nw <- mvQuad::createNIGrid(dim=1, type="GHe", level=quad.n)

  t=shrinkage

  longi_id=longi_idset[,1]
  longi_order=longi_idset[,2]
  logistic_order=logistic_idset[,2]   # longi_order, logistic_order should be created in wrap function, using Input #2 pair indicator


  longi.size1=sum(longi_order==1)
  longi.size2=sum(longi_order==2)
  logistic.size=sum(logistic_order==1)
  gh.nodes1=matrix(rep(nw$nodes,longi.size1),nrow = longi.size1,byrow = T)
  gh.nodes2=matrix(rep(nw$nodes,longi.size2),nrow = longi.size2,byrow = T)

  gh.weights=matrix(rep(nw$weights,logistic.size),nrow = logistic.size,byrow = T)

  gh.nodes1_sub=gh.nodes1[1:logistic.size,]
  gh.nodes2_sub=gh.nodes2[1:logistic.size,]

  taxa_1=taxa[longi_order==1]
  taxa_2=taxa[longi_order==2]
  longi_design_all_1=longi_design_all[longi_order==1,]
  longi_design_all_2=longi_design_all[longi_order==2,]
  X_rand1=longi_design_all_1[,rand.var]
  X_rand2=longi_design_all_2[,rand.var]

  outcome_1=outcome[logistic_order==1]
  outcome_2=outcome[logistic_order==2]
  logistic_design_all_1=logistic_design_all[logistic_order==1,]
  logistic_design_all_2=logistic_design_all[logistic_order==2,]




  prod.mat1<-t(sapply(unique(longi_id[longi_order==1]),FUN = function(x){
    ifelse(longi_id[longi_order==1]==x,1,0)
  }))
  prod.mat2<-t(sapply(unique(longi_id[longi_order==2]),FUN = function(x){
    ifelse(longi_id[longi_order==2]==x,1,0)
  }))

  joint.logl_1=prod.mat1%*%log(lh1_single(taxa_1,phi,longi_design_all_1,beta1,beta2,lamda1,lamda2,gh.nodes1*sig_sub*X_rand1*sqrt(2)))+
    log(lh2_single(outcome_1,alpha,logistic_design_all_1,gh.nodes1_sub*sig_sub*sqrt(2)))
  joint.logl_2=prod.mat2%*%log(lh1_single(taxa_2,phi,longi_design_all_2,beta1,beta2,lamda1,lamda2,gh.nodes2*sig_sub*X_rand2*sqrt(2)))+
    log(lh2_single(outcome_2,alpha,logistic_design_all_2,gh.nodes2_sub*sig_sub*sqrt(2)))


  l=-log(rowSums(gh.weights*exp(joint.logl_1+joint.logl_2)))


  l=ifelse(is.na(l)|is.infinite(l),0,l)
  penalty=t*sum(par^2)
  return((sum(l)+penalty))

}
