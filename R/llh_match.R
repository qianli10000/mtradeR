
#' @importFrom stats  model.matrix p.adjust pchisq quantile rnorm

llh_match<-function(par_all,input,shrinkage){
  
  if(is.null(input$others_abun_1) & is.null(input$others_pres_1)){
    n_others_abun=0
    n_others_pres=0
    par=par_all
    beta01=beta02<-NULL
  }else if(is.null(input$others_abun_1) & (!is.null(input$others_pres_1)) ){
    n_others_abun=0
    n_others_pres=ncol(input$others_pres_1)
    par=par_all[-(1:(n_others_pres))]
    beta01<-NULL
    beta02<-par_all[1:n_others_pres]
  }else if((!is.null(input$others_abun_1)) & is.null(input$others_pres_1)){
    n_others_abun=ncol(input$others_abun_1)
    n_others_pres=0
    par=par_all[-(1:(n_others_abun))]
    beta01<-par_all[1:n_others_abun]
    beta02<-NULL
  } else{
    n_others_abun=ncol(input$others_abun_1)
    n_others_pres=ncol(input$others_pres_1)
    par=par_all[-(1:(n_others_abun+n_others_pres))]
    beta01<-par_all[1:n_others_abun]
    beta02<-par_all[n_others_abun+(1:n_others_pres)]
  }
  
  n_par=length(par)
 
  long_dim<-ncol(input$long_design_1)
  logistic_dim<-ncol(input$logistic_design_1)
  
  beta1<-par[1:long_dim]
  lamda1<-par[(long_dim+1)]
  gamma1<-par[(long_dim+2)]
  beta2<-par[(long_dim+3):(2*long_dim+2)]
  lamda2<-par[(2*long_dim+3)]
  gamma2<-par[(2*long_dim+4)]
  alph<-par[(2*long_dim+5):(2*long_dim+4+logistic_dim)]
  phi<-par[n_par-2]
  sig_sub<-par[n_par-1]
  sig_set<-par[n_par]
  
  llkh1_1=log(lh1_match(input$taxa_1,input$long_design_1,input$others_abun_1,input$others_pres_1,input$gh.nodes1*sig_sub*input$X_rand1*sqrt(2),
                        input$gh.nodes_set1*sig_set*input$X_rand1*sqrt(2),beta1,beta01,beta2,beta02,lamda1,lamda2,gamma1,gamma2,phi))
  llkh2_1=log(lh2_match(input$outcome_1,input$logistic_design_1,input$gh.nodes1_sub*sig_sub*sqrt(2),input$gh.nodes_set1_sub*sig_set*sqrt(2),alph))
  llkh1_2=log(lh1_match(input$taxa_2,input$long_design_2,input$others_abun_2,input$others_pres_2,input$gh.nodes2*sig_sub*input$X_rand2*sqrt(2),
                        input$gh.nodes_set2*sig_set*input$X_rand2*sqrt(2),beta1,beta01,beta2,beta02,lamda1,lamda2,gamma1,gamma2,phi))
  llkh2_2=log(lh2_match(input$outcome_2,input$logistic_design_2,input$gh.nodes2_sub*sig_sub*sqrt(2),input$gh.nodes_set2_sub*sig_set*sqrt(2),alph))
  
  joint.logl_1=tcrossprod(input$prod.mat1,t(llkh1_1))+llkh2_1
  joint.logl_2=tcrossprod(input$prod.mat2,t(llkh1_2))+llkh2_2
    
  
  l=-log(rowSums(input$gh.weights*exp(joint.logl_1+joint.logl_2)))  
  
  
  l=ifelse(is.na(l)|is.infinite(l),0,l)
  penalty=shrinkage*sum(par^2)
  return((sum(l)+penalty))
}


prep<-function(taxa,others_abun,others_pres,long_design,logistic_design,outcome,long_idset,logistic_idset,rand.var){
  quad.n<-5
  nw <- mvQuad::createNIGrid(dim=3, type="GHe", level=quad.n)
  
  long_dim<-ncol(long_design)
  logistic_dim<-ncol(logistic_design)
  
  long_id<-long_idset[,1]
  long_order<-long_idset[,3]
  logistic_order<-logistic_idset[,3] 
  
   if(is.null(others_abun) & is.null(others_pres)){
    others_abun_1=others_abun_2=others_abun
    others_pres_1=others_pres_2=others_pres
  }else if( is.null(others_abun) & (!is.null(others_pres)) ){ 
    others_abun_1=others_abun_2=others_abun
    others_pres_1=as.matrix(others_pres[long_order==1,])
    others_pres_2=as.matrix(others_pres[long_order==2,])
  }else if( (!is.null(others_abun)) & is.null(others_pres) ){
    others_pres_1=others_pres_2=others_pres
    others_abun_1=as.matrix(others_abun[long_order==1,])
    others_abun_2=as.matrix(others_abun[long_order==2,])
  }else{
    others_pres_1=as.matrix(others_pres[long_order==1,])
    others_pres_2=as.matrix(others_pres[long_order==2,])
    others_abun_1=as.matrix(others_abun[long_order==1,])
    others_abun_2=as.matrix(others_abun[long_order==2,])
  }
  
  
  long.size1=sum(long_order==1)
  long.size2=sum(long_order==2)
  logistic.size1=sum(logistic_order==1)
  logistic.size2=sum(logistic_order==2)
  gh.nodes1=matrix(rep(nw$nodes[,1],long.size1),nrow = long.size1,byrow = T)
  gh.nodes2=matrix(rep(nw$nodes[,2],long.size2),nrow = long.size2,byrow = T)
  gh.nodes_set1=matrix(rep(nw$nodes[,3],long.size1),nrow = long.size1,byrow = T)
  gh.nodes_set2=matrix(rep(nw$nodes[,3],long.size2),nrow = long.size2,byrow = T)
  gh.weights=matrix(rep(nw$weights,logistic.size1),nrow = logistic.size1,byrow = T)
  
  
  prod.mat1=t(sapply(unique(long_id[long_order==1]),FUN = function(x){
    ifelse(long_id[long_order==1]==x,1,0)
  }))
  #prod.mat1=bdiag(prod.mat1)
  
  prod.mat2=t(sapply(unique(long_id[long_order==2]),FUN = function(x){
    ifelse(long_id[long_order==2]==x,1,0)
  }))
  #prod.mat2=bdiag(prod.mat2)
  
  long_design_1=long_design[long_order==1,]
  long_design_2=long_design[long_order==2,]
  X_rand1<-long_design_1[,rand.var]
  X_rand2<-long_design_2[,rand.var]
  
  outcome_1=outcome[logistic_order==1]
  outcome_2=outcome[logistic_order==2]
  logistic_design_1=logistic_design[logistic_order==1,]
  logistic_design_2=logistic_design[logistic_order==2,]
  
  
  taxa_1=as.numeric(taxa[long_order==1])
  taxa_2=taxa[long_order==2]
  
  
  gh.nodes1_sub=gh.nodes1[1:logistic.size1,]
  gh.nodes2_sub=gh.nodes2[1:logistic.size2,]
  gh.nodes_set1_sub=gh.nodes_set1[1:logistic.size1,]
  gh.nodes_set2_sub=gh.nodes_set2[1:logistic.size2,]
  
  output=list(
    taxa_1=taxa_1,logistic_design_1=logistic_design_1,long_design_1=long_design_1,others_abun_1=others_abun_1,others_pres_1=others_pres_1,
    gh.nodes1=gh.nodes1,gh.nodes_set1=gh.nodes_set1,X_rand1=X_rand1,outcome_1=outcome_1,
    taxa_2=taxa_2,logistic_design_2=logistic_design_2,long_design_2=long_design_2,others_abun_2=others_abun_2,others_pres_2=others_pres_2,
    gh.nodes2=gh.nodes2,gh.nodes_set2=gh.nodes_set2,X_rand2=X_rand2,outcome_2=outcome_2,
    prod.mat1= prod.mat1, prod.mat2= prod.mat2, gh.weights=gh.weights, gh.nodes1_sub=gh.nodes1_sub,gh.nodes2_sub=gh.nodes2_sub,
    gh.nodes_set1_sub=gh.nodes_set1_sub,gh.nodes_set2_sub=gh.nodes_set2_sub)
  
  return(output)
  
}