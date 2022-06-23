

JMR_core<-function(taxa,others_abun,others_pres,long_design,logistic_design,outcome,long_idset,logistic_idset,rand.var,shrinkage,trace=T,n.cores=NULL){
 
  long_dim<-ncol(long_design)
  logistic_dim<-ncol(logistic_design)

    if(is.null(others_abun)){
      P_others_abun=0
    }else{ 
      P_others_abun=ncol(others_abun)
    }
    
    if(is.null(others_pres)){
      P_others_pres=0
    }else{ 
      P_others_pres=ncol(others_pres)
    }
    
    P<-P_others_abun+P_others_pres+2*long_dim+logistic_dim+8
    
<<<<<<< HEAD
    # if(is.null(n.cores)){
    #   n.cores=min(parallel::detectCores()-1,16)
    # }
    # cl <- parallel::makeCluster(n.cores)     # set the number of processor cores

  #parallel::setDefaultCluster(cl=cl) # set 'cl' as default cluster
=======
    if(is.null(n.cores)){
    n.cores=min(parallel::detectCores()-1,16)
    }
    cl <- parallel::makeCluster(n.cores)     # set the number of processor cores
    
  parallel::setDefaultCluster(cl=cl) # set 'cl' as default cluster
>>>>>>> 54a163a8764e7463463f5d3cb7f735d320d7409e
  par.ini=c(composition_others=rep(0,P_others_abun),presence_others=rep(0,P_others_pres),composition_beta=rep(0,long_dim),composition_lambda=0,composition_gamma=0,presence_beta=rep(0,long_dim),presence_lambda=0,presence_gamma=0,outcome_alpha=rep(0,logistic_dim),dispersion=1,variance_a=1,variance_b=1)
  prep.data=prep(taxa,others_abun,others_pres,long_design,logistic_design,outcome,long_idset,logistic_idset,rand.var)
  
  m=optimParallel::optimParallel(par.ini,llh_match,input=prep.data,shrinkage=shrinkage,#taxa=taxa,others_abun=others_abun,others_pres=others_pres,long_design=long_design,logistic_design=logistic_design,outcome=outcome,long_idset=long_idset,logistic_idset=logistic_idset,rand.var=rand.var,
                                 lower = c(rep(-Inf,(P-3)),1e-05,1e-05,1e-05),upper=rep(Inf,P),method = 'L-BFGS-B',hessian = T,parallel = list(forward=T),control = list(ndeps=rep(5*1e-7,length(par.ini))))
  
  pn.llh=-m$value
  
 #  parallel::stopCluster(cl)
  
  se.coef = HelpersMG::SEfromHessian(m$hessian)
  se.coef[se.coef<1e-6]<-1e-6
  tval = m$par/se.coef
  p.value=pchisq(tval^2,1,lower.tail = F)
  significance=ifelse(p.value<0.1 & p.value>0.05,'.',ifelse(p.value<0.05 & p.value>0.01,'*',ifelse(p.value<0.01 & p.value>0.001,'**',
                                                                                                   ifelse(p.value<0.001 & p.value>0.0001,'***',ifelse(p.value<0.0001,'****','')))))
  matcoef = cbind(m$par, se.coef, tval, p.value,significance)
  
  dimnames(matcoef) = list(names(par.ini), c("Estimate", "Std. Error", "t value", "Wald P-value"," "))
  
  if(trace){
    print(matcoef,quote = F)
  }
  
  if(is.null(others_abun) & is.null(others_pres)){
    return(list(Main_coef=as.data.frame(matcoef),LogLikelihood=-m$value))
  }else{            
    return(list(Others_coef=as.data.frame(matcoef[1:(P_others_abun+P_others_pres),]),Main_coef=as.data.frame(matcoef[-(1:(P_others_abun+P_others_pres)),]),LogLikelihood=-m$value))
  }
  }  