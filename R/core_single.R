#' @export

#############################################################################################
#############################################################################################

# pre-defined likelihood function

core_single<-function(otu,longi_design_all,logistic_design_all,outcome,longi_idset,logistic_idset,rand.var,shrinkage,trace=T){
  options(stringsAsFactors = F)
  require(optimParallel)
  require(mvQuad)
  require(HelpersMG)

  longi_dim<<-ncol(longi_design_all)
  logistic_dim<<-ncol(logistic_design_all)
  P<<-2*longi_dim+logistic_dim+4

  Output.list=vector(mode='list',length=ncol(otu))
  Test.table=matrix(nrow = ncol(otu),ncol=8)


  cl <-parallel::makeCluster(min(16,parallel::detectCores()-1))     # set the number of processor cores
  parallel::clusterExport(cl, list( "longi_dim", "logistic_dim","llh_single", "lh1_single", "lh2_single"))

  parallel::setDefaultCluster(cl=cl) # set 'cl' as default cluster

  for(i in 1:ncol(otu)){
    taxa=otu[,i]
    cat('Taxon', i, ':',colnames(otu)[i],'\n')
    cat('Mean relative abundance:',mean(taxa),'\t')
    cat('Prevalence:',sum(taxa>0)/length(taxa),'\n')

    par.ini=c(composition_cov=rep(0,longi_dim),composition_lambda=0,presence_cov=rep(0,longi_dim),presence_lambda=0,outcome_cov=rep(0,logistic_dim),dispersion=1,variance_a=1)
    names(par.ini)=c(paste('nzAbundance',colnames(longi_design_all),sep =':'), c('nzAbundance_sbj'),
                     paste('Presence',colnames(longi_design_all),sep = ':'), c('Presence_sbj'),
                     paste('Disease',colnames(logistic_design_all),sep = ':'), c('Overdispersion','Variance_sbj'))

    m=optimParallel::optimParallel(par.ini,llh_single,taxa=taxa,longi_design_all=longi_design_all,logistic_design_all=logistic_design_all,outcome=outcome,longi_idset=longi_idset,logistic_idset=logistic_idset,rand.var=rand.var,shrinkage=shrinkage,
                                   lower = c(rep(-Inf,(P-2)),1e-04,1e-04),upper=rep(Inf,P),method = 'L-BFGS-B',hessian = T,parallel = list(forward=T),control = list(ndeps=rep(5*1e-7,length(par.ini))))


    se.coef = HelpersMG::SEfromHessian(m$hessian)
    se.coef[se.coef<1e-4]<-1e-4
    tval = m$par/se.coef
    p.value=stats::pchisq(tval^2,1,lower.tail = F)
    significance=ifelse(p.value<0.1 & p.value>0.05,'.',ifelse(p.value<0.05 & p.value>0.01,'*',ifelse(p.value<0.01 & p.value>0.001,'**',
                                                                                                     ifelse(p.value<0.001 & p.value>0.0001,'***',ifelse(p.value<0.0001,'****','')))))
    matcoef = cbind(m$par, se.coef, tval, p.value,significance)

    dimnames(matcoef) = list(names(par.ini), c("Estimate", "Std. Error", "t value", "p-value"," "))

    if(trace){
      cat(colnames(otu)[i],'\n')
      print(matcoef,quote = F)
    }
    Output=as.data.frame(matcoef)
    Output.list[[i]]=Output
    Test.table[i,]=c(matcoef['nzAbundance_sbj',1:4],matcoef['Presence_sbj',1:4])
  }
  parallel::stopCluster(cl)
  names(Output.list)=colnames(otu)
  Test.table=as.data.frame(Test.table)
  colnames(Test.table)[1:4]=paste('nzAbun', c("Estimate", "Std. Error", "t value", "p-value"),sep = '_')
  colnames(Test.table)[5:8]=paste('Pres',c("Estimate", "Std. Error", "t value", "p-value"),sep = '_')
  rownames(Test.table)=colnames(otu)

  joint.stat=as.numeric(Test.table[,3])^2+as.numeric(Test.table[,7])^2
  joint.p=stats::pchisq(joint.stat, df = 2,lower.tail = F)
  joint.fdr=stats::p.adjust(joint.p,method = 'BH')
  JointTest=data.frame(W.stat=joint.stat,W.pvalue=joint.p,W.adj.pvalue=joint.fdr)
  rownames(JointTest)=colnames(otu)

  return(list(OutputPrt=Output.list,IndividualTest=Test.table,JointTest=JointTest))

}
