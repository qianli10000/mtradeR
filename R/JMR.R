
#' @title Joint model with Matching and Regularization (JMR)
#' @description This function test trajectory (either intercept or slope) association with disease outcome in matched sets.
#' @param otu_tab A table of relative abundance with rows as OTUs and columns as samples. 
#' @param long_design Design matrix for longitudinal variables in the abundance-presence sub-model.
#' @param logistic_design Design matrix for risk factors in the disease sub-model.
#' @param outcome A vector of disease outcome per subject.
#' @param long_idset A dataframe of subject and set identifiers mapping to long_design, in the order of subjectID, setID, order (within-set indicator).
#' @param logistic_idset A dataframe of subject and set identifier mapping to logistic_design, in the order of subjectID, setID, order (within-set indicator).
#' @param rand.var The type of trajectory analysis, intercept: '(Intercept)', or slope: the variable for time point or age in long_design.
#' @param tune A scalar or vector of tuning parameter for L2 regularization. If otu_tab contains <10 rows (OTUs), tune must be a scalar. 
#' @param cov.taxa Whether to adjust for interdependence between taxa, default as cov.taxa=TRUE. If otu_tab contains only one row (OTU), cov.taxa must be FALSE. 
#' @param n.cores Number of workers registered in parallel computing. 
#' @return  \item{$test.result}{The result of joint test on relative abundance and presence.}
#'          \item{$rho}{Tuning parameter value.}
#' @import optimParallel
#' @export
#' @examples
#' data("DM_MLE")
#' 
#' #Generate set indicator and disease outcome
#' meta_data<-StatSim(n=4,tps=2)
#' meta_data<-meta_data[order(meta_data$set,meta_data$id),]
#' subj_data<-unique(meta_data[,c('id','outcome')])
#' outcome<-subj_data$outcome
#' names(outcome)<-subj_data$id
#' long_design <- model.matrix(~age,meta_data)
#' logistic<-unique(subset(meta_data,select = -c(age,ageset.id)))
#' rownames(logistic)=logistic$id
#' logistic_design <- model.matrix(~genetic,logistic)
#' long_idset <- meta_data[,c('id','set','order')]
#' logistic_idset <- logistic[,c('id','set','order')]
#' 
#' #Generate a toy example of metagenomic raw counts table with dimension P=20.
#' raw.counts=TaxaSim(DM_MLE[1:20],StatSim = meta_data,shift_subject = 0,trace =FALSE)
#' rel.abun=t(t(raw.counts)/colSums(raw.counts))
#' 
#' #Filter taxa by relative abundance and prevalence
#' mean.rel.abun=rowMeans(rel.abun)
#' filter=mean.rel.abun>1e-6 & rowSums(rel.abun==0)<0.95*ncol(rel.abun)
#' input_tab=rel.abun[filter,]
#' 
#' #Run JMR without tuning and covariate taxa for the first OTU. 
#' JMR.res=JMR(otu_tab = input_tab[1,],long_design = long_design,logistic_design = logistic_design,
#' outcome = outcome, long_idset = long_idset,logistic_idset = logistic_idset,rand.var = '(Intercept)',
#' tune=0.1,cov.taxa=FALSE,n.cores=1)



JMR<-function(otu_tab,long_design,logistic_design,outcome,long_idset,logistic_idset,rand.var,
                  tune=0.15,cov.taxa=T,n.cores=NULL){


  
if(is.null(nrow(otu_tab))){
  otu_tab=matrix(otu_tab,nrow=1)
}
  otu_tab=as.matrix(otu_tab)
  mean.abun=rowMeans(otu_tab)
 
   if(is.null(n.cores)){
    n.cores=min(parallel::detectCores()-1,16)
  }
  
  cl <- parallel::makeCluster(n.cores)     # set the number of processor cores

  parallel::setDefaultCluster(cl=cl)
 
# cross-validation
if(length(tune)>1){
  filter=which(mean.abun>0.1)
  if(length(filter)>0){
  otu_id1=sample(filter,size = 1)
  } else{otu_id1=NULL}
  filter=which(mean.abun>0.01 & mean.abun<0.1)
  if(length(filter)>0){
  otu_id2=sample(filter,size = 1)
  }else{otu_id2=NULL}
  filter=which(mean.abun>0.001 & mean.abun<0.01)
  otu_id3=sample(filter,size = 1)
  filter=which(mean.abun>0.0001 & mean.abun<0.001)
  otu_id4=sample(filter,size = 1)
  filter=which(mean.abun>0.00001 & mean.abun<0.0001)
  otu_id5=sample(filter,size = 1)
  otu_id=c(otu_id1,otu_id2,otu_id3,otu_id4,otu_id5)

  cv=NULL
  for(rho in tune){
    cat('shrinkage=',rho,'\n')
    cv_val=cv.JMR(otu_tab,otu_id,long_design,logistic_design,outcome,long_idset,logistic_idset,rand.var,cov.taxa,shrinkage=rho,n.cores)
    cv=rbind(cv,cv_val)
  }

  mean.pllk=rowMeans(cv)
  names(mean.pllk)=tune
  cvg=mean.pllk[-1]/mean.pllk[-length(mean.pllk)]
  if(max(cvg>1)){
   rho_opt=tune[mean.pllk==max(mean.pllk)]
   shrinkage=rho_opt
  }else if(sum(cvg>0.98 & cvg<1)>0){
    sl=min(which(cvg>0.98 & cvg<1))
    rho_opt=as.numeric(names(cvg)[sl])
    shrinkage=rho_opt
  }else{shrinkage=max(tune)}
}else{shrinkage=tune}
cat('shrinkage=',shrinkage,'\n')



n_otu=nrow(otu_tab)
res_lambda=NULL
for(i in 1:n_otu){
  if(cov.taxa==F){
    limit_abun=0
    others_abun=NULL
    limit_pres=0
    others_pres=NULL
  }else{ 
  
  dist.abun <- as.matrix(vegan::vegdist(otu_tab,method = "bray"))
  dist.pres <- as.matrix(vegan::vegdist(otu_tab,method = "bray",binary = T))
  d_abun=quantile(dist.abun[dist.abun!=0],probs = 0.1)
  d_pres=quantile(dist.pres[dist.pres!=0],probs = 0.1)  
  dist_abun_i=dist.abun[i,][-i]
  dist_pres_i=dist.pres[i,][-i]
  if(min(dist_abun_i)<d_abun){
  limit_abun=min(25,sum(dist_abun_i<d_abun))
  id_abun=order(dist.abun[i,],decreasing = F)[1+(1:limit_abun)]
    
  if(length(id_abun)>1){
    others_abun=t(otu_tab[id_abun,])
  }else{others_abun=as.matrix(otu_tab[id_abun,])}
  }else{
    limit_abun=0
    others_abun=NULL
  }
  
  if(min(dist_pres_i)<d_pres){
    limit_pres=min(25,sum(dist_pres_i<d_pres))
    id_pres=order(dist.pres[i,],decreasing = F)[1+(1:limit_pres)]
    if(length(id_pres)>1){
    others_pres=t(otu_tab[id_pres,])
    }else{others_pres=as.matrix(otu_tab[id_pres,])}
  } else{
    limit_pres=0
    others_pres=NULL
  }
  
  }
  
  taxa=otu_tab[i,]
  cat('Taxon', i, ':',rownames(otu_tab)[i],'Correlated taxa:',limit_abun,',',limit_pres,'\n')
  cat('Mean relative abundance:',mean(taxa),'\t')
  cat('Prevalence:',sum(taxa>0)/length(taxa),'\n')
  
  if(limit_abun>4){
  taxa_tr=asin(sqrt(taxa))
  set.seed(321)
  cvfit.abun = glmnet::cv.glmnet(others_abun, taxa_tr,alpha=0.05,nfolds = 10,family='gaussian')
 
  select=cvfit.abun$lambda==cvfit.abun$lambda.min
  beta_abun=cvfit.abun$glmnet.fit$beta[,select]
  others_abun_input=as.matrix(others_abun[,beta_abun!=0])
  if(ncol(others_abun_input)==0){others_abun_input=NULL}
  }else{others_abun_input=others_abun}

  if(limit_pres>4){
  if(sum(taxa==0)>0.02*length(taxa)){
  taxa_bin=ifelse(taxa==0,0,1)
  others_pres=ifelse(others_pres==0,0,1)
  
  if(sum(rowSums(others_pres==0)>0)>0.1*nrow(others_pres)){
  set.seed(321)
  cvfit.pres = glmnet::cv.glmnet(others_pres, taxa_bin,alpha=0.05,nfolds = 10,family='binomial')
  
  select=cvfit.pres$lambda==cvfit.pres$lambda.min
  beta_pres=cvfit.pres$glmnet.fit$beta[,select]
  others_pres_input=as.matrix(others_pres[,beta_pres!=0])
  if(ncol(others_pres_input)==0){others_pres_input=NULL}
  }else{others_pres_input=NULL}
  }else{others_pres_input=NULL}
  }else{others_pres_input=others_pres}

  
  res=JMR_core(taxa,others_abun_input,others_pres_input,long_design,logistic_design,outcome,long_idset,logistic_idset,rand.var,shrinkage,trace=F)
  
  res_lambda=rbind(res_lambda,as.numeric(c(res$Main_coef['composition_lambda',-5],res$Main_coef['presence_lambda',-5])))

  }

parallel::stopCluster(cl)

colnames(res_lambda)=c('rabun_coef','rabun_se','rabun_t','rabun_p','pres_coef','pres_se','pres_t','pres_p')
rownames(res_lambda)=rownames(otu_tab)
joint.pvalue=pchisq(as.numeric(res_lambda[,'rabun_t'])^2+as.numeric(res_lambda[,'pres_t'])^2,df = 2,lower.tail = F)
FDR=p.adjust(joint.pvalue,method = 'BH')
res_lambda=as.data.frame(res_lambda)
res_lambda$joint.pvalue=joint.pvalue
res_lambda$FDR=FDR
output=list(test.result=res_lambda,rho=shrinkage)
return(output)

}