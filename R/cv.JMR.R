
 
cv.JMR<-function(otu_tab,otu_id,long_design,logistic_design,outcome,long_idset,logistic_idset,rand.var,cov.taxa=T,shrinkage=0.1,n.cores=NULL){

  set_ids=unique(logistic_idset[,2])
  n_set=length(set_ids)
  foldsize=floor(n_set/5)
  
  cv=NULL
  for(i in otu_id){
    cat('Tuning OTU',i,':','\n')
  if(cov.taxa==T){
    mean.abun=rowMeans(otu_tab)
    dist.abun <- as.matrix(vegan::vegdist(otu_tab,method = "bray"))
    dist.pres <- as.matrix(vegan::vegdist(otu_tab,method = "bray",binary = T))
    
    
      dist_abun_i=dist.abun[i,][-i]
      dist_pres_i=dist.pres[i,][-i]
      if(min(dist_abun_i)<0.5){
        limit_abun=min(sum(dist_abun_i<0.5),5)
        id_abun=order(dist.abun[i,],decreasing = F)[1+(1:limit_abun)]
      }else{
        id_abun=NULL
      }
      
      if(min(dist_pres_i)<0.5){
        limit_pres=min(sum(dist_pres_i<0.5),5)
        id_pres=order(dist.pres[i,],decreasing = F)[1+(1:limit_pres)]
      } else{
        id_pres=NULL
      }
      
      llh_t=NULL
      for(f in 1:5){
        
        tfold.set.id=set_ids[((f-1)*foldsize+1):(f*foldsize)]
        fold.smp.id=rownames(long_idset)[!long_idset[,2] %in% tfold.set.id]
        fold.sub.id=rownames(logistic_idset)[!logistic_idset[,2] %in% tfold.set.id]
        tr_otu_tab=otu_tab[,fold.smp.id]
        tr_long_design=long_design[fold.smp.id,]
        tr_logistic_design=logistic_design[fold.sub.id,]
        tr_outcome=outcome[fold.sub.id]
        tr_long_idset=long_idset[fold.smp.id,]
        tr_logistic_idset=logistic_idset[fold.sub.id,]
        
        t_otu_tab=otu_tab[,!colnames(otu_tab) %in% fold.smp.id]
        t_long_design=long_design[!rownames(long_design) %in% fold.smp.id,]
        t_logistic_design=logistic_design[!rownames(logistic_design) %in% fold.sub.id,]
        t_outcome=outcome[!names(outcome) %in% fold.sub.id]
        t_long_idset=long_idset[!rownames(long_idset) %in% fold.smp.id,]
        t_logistic_idset=logistic_idset[!rownames(logistic_idset) %in% fold.sub.id,]
        tr_taxa=tr_otu_tab[i,]
        t_taxa=t_otu_tab[i,]
        
        if(length(id_abun)>1){
          tr_others_abun=t(tr_otu_tab[id_abun,])
          t_others_abun=t(t_otu_tab[id_abun,])
        }else if(length(id_abun)==1){
          tr_others_abun=as.matrix(tr_otu_tab[id_abun,])
          t_others_abun=as.matrix(t_otu_tab[id_abun,])
        }else {
          tr_others_abun=t_others_abun=NULL
        }
        
        if(length(id_pres)>1){
          tr_others_pres=t(tr_otu_tab[id_pres,])
          t_others_pres=t(t_otu_tab[id_pres,])
        }else if(length(id_pres)==1){
          tr_others_pres=as.matrix(tr_otu_tab[id_pres,])
          t_others_pres=as.matrix(t_otu_tab[id_pres,])
        } else {
          tr_others_pres=t_others_pres=NULL
        }
        
        
        joint.res=JMR_core(taxa = tr_taxa,others_abun=tr_others_abun,others_pres=tr_others_pres,long_design = tr_long_design,logistic_design = tr_logistic_design,
                      outcome = tr_outcome,long_idset = tr_long_idset,logistic_idset = tr_logistic_idset,rand.var,shrinkage ,trace = F,n.cores)
        
        if(is.null(tr_others_abun) & is.null(tr_others_pres)){
          all_coef=as.numeric(c(joint.res$Main_coef$Estimate))
        }else{
          all_coef=as.numeric(c(joint.res$Others_coef$Estimate,joint.res$Main_coef$Estimate))
        }
        
        test.llh=-llh_match(par_all = all_coef,t_taxa,t_others_abun,t_others_pres,t_long_design,t_logistic_design,
                            t_outcome,t_long_idset,t_logistic_idset,rand.var,shrinkage)
        
        llh_t=c(llh_t,test.llh)
        cat('Validation Fold',f,';','Penalized likelihood:',test.llh,'\n')
      }
      cv=c(cv,llh_t)
    }
  else{
    
    llh_t=NULL
    for(f in 1:5){
      
      tfold.set.id=set_ids[((f-1)*foldsize+1):(f*foldsize)]
      fold.smp.id=rownames(long_idset)[!long_idset[,2] %in% tfold.set.id]
      fold.sub.id=rownames(logistic_idset)[!logistic_idset[,2] %in% tfold.set.id]
      tr_otu_tab=otu_tab[,fold.smp.id]
      tr_long_design=long_design[fold.smp.id,]
      tr_logistic_design=logistic_design[fold.sub.id,]
      tr_outcome=outcome[fold.sub.id]
      tr_long_idset=long_idset[fold.smp.id,]
      tr_logistic_idset=logistic_idset[fold.sub.id,]
      
      t_otu_tab=otu_tab[,!colnames(otu_tab) %in% fold.smp.id]
      t_long_design=long_design[!rownames(long_design) %in% fold.smp.id,]
      t_logistic_design=logistic_design[!rownames(logistic_design) %in% fold.sub.id,]
      t_outcome=outcome[!names(outcome) %in% fold.sub.id]
      t_long_idset=long_idset[!rownames(long_idset) %in% fold.smp.id,]
      t_logistic_idset=logistic_idset[!rownames(logistic_idset) %in% fold.sub.id,]
      tr_taxa=tr_otu_tab[i,]
      t_taxa=t_otu_tab[i,]
      
      
      tr_others_abun=t_others_abun=NULL
      tr_others_pres=t_others_pres=NULL
      
      joint.res=JMR_core(taxa = tr_taxa,others_abun=tr_others_abun,others_pres=tr_others_pres,long_design = tr_long_design,logistic_design = tr_logistic_design,
                    outcome = tr_outcome,long_idset = tr_long_idset,logistic_idset = tr_logistic_idset,rand.var,shrinkage ,trace = F,n.cores)
      
      if(is.null(tr_others_abun) & is.null(tr_others_pres)){
        all_coef=as.numeric(c(joint.res$Main_coef$Estimate))
      }else{
        all_coef=as.numeric(c(joint.res$Others_coef$Estimate,joint.res$Main_coef$Estimate))
      }
      
      test.llh=-llh_match(par_all = all_coef,t_taxa,t_others_abun,t_others_pres,t_long_design,t_logistic_design,
                          t_outcome,t_long_idset,t_logistic_idset,rand.var,shrinkage)
      
      llh_t=c(llh_t,test.llh)
      cat('Validation Fold',f,';','Penalized likelihood:',test.llh,'\n')
    }
    cv=c(cv,llh_t)
    
  }
  }
  
  return(as.numeric(cv))
}