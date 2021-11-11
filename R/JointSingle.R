#' @title Joint Single random effect model
#' @description This function performs trajectory analysis (either intercept or slope) with disease outcome, not using matched set indicator.
#' @param DataPrep An object returned by DataPrep function.
#' @param filters A numeric vector for taxa filters, specified in the order of c(minimum relative abundance, minimum prevalence )
#' @param identifier A character vector for identifiers in DataPrep$meta_data, in the order of c(set indicator, subject ID).
#' @param disease_status A character for disease outcome variable in DataPrep$meta_data.
#' @param sample_age Name of the covariate representing age (or time point) for each sample.
#' @param disease_cov Names of the covariates in disease sub-model.
#' @param OTU_cov Names of the covariates used in OTU sub-models.
#' @param trajectory_type The type of trajectory analysis: intercept or slope.
#' @param shrinkage Regularization parameter defaults to 0.2.
#' @param trace If trace=FALSE, JointMatch does not print out parameter estimate per taxon. Default to TRUE.
#' @export


JointSingle <- function(DataPrep,filters,identifiers,disease_status,sample_age=NULL,disease_cov,
                       OTU_cov,trajectory_type,shrinkage=0.2,matching_type=NULL,matching_factors=NULL,trace=TRUE)

{
  if(!trajectory_type %in% c('intercept','slope')){
    stop('Error: trajectory type is not correct')
  }

  if(trajectory_type=='slope' & is.null(sample_age)){
    stop('Error: sample_age must be provided if trajectory_type is slope')
  }

  relabun=DataPrep$relabun
  metadata=DataPrep$meta_data

  metadata = metadata[order(metadata[,identifiers[1]],metadata[,identifiers[2]]),]

  id.only=unique(subset(metadata,select = c(identifiers[1],identifiers[2])))
  set.sum=stats::aggregate(id.only[,identifiers[2]],by=list(id.only[,1]),FUN=length)

  sum(set.sum$x==2)
  keep.set=set.sum$Group.1[set.sum$x==2]

  new.longi=subset(metadata, metadata[,identifiers[1]] %in% keep.set)

  clean.id=subset(id.only,id.only[,identifiers[1]] %in% keep.set)
  clean.id$order=rep(1:2,length(keep.set))

  new.longi=merge(new.longi,clean.id,by = c(identifiers))

  if(!is.null(matching_type)){
    if (matching_type=='ncc'){
      disease_cov=disease_cov[!disease_cov %in% matching_factors]
    } else {stop('Matching type was not correct')}
  }
  logistic_only=unique(subset(new.longi,select=c(identifiers,'order',disease_cov,disease_status)))


  longi_only=unique(subset(new.longi,select=c(identifiers,'order',OTU_cov)))

  logistic_formula='~'
  for(var in disease_cov){
    if(which(var==disease_cov)==1)
    {
      logistic_formula=paste(logistic_formula,var,sep = '')
    } else {
      logistic_formula=paste(logistic_formula,var,sep = '+')
    }
  }

  longi_formula='~'
  for(var in OTU_cov){
    if(which(var==OTU_cov)==1)
    {
      longi_formula=paste(longi_formula,var,sep = '')
    } else {
      longi_formula=paste(longi_formula,var,sep = '+')
    }
  }
  #}

  logistic_design=stats::model.matrix(stats::as.formula(logistic_formula),data=logistic_only)

  longi_design=stats::model.matrix(stats::as.formula(longi_formula),data=longi_only)

  logistic_idset=logistic_only[,c(identifiers[2],'order')]

  longi_idset=longi_only[,c(identifiers[2],'order')]

  rand.var=ifelse(trajectory_type=='intercept','(Intercept)',sample_age)

  outcome=logistic_only[,colnames(logistic_only)==disease_status]

  filter=rowMeans(relabun)>filters[1] & rowSums(relabun>0)>filters[2]*ncol(relabun)
  taxa_input=relabun[filter,]


  res=core_single(otu=t(taxa_input),longi_design_all=longi_design,logistic_design_all=logistic_design,outcome,longi_idset,logistic_idset,rand.var,shrinkage,trace)

  return(res)

}
