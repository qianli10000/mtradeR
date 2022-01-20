
#' @title Joint nested random effect model with Matched sets indicator
#' @description This function performs trajectory analysis (either intercept or slope) with disease outcome in matched sets.
#' @param DataPrep A list of object returned by DataPrep function or prepared by users. This object must contain $relabun: a relative abundance
#'                 matrix or data frame with rows as taxa and columns as samples; $meta_data: a matrix or data frame of subject ID, set ID, sample ID, disease outcome, time-variant and time-invariant covariates for all submodels,
#'                 with columns as variables, and rows as samples in the order of $relabun's samples (columns).
#' @param filters A numeric vector for taxa filters, specified in the order of c(minimum relative abundance, minimum prevalence )
#' @param identifiers A character vector for set and subject identifiers in DataPrep$meta_data, in the order of c(set ID, subject ID)
#' @param disease_status A character for disease outcome variable in DataPrep$meta_data.
#' @param sample_age Name of the covariate representing age (or time point) for each sample.
#' @param disease_cov Names of the covariates in disease sub-model.
#' @param OTU_cov Names of the covariates used in OTU sub-models.
#' @param trajectory_type The type of trajectory analysis: intercept or slope.
#' @param shrinkage Regularization parameter defaults to 0.2.
#' @param matching_type If participants are matched by a nested case-control design, set matching_type='ncc'.Otherwise, matching_type is NULL as default.
#' @param matching_factors If matching_type='ncc', matching_factors must be provided.
#' @param trace If trace=FALSE, JointMatch does not print out parameter estimate per taxon. Default to TRUE.
#' @return  \item{$OutputPrt}{A list of analysis result for each submodel per taxon. Rows annotated as 'nzAbundance:' is the result for the submodel of non-zero abundance;
#'          rows annotated as 'Presence:' is the result for the submodel of presence; rows annotated as 'Disease:' is the result for the submodel of disease risk.}
#'          \item{$IndividualTest}{The test of taxon-disease association in non-zero abundance and presence, individually}
#'          \item{$JointTest}{Jointly test if a taxon is associated with disease outcome in either non-zero abundance or presence}
#' @export


JointMatch <- function(DataPrep,filters,identifiers,disease_status,sample_age=NULL,disease_cov,
                       OTU_cov,trajectory_type,shrinkage=seq(0.17,0.22,0.01),matching_type=NULL,matching_factors=NULL,trace=TRUE)

{
  if(!trajectory_type %in% c('intercept','slope')){
    stop('Error: trajectory type is not correct')
  }

  if(trajectory_type=='slope' & is.null(sample_age)){
    stop('Error: sample_age must be provided if trajectory_type is slope')
  }


  relabun=DataPrep$relabun
  metadata=DataPrep$meta_data

  relabun=relabun[,order(metadata[,identifiers[1]],metadata[,identifiers[2]])]
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

  sum.rm=sum(!filter)
  cat(sum.rm, 'of',nrow(relabun), 'taxa removed by filters','\n')

  res=core_match(otu=t(taxa_input),longi_design_all=longi_design,logistic_design_all=logistic_design,outcome,longi_idset,logistic_idset,rand.var,shrinkage,trace)

  return(res)

}
