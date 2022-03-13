# mtradeR 
Metagenomic TRajectory Analysis with Disease Endpoint and Risk factors. 

We proposed a joint mixed effect model (JointMatch) to detect OTUs predictive of hosts’ disease status in matched sets. We also designed and implemented a semi-parametric simulation pipeline was built to generate disease outcome and temporal high-dimensional metagenomic counts in matched sets. 

# An example

#Below is an exmaple code for running the simulation pipeline and test JointMatch:

install.packages("devtools")

library(devtools)

install_github('qianli10000/mtradeR')

library('mtradeR')

data("DM_MLE")

meta_data<-meta_data[order(meta_data$set,meta_data$id),]

outcome<-meta_data$outcome

names(outcome)<-meta_data$id

long_design <- model.matrix(~age,meta_data)

logistic<-unique(subset(meta_data,select = -c(age,ageset.id)) )

logistic_design <- model.matrix(~genetic,logistic)

long_idset <- meta_data[,c('id','set','order')]

logistic_idset <- logistic[,c('id','set','order')]

JMR.res=JMR_Tab(otu_tab = taxa_input,long_design = long_design,logistic_design = logistic_design,outcome = outcome,
                          long_idset = long_idset,logistic_idset = logistic_idset,rand.var = '(Intercept)',
                          tune=seq(0.05,0.15,0.05),cov.taxa=T)
