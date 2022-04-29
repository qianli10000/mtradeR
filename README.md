# mtradeR 
Metagenomic TRajectory Analysis with Disease Endpoint and Risk factors. 

We proposed a Joint model with Matching and Regularization (JMR) to detect OTUs predictive of hosts’ disease status in matched sets. We also designed and implemented a simulation pipeline to generate disease outcome and temporal high-dimensional metagenomic counts in matched sets. A preprint is available at https://www.biorxiv.org/content/10.1101/2022.04.19.488854v1.full.

# An example

#Below is an exmaple code for running the simulation pipeline and test JointMatch:

install.packages("devtools")

library(devtools)

install_github('qianli10000/mtradeR')

library('mtradeR')

data("DM_MLE")

#Generate set indicator and disease outcome

meta_data<-StatSim(n=150)

meta_data<-meta_data[order(meta_data$set,meta_data$id),]

outcome<-meta_data$outcome

names(outcome)<-meta_data$id

long_design <- model.matrix(~age,meta_data)

logistic<-unique(subset(meta_data,select = -c(age,ageset.id)))

rownames(logistic)=logistic$id

logistic_design <- model.matrix(~genetic,logistic)

long_idset <- meta_data[,c('id','set','order')]

logistic_idset <- logistic[,c('id','set','order')]

#Generate metagenomic raw counts table 

raw.counts=TaxaSim(DM_MLE,StatSim = meta_data,shift_subject = 0,trace =F)

rel.abun=t(t(raw.counts)/colSums(raw.counts))

#Filter taxa by relative abundance and prevalence

mean.rel.abun=rowMeans(rel.abun)

filter=mean.rel.abun>1e-6 & rowSums(rel.abun==0)<0.95*ncol(rel.abun)

input_tab=rel.abun[filter,]

#Run JMR with tuning

JMR.res=JMR_Tab(otu_tab = input_tab,long_design = long_design,logistic_design = logistic_design,outcome = outcome,
                          long_idset = long_idset,logistic_idset = logistic_idset,rand.var = '(Intercept)',
                          tune=seq(0.05,0.15,0.05),cov.taxa=T)
