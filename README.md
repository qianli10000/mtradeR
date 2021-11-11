# mtradeR: metagenomic trajectory analysis with disease endpoint and risk factors
# Below is an exmaple code for running simulation pipeline and test JointMatch:

data("DM_MLE")

meta_data<-StatSim(n=150)

otu_counts<-TaxaSim(base_par = DM_MLE,StatSim = meta_data)
rel_abun<-t(t(otu_counts)/colSums(otu_counts))

#Use the first 100 taxa to test 
dp<-list(relabun=rel_abun[1:100,],meta_data=meta_data)

head(dp$meta_data)

##                    id  set  Rand_sub.cut  Rand_set.cut genetic outcome age ageset.id
##subject_1_1   subject_1 set1 [-4.49,0.516] [-4.43,0.333]       L       0   1         1
##subject_10_1 subject_10 set5 [-4.49,0.516] [-4.43,0.333]       H       1   1         1


taxa_filters=c(0.000001,0.05)
identifiers=c('set','id')
disease_s='outcome'
disease_c='genetic'
otu_c=c('age','genetic')


test.run<-JointMatch(DataPrep = dp,filters = taxa_filters,identifiers = identifiers,disease_status = disease_s,
disease_cov = disease_c,OTU_cov =otu_c,trajectory_type = 'intercept' )
