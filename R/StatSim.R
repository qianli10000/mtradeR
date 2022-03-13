#' @title Simulation for disease status, risk factor(s) and time points in matched sets
#' @param n Number of subjects
#' @param omegas The values of omegas for subject-level and set-level effect on disease
#' @param tps Number of time points.
#' @export

StatSim<-function(n,omegas=c(1.5,1),tps=5){

omega1=omegas[1]
omega2=omegas[2]

sub.id=paste('subject',1:n,sep = '_')
set.id=paste('set',rep(1:(n/2),each=2),sep = '')
sub.order.in.set=rep(1:2,n/2)
#Subject-level predictors
X_sub=data.frame(
  id=sub.id,set=set.id,order=sub.order.in.set)

X_long=data.frame(id=rep(sub.id,each=tps),age=rep(1:tps,times=n))

rownames(X_long)=paste(X_long$id,X_long$age,sep ='_' )
X_all=merge(X_sub,X_long,by = 'id')
rownames(X_all)=paste(X_all$id,X_all$age,sep ='_' )
X_all=X_all[rownames(X_long),]
X=model.matrix(~age,X_all)
rand.var='(Intercept)'
X_rand=X[,rand.var]
X_fix=X
Sigma1=2
Rand_sub=rnorm(n,0,sd = Sigma1)
names(Rand_sub)=sub.id
Rand_sub.cut=cut(Rand_sub,breaks = quantile(Rand_sub,probs = c(0,0.6,1)),include.lowest = T)
names(Rand_sub.cut)=sub.id
Rand_a=Rand_sub[X_all$id]
names(Rand_a)=rownames(X_fix)
Rand_a.cut=cut(Rand_a,breaks = quantile(Rand_a,probs = c(0,0.6,1)),include.lowest = T)
names(Rand_a.cut)=names(Rand_a)

Sigma2=2
Rand_set=rnorm(n/2,0,sd = Sigma2)
names(Rand_set)=unique(set.id)
Rand_set.cut=cut(Rand_set,breaks = quantile(Rand_set,probs = c(0,0.6,1)),include.lowest = T)
names(Rand_set.cut)=unique(set.id)
Rand_b=Rand_set
Rand_b.cut=cut(Rand_b,breaks = quantile(Rand_b,probs = c(0,0.6,1)),include.lowest = T)
names(Rand_b.cut)=names(Rand_b)

Z_all=data.frame(id=sub.id,set=set.id,order=sub.order.in.set,Rand_sub.cut=Rand_sub.cut,Rand_set.cut=Rand_set.cut[set.id],
                 genetic=factor(sample(c(rep(1,n/2),rep(2,n/2)),size = n),levels = 1:2,labels = c('H','L')))
rownames(Z_all)=Z_all$id
Z=model.matrix(~genetic+Rand_sub.cut+Rand_set.cut,Z_all)
Alpha=c(0.5,-2,omega1,omega2)

logit_prob=Z%*%as.matrix(Alpha)
prob=round(exp(logit_prob)/(1+exp(logit_prob)),digits = 2 )

O=rep(0,n)
for(p in unique(prob)){
  O[prob==p]=stats::rbinom(sum(prob==p),size =1 ,prob = p)
}
Z_all$outcome=O

var_all=merge(Z_all,X_all,by=c('id','set','order'),all=F)
var_all$ageset.id=as.factor(var_all$age)
var_all$ageset.id=as.numeric(var_all$ageset.id)
rownames(var_all)=paste(var_all$id,var_all$age,sep = '_')
meta_data=var_all[order(var_all$ageset.id),]

meta_data=subset(meta_data)

return(meta_data)

}
