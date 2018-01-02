rm(list=ls())
#install.packages("TreeSim")
library(TreeSim)
library(EasyABC)
library(coda)
n<-30
numbsim<-1;lambda<-2.0;mu<-0.5;frac<-0.6;age<-2
#?TreeSim
tree<-sim.bd.taxa.age(n=n,numbsim=1,lambda=lambda,mu=mu,frac=frac,age=age,mrca=TRUE)[[1]]
tree<-reorder(tree,"postorder")
tree$edge
plot(tree)
nodelabels()
tiplabels()

root<-0
true.sigma.sq<-1
nodestates<-array(0,c(2*n-1))
nodestates[n+1]<-root
N<-dim(tree$edge)[1]
anc<-tree$edge[,1]
des<-tree$edge[,2]
treelength<-tree$edge.length
for(index in N:1){

#  print(c(anc[index],des[index]))
  nodestates[des[index]]<-rnorm(n=1,mean=nodestates[anc[index]],sd= sqrt(true.sigma.sq*treelength[index]))
  }
tipstates<-nodestates[1:n]
print(tipstates)

truetrait<-tipstates

s1=mean(truetrait)
s2=min(truetrait)
s3=max(truetrait)
s4=median(truetrait)
s5=sd(truetrait)
s6=IQR = summary(truetrait)[4]-summary(truetrait)[2]
s7=MAD = mad(truetrait)

c(s1,s2,s3,s4,s5,s6,s7)/mad(c(s1,s2,s3,s4,s5,s6,s7))

x=c(1,1,2,2,4,6,9)
median(x)
median(abs(x-median(x)))
summarytrait=c(mean(truetrait),sd(truetrait))



model<-function(sigmasq,root=root,tree=tree){
  #sigmasq<-1
  n<-Ntip(tree)
  nodestates<-array(0,c(2*n-1))
  nodestates[n+1]<-root
  N<-dim(tree$edge)[1]
  anc<-tree$edge[,1]
  des<-tree$edge[,2]
  treelength<-tree$edge.length
  for(index in N:1){
  nodestates[des[index]]<-rnorm(n=1,mean=nodestates[anc[index]],sd=sqrt(sigmasq*treelength[index]))
  }
  simtrait<-nodestates[1:n]
  return(c(mean(simtrait),sd(simtrait)))
  }
#model(1,root=root,tree=tree)
# sigmasqprior<- list(c("unif",0,100))
# ABC_Marjoram_original<-ABC_mcmc(method="Marjoram",model=model,root=root,tree=tree,prior=sigmasqprior,summary_stat_target=summarytrait,n_rec=10000)
# str(ABC_Marjoram_original)
# hist(ABC_Marjoram_original$param[5000:10000],main="Posterior of sigma.sq")
#cannot pass through argument
ABC_acc<-function(sigmasq,root=root,tree=tree,summarytrait=summarytrait){
    if(sigmasq<=0){return(FALSE)}
    summarysamples<-model(sigmasq,root=root,tree=tree)
    diffmean<-abs(summarysamples[1]-summarytrait[1])
    diffsd<-abs(summarysamples[2]-summarytrait[2])
    if( (diffmean<0.1) & (diffsd<0.2)){
      return(TRUE)}else{return(FALSE)}
      }


run_MCMC_ABC<-function(startvalue,iterations,root=root,tree=tree,summarytrait=summarytrait){
  chain<-array(0,dim=c(iterations+1,1))
  chain[1,]<-startvalue
  for(i in 1:iterations){
    proposal<-rnorm(n=1,mean=chain[i,],sd=c(1))# it has to be a prior distribution, may not be proposal
    if(ABC_acc(proposal, root=root,tree=tree,summarytrait=summarytrait)){
    chain[i+1,]<-proposal
    }else{
    chain[i+1,]<-chain[i,]
    }
  }
  return(mcmc(chain))
  }

posterior<-run_MCMC_ABC(c(1),50000,root=root,tree=tree,summarytrait=summarytrait)
plot(posterior)
mean(unique(posterior))
hist(unique(posterior))
boxplot(unique(posterior))
print(mean(posterior))


#https://theoreticalecology.wordpress.com/2012/07/15/a-simple-approximate-bayesian-computation-mcmc-abc-mcmc-in-r/
# data<-rnorm(10,mean=5.3,sd=2.7)
# meandata<-mean(data)
# sddata<-sd(data)
#
# ABC_acc<-function(par){
#   if(par[2]<=0){return(FALSE)}
#   samples<-rnorm(10,mean=par[1],sd=par[2])
#   diffmean<- abs(mean(samples)-meandata)
#   diffsd<-abs(sd(samples)-sddata)
#   if((diffmean<0.1) & (diffsd<0.2)){
#     return(T)}else{
#       return(F)}
#   }
#
# run_MCMC_ABC<-function(startvalue,iterations){
#   chain<-array(0,dim=c(iterations+1,2))
#   chain[1,]=startvalue
#   for(i in 1:iterations){
#     proposal <- rnorm(2,mean=chain[i,], sd=c(0.7,0.7))
#  if(ABC_acc(proposal)){
#    chain[i+1,]=proposal
#     }else{
#       chain[i+1,]=chain[i,]
#       }
#     }
#   return(mcmc(chain))
#   }
#
# posterior<-run_MCMC_ABC(c(4,2.3),300000)
#
# plot(posterior)




# The following can only work for model without argument
#https://theoreticalecology.wordpress.com/2012/12/02/the-easyabc-package-for-approximate-bayesian-computation-in-r/
 # data<-rnorm(10,mean=5.3,sd=2.7)
 # summarydata<-c(mean(data),sd(data))
 #
 # model<-function(par){
 #   samples<-rnorm(10,mean=par[1],sd=par[2])
 #   #print(hi)
 #   return(c(mean(samples),sd(samples)))
 # }
 # ABC_Marjoram_original<- ABC_mcmc(method="Marjoram_original",model=model, prior=list(c("unif",0,10),c("unif",1,5)),summary_stat_target=summarydata, n_rec=10000)
 # str(ABC_Marjoram_original)
 # par(mfrow=c(2,1))
 # hist(ABC_Marjoram_original$param[5000:10000,1],main="Posterior for mean")
 # hist(ABC_Marjoram_original$param[5000:10000,2],main="Posterior for standard deviation")


#rm(list=ls())
#install.packages("EasyABC")
# library(EasyABC)
#my_prior=list(c("unif",0,1),c("normal",1,2))
#my_prior=list(list(c("runif",1,0,1), c("dunif",0,1)))
# toy_model<-function(x){
#   c(x[1]+x[2]+rnorm(1,0,0.1),x[1]*x[2]+rnorm(1,0,0.1))
# }
# toy_prior = list(c("unif",0,100),c("normal",1,2))
# sum_stat_obs=c(1.5,0.5)
# set.seed(1)
# n=10
# p=0.5
# ABC_rej<-ABC_rejection(model=toy_model,prior=toy_prior,nb_simul=n,summary_stat_target=sum_stat_obs,tol=p)
