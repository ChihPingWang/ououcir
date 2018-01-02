rm(list=ls())
#install.packages("TreeSim")
library(TreeSim)
library(EasyABC)
library(coda)
library(Sim.DiffProc)
#oubmbmintegrand<-function(s, alpha.y=alpha.y){
#    alpha.y*exp(alpha.y*s)*rnorm(n=1,mean=0,sd=s)
#    }

n<-20
numbsim<-1;lambda<-2.0;mu<-0.5;frac<-0.6;age<-2
tree<-sim.bd.taxa.age(n=n,numbsim=1,lambda=lambda,mu=mu,frac=frac,age=age,mrca=TRUE)[[1]]
tree<-reorder(tree,"postorder")
tree$edge
plot(tree)
nodelabels()
tiplabels()

root<-0
true.alpha.y<-1
true.alpha.x<-0.2
true.optim.x<-0.5
true.sigma.sq.x<-2
true.tau<- 1

true.b0 <- 0
true.b1 <- 1
true.b2 <- 0.2

##################NOT YET

model<-function(model.params,reg.params,root=root,tree=tree){
  alpha.y<<-model.params[1]
  alpha.x<-model.params[2]
  optim.x<-model.params[3]
  sigma.sq.x<-model.params[4]
  tau<-model.params[5]

  b0<-reg.params[1]
  b1<-reg.params[2]
  b2<-reg.params[3]

  n<-Ntip(tree)
  x1nodestates<-array(0,c(2*n-1))
  x1nodestates[n+1]<-root
  x2nodestates<-array(0,c(2*n-1))
  x2nodestates[n+1]<-root
  optimnodestates<-array(0,c(2*n-1))
  optimnodestates[n+1]<-root
  sigmasqnodestates<-array(0,c(2*n-1))
  sigmasqnodestates[n+1]<-root
  ynodestates<-array(0,c(2*n-1))
  ynodestates[n+1]<-root

  N<-dim(tree$edge)[1]
  anc<-tree$edge[,1]
  des<-tree$edge[,2]
  treelength<-tree$edge.length

  for(index in N:1){
    
    x1ou.mean<-x1nodestates[anc[index]]*exp(-alpha.x*treelength[index]) + optim.x*(1-exp(-alpha.x*treelength[index]))
    x1ou.sd<-sqrt((sigma.sq.x/(2*alpha.x))*(1-exp(-2*alpha.x*treelength[index])))
    x1nodestates[des[index]]<-rnorm(n=1,mean=x1ou.mean,sd=x1ou.sd)

    x2ou.mean<-x2nodestates[anc[index]]*exp(-alpha.x*treelength[index]) + optim.x*(1-exp(-alpha.x*treelength[index]))
    x2ou.sd<-sqrt((sigma.sq.x/(2*alpha.x))*(1-exp(-2*alpha.x*treelength[index])))
    x2nodestates[des[index]]<-rnorm(n=1,mean=x2ou.mean,sd=x2ou.sd)

    optimnodestates[des[index]]<- b0 + b1*x1nodestates[des[index]] + b2*x2nodestates[des[index]]

    sigmasqnodestates[des[index]]<-rnorm(n=1,mean=sigmasqnodestates[anc[index]],sd=sqrt(tau*treelength[index]))
    sigma.sq.theta<- b1^2*sigma.sq.x + b2^2*sigma.sq.x
    #inttheta<-integrate(oubmbmintegrand,lower=0 ,upper=treelength[index], alpha.y=true.alpha.y)

    A1<-(alpha.y*optimnodestates[des[index]]/ (alpha.y-alpha.x)) *(exp((alpha.y-alpha.x)*treelength[index]) -1)
    theta1<-0
    A2<- theta1*(exp(alpha.y*treelength[index])-1) - (alpha.y*theta1/(alpha.y-alpha.x))*(exp((alpha.y-alpha.x)*treelength[index]) -1)
    A3<-1/(sigma.sq.theta*alpha.y^2)*rnorm(1,mean=0,sd=sqrt((sigma.sq.theta*alpha.y^2)^3 *treelength[index]^3 /3))
    INTtime<- A1+A2+A3
    INT1<-exp(-alpha.y*treelength[index])*INTtime

    fexpr<-expression(alpha.y*exp(alpha.y*t)*w)
    res<-st.int(fexpr,type="ito",M=1,lower=0,upper=treelength[index])
    INT2<-median(res$X)

    ynodestates[des[index]]<-ynodestates[anc[index]] + INT1 + INT2
    }

  simtrait<-ynodestates[1:n]
  return(list(y=simtrait,x1=x1nodestates[1:n],x2=x2nodestates[1:n]))
  #return(c(mean(simtrait),sd(simtrait)))
  }

print(model(model.params=c(5,1,1,3,2),reg.params=c(1,1,1),root=root,tree=tree))

#model.params=c(5,1,1,3,2);reg.params=c(1,1,1)
