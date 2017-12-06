#OUBMCIR
rm(list=ls())
library(ape)
library(phytools)
library(corpcor)
library(stringi)
library(MCMCpack)
##given tree of n taxa,response,predictor.dataset,interaction,model
##output everything
#--------------------------------------------------------------------
CovRes<-function(x,n=n,G=G,distG=distG,one=one,response=response,pred.dataset=pred.dataset){
  alp<-x[1]
  bet<-x[2]
  th_1<-x[3]
  xi<-x[4]
  y_0<-0
  th_0<-0
  sig_0<-0
  sig1<-compute.optimal.sigma(model.params = x , response=response, pred.dataset=pred.dataset,G=G,opt.model=NULL)
  cov_yi_yj<-array(0,c(n,n))
  cov_y_th<-array(0,c(n,n))
  cov_thi_thj<-array(0,c(n,n))
  b1_hat<-array(0,c(n,n))
  cov_ri_rj<-array(0,c(n,n))
  for(i in 1:n){
    for(j in 1:n){



    if(bet!=2*alp){
    cov_yi_yj[i,j] <-
    sig1^2*G[i,j]*(exp(1/2*alp*distG[i,j]) - 1)^2*exp(-alp*distG[i,j]) -
    2*((th_0*(exp(alp*G[i,j]) - 1) + y_0)*th_0*exp(-alp*G[i,j]) -
    ((alp*exp(alp*G[i,j]) - alp)*th_0^2 + alp*th_0*y_0 + sig1^2 +
    (alp*sig1^2*G[i,j] - sig1^2)*exp(alp*G[i,j]))*exp(-alp*G[i,j])/alp)*(exp(1/2*alp*distG[i,j]) - 1) -
    1/4*(4*(th_0*(exp(alp*G[i,j]) - 1) + y_0)^2*exp(-2*alp*G[i,j]) -
    (4*(2*(2*alp^2*bet^2 - alp*bet^3)*th_1^2*exp(bet*G[i,j]) -
    (4*alp^3*bet - alp*bet^3)*th_1*exp(bet*G[i,j]) + (2*alp^3*bet -
    alp^2*bet^2)*exp(bet*G[i,j]))*y_0^2*exp(2*bet*G[i,j]*th_1) -
    8*(2*(2*alp^2*bet^2 - alp*bet^3 - (2*alp^2*bet^2 -
    alp*bet^3)*exp(alp*G[i,j]))*th_0*th_1^2*exp(bet*G[i,j]) - (4*alp^3*bet
    - alp*bet^3 - (4*alp^3*bet -
    alp*bet^3)*exp(alp*G[i,j]))*th_0*th_1*exp(bet*G[i,j]) + (2*alp^3*bet -
    alp^2*bet^2 - (2*alp^3*bet -
    alp^2*bet^2)*exp(alp*G[i,j]))*th_0*exp(bet*G[i,j]))*y_0*exp(2*bet*G[i,j]*th_1
    ) - 2*(2*alp^2*bet - alp*bet^2 + 2*(2*alp^2*bet -
    alp*bet^2)*sig_0^2 - 2*(2*alp^2*bet -
    alp*bet^2)*sig_0)*th_1*exp(2*alp*G[i,j] + bet*G[i,j]) -
    (4*alp^2*sig_0*xi^2*exp(2*alp*G[i,j]) +
    4*(2*alp*bet^2*exp(2*alp*G[i,j]) - (bet^3 + (2*alp*bet^2 -
    bet^3)*exp(2*alp*G[i,j]))*exp(bet*G[i,j]))*th_1^3 - 4*(2*alp^3*bet -
    alp^2*bet^2 + (2*alp^3*bet - alp^2*bet^2)*exp(2*alp*G[i,j]) -
    2*(2*alp^3*bet - alp^2*bet^2)*exp(alp*G[i,j]))*th_0^2*exp(bet*G[i,j])
    - 2*(4*(2*alp^2*bet^2 - alp*bet^3 + (2*alp^2*bet^2 -
    alp*bet^3)*exp(2*alp*G[i,j]) - 2*(2*alp^2*bet^2 -
    alp*bet^3)*exp(alp*G[i,j]))*th_0^2*exp(bet*G[i,j]) +
    2*(2*alp*bet^2*sig_0 - alp*bet*xi^2 +
    2*alp^2*bet)*exp(2*alp*G[i,j]) - (4*alp*bet^2*sig_0 -
    bet^2*xi^2 - 8*(2*alp*bet^2 - bet^3)*sig1^2*exp(alp*G[i,j]) +
    bet^3 + 2*(2*alp*bet^2 - bet^3)*sig1^2 - (4*(2*alp^2*bet^2
    - alp*bet^3)*sig1^2*G[i,j] - 4*alp^2*bet + bet^3 -
    6*(2*alp*bet^2 - bet^3)*sig1^2 + (2*alp*bet -
    bet^2)*xi^2)*exp(2*alp*G[i,j]))*exp(bet*G[i,j]))*th_1^2 +
    (4*(4*alp^3*bet - alp*bet^3 + (4*alp^3*bet -
    alp*bet^3)*exp(2*alp*G[i,j]) - 2*(4*alp^3*bet -
    alp*bet^3)*exp(alp*G[i,j]))*th_0^2*exp(bet*G[i,j]) +
    4*(2*alp^2*bet*sig_0 - (alp*bet*sig_0 +
    alp^2)*xi^2)*exp(2*alp*G[i,j]) - (4*alp*bet^2*sig_0 -
    8*(4*alp^2*bet - bet^3)*sig1^2*exp(alp*G[i,j]) + 4*(2*alp^2*bet
    - alp*bet^2)*sig_0^2 + 2*(4*alp^2*bet - bet^3)*sig1^2 -
    (4*alp*bet*sig_0 + bet^2)*xi^2 - (4*(4*alp^3*bet -
    alp*bet^3)*sig1^2*G[i,j] - 4*alp^2*bet + 2*alp*bet^2 -
    6*(4*alp^2*bet - bet^3)*sig1^2 + (4*alp^2 -
    bet^2)*xi^2)*exp(2*alp*G[i,j]))*exp(bet*G[i,j]))*th_1 -
    (2*alp*bet*sig_0*xi^2 + 8*(2*alp^2*bet -
    alp*bet^2)*sig1^2*exp(alp*G[i,j]) - 2*(2*alp^2*bet -
    alp*bet^2)*sig_0^2 - 2*(2*alp^2*bet - alp*bet^2)*sig1^2 +
    (4*(2*alp^3*bet - alp^2*bet^2)*sig1^2*G[i,j] - 6*(2*alp^2*bet -
    alp*bet^2)*sig1^2 + (2*alp^2 -
    alp*bet)*xi^2)*exp(2*alp*G[i,j]))*exp(bet*G[i,j]))*exp(2*bet*G[i,j]*th_1)
    + (2*(2*alp^2*bet - alp*bet^2)*sig_0^2 - (2*alp^2 -
    alp*bet - 2*(2*alp^2 - alp*bet)*sig_0)*xi^2)*exp(2*alp*G[i,j]
    + bet*G[i,j]))*exp(-2*bet*G[i,j]*th_1)/(2*(2*alp^2*bet^2 -
    alp*bet^3)*th_1^2*exp(2*alp*G[i,j] + bet*G[i,j]) - (4*alp^3*bet -
    alp*bet^3)*th_1*exp(2*alp*G[i,j] + bet*G[i,j]) + (2*alp^3*bet -
    alp^2*bet^2)*exp(2*alp*G[i,j] + bet*G[i,j]))) *exp(alp*distG[i,j])
}


    cov_y_th[i,j] <-
    -(th_0*(exp(alp*G[i,j]) - 1) + y_0)*th_0*exp(-alp*G[i,j]) +
    ((alp*exp(alp*G[i,j]) - alp)*th_0^2 + alp*th_0*y_0 + sig1^2 +
    (alp*sig1^2*G[i,j] - sig1^2)*exp(alp*G[i,j]))*exp(-alp*G[i,j])/alp

    cov_thi_thj[i,j] <-
    sig1^2*G[i,j]
if(distG[i,j]<10^(-10)){b1_hat[i,j]<-0}else{
    b1_hat[i,j] <-
    -((th_0*(exp(alp*distG[i,j]/2) - 1) + y_0)*th_0*exp(-alp*distG[i,j]/2) -
    ((alp*exp(alp*distG[i,j]/2) - alp)*th_0^2 + alp*th_0*y_0 + sig1^2 +
    (alp*distG[i,j]/2*sig1^2 -
    sig1^2)*exp(alp*distG[i,j]/2))*exp(-alp*distG[i,j]/2)/alp)/(distG[i,j]/2*sig1^2)
}
    cov_ri_rj[i,j]<-
      cov_yi_yj[i,j]-2*b1_hat[i,j]*cov_y_th[i,j]+(b1_hat[i,j])^2*cov_thi_thj[i,j]
      }#end of j loop
  }#end of i loop
    return(cov_ri_rj)
}# end of CovRes function
#-------------------------------------------------------------------------------
NegLogLike<-function(x,regb=regb,n=n,G=G,distG=distG,one=one,response=response,pred.dataset=pred.dataset,dsm=dsm){
  badval<-(0.5)*.Machine$double.xmax
  alp<-x[1]
  bet<-x[2]
  th_1<-x[3]
   xi<-x[4]

  p_est<-c(alp,bet,th_1,xi)
  #print(p_est)
  V<-CovRes(p_est,n=n,G=G,distG=distG,one=one,response=response,pred.dataset=pred.dataset)
  #V<-as.matrix(nearPD(V)$mat)
  #CovRes(p_est,n=n,G=G,distG=distG,one=one,response=response,predictor=predictor)

  #regb[2]<-b1
  #put reg b here modified
  negloglike<-n/2*log(2*pi)+1/2*log(abs(det(V)))
  negloglike<-negloglike+1/2*t(response-dsm%*%regb)%*%pseudoinverse(V)%*%(response-dsm%*%regb)
  #print(negloglike)
  #print(c(alp,sig,negloglike))
  if(alp<0 || alp>100 || xi<0 || !is.finite(negloglike) || bet<0 || bet>100 || negloglike <= -1000){
      return(badval)
  }
  return(negloglike[1])
}#Neglog
#-------------------------------------------------------------------------------
#sim.dsm.oubmcir
sim.dsm.oubmcir<-function(model.params,pred.dataset=pred.dataset,x.tre=x.tre){
  alp<-model.params[1]
  bet<-model.params[2]
  th_1<-model.params[3]
  xi<-model.params[4]
  #b1<-x[3]

  n<-dim(pred.dataset)[1]
  one<-array(1,c(n,1))
  #print(one)
  second<-phyloCorrectFactor.OUBMCIR(c(alp,bet,th_1,xi),pred.dataset = pred.dataset,x.tre=x.tre)*pred.dataset
  return(second)
}
#-------------------------------------------------------------------------------
phyloCorrectFactor.OUBMCIR<-function(model.params, pred.params=pred.params, pred.dataset=pred.dataset,x.tre=x.tre){
  alp<-model.params[1]
  bet<-model.params[2]
  th_1<-model.params[3]
  xi<-model.params[4]

  rho_y_th<-0
  y_0<-0
  th_0<-0
  tree.vcv<-vcv(x.tre)
  sig1<-compute.optimal.sigma(model.params, pred.int.params = predictor.dataset.estimate$b_ini[1:dim(pred.dataset)[2]], pred.params.mtx = predictor.dataset.estimate$pred.params.mtx,pred.int.dataset=predictor.dataset.estimate$pred.int.dataset, G=G , opt.model = NULL)
  trh<-max(tree.vcv) #tree height
  cov_y_th <- -(th_0*(exp(alp*trh) - 1) + y_0)*th_0*exp(-alp*trh) +
    ((alp*exp(alp*trh) - alp)*th_0^2 + alp*th_0*y_0 + sig1^2 + (alp*sig1^2*trh - sig1^2)*exp(alp*trh))*exp(-alp*trh)/alp
  var_th <- sig1^2*trh
  return(c(cov_y_th/var_th))
}
#--------------------------------------r_sq&AICc-------------------------------------------
r_sq<-function(b_est=b_est,est_V=est_V,dsm=dsm,pred.dataset=pred.dataset,response=response){
  one<-array(1,c(length(response),1))
  inv_V_est<-pseudoinverse(est_V)
  mean_response<-(t(one)%*%inv_V_est%*%response)/sum(inv_V_est)
  SST<-t(response-c(mean_response))%*%inv_V_est%*%(response-c(mean_response))
  SSE<-t(response-dsm%*%b_est)%*%inv_V_est%*%(response-dsm%*%b_est)
  return((SST-SSE)/SST)
}

AICc<-function(n,k,negloglike){
  AIC<-2*k+negloglike
  AICc<-AIC+ 2*k*(k+1)/(n-k-1)
  return(AICc)
}
full.pred.dataset.estimate<-function(pred.dataset=pred.dataset,int.vars=int.vars, response=response){
  pred.names<-colnames(pred.dataset)
  for(j in 1 : dim(int.vars)[2]){
    int_j1 <- which(colnames(pred.dataset) == int.vars[j,][1])
    int_j2 <- which(colnames(pred.dataset) == int.vars[j,][2])
    pred.int_j12 <- matrix(pred.dataset[,int_j1]*pred.dataset[,int_j2],ncol=1)
    pred.int.dataset <- cbind(pred.int.dataset,pred.int_j12)
    temp.names <- paste(colnames(pred.dataset)[int_j1],colnames(pred.dataset)[int_j2],sep ="")
    pred.names = c(pred.names,temp.names)
  }
  colnames(pred.int.dataset)<-pred.names
  b_ini <- pseudoinverse(t(pred.int.dataset)%*%pred.int.dataset)%*%t(pred.int.dataset)%*%response
  rownames(b_ini) <- pred.names
  prdc <- sum(stri_length(colnames(pred.int.dataset)) == 1)
  pred.params.mtx <- array(0,c(prdc,prdc))
  colnames(pred.params.mtx)<-rownames(pred.params.mtx)<-pred.names[1:prdc]
  for(i in 1:(prdc-1)){
    for(j in i:prdc){
      if (i != j){
        cri.int<- rownames(b_ini) == paste(pred.names[i],pred.names[j],sep = "")
        if(any(cri.int)){
          pred.params.mtx[i,j]<-b_ini[which(cri.int),1]
          pred.params.mtx[j,i]<-pred.params.mtx[i,j]
        }
      }
    }
  }
  return( list(pred.int.dataset=pred.int.dataset, b_ini=b_ini, pred.params.mtx=pred.params.mtx))
}

#-------------------------------------------------------------------

compute.optimal.sigma <- function(model.params = model.params , response=response, pred.dataset=pred.dataset,G=G,opt.model=NULL){

    bet <- model.params[2]
    V<- (1/(2*bet))*(1-exp(-2*bet*max(G)))*exp(-2*bet*(max(G)-G))
    inv.V<-pseudoinverse(V)

  #prdc <- sum(stri_length(colnames(pred.dataset)) == 1)
  n<-dim(G)[1]
  one<-array(1,c(n,1))
  inv.G<-pseudoinverse(G)
  sig1<-0
  v.pred.array<-NULL
  b_est<- pseudoinverse(t(pred.dataset)%*%pred.dataset)%*%t(pred.dataset)%*%response

  for(bIndex in 1 : dim(pred.dataset)[2]){

      temp.m.pred <- c(t(one)%*%inv.V%*%pred.dataset[,bIndex])
      temp.m.pred <- c(temp.m.pred/(t(one)%*%inv.V%*%one))
      temp.v.pred <- t(pred.dataset[,bIndex] - temp.m.pred*one)%*%inv.V%*%(pred.dataset[,bIndex] - temp.m.pred*one)/n

    sig1<-sig1 + b_est[bIndex]*temp.v.pred
    }
    return(sig1)
}
#--------------------------------------------------------------------------------
prior<-function(model.params, pred.params){
  alp<-model.params[1]
  bet<-model.params[2]
  th_1<-model.params[3]
  xi<-model.params[4]
  #alpprior<- -dgamma(x=alp,shape=1,log=TRUE)
  #sigprior<- -dinvgamma(x=sig, shape=1,log=TRUE)
  #predprior<-dnorm(pred.params, mean=0,sd=0.001 ,log=TRUE)
  alpprior<- -dunif(x=alp,min=0,max=100,log=TRUE)
  betprior<- -dunif(x=bet,min=0,max=100,log=TRUE)
  th_1prior<- -dunif(x=th_1,min=0,max=100,log=TRUE)
  xiprior<- -dunif(x=xi,min=0,max=100,log=TRUE)
  predprior<- sum(-dunif(pred.params, min=-50,max=50 ,log=TRUE))
  return(alpprior + betprior + th_1prior + xiprior + predprior)
  }

posterior <- function(model.params=model.params, pred.params=pred.params, n=n, G=G, distG=distG, one=one,response=response,pred.dataset=pred.dataset,dsm=dsm){

#  print( prior(model.params,pred.params)  )
  return(NegLogLike(x=model.params,regb=b_ini, n=n, G=G,distG=distG,one=one,response=response,pred.dataset=pred.dataset,dsm=dsm)+ prior(model.params,pred.params))
  }


proposalfunction<-function(model.params, pred.params){
  alp<-model.params[1]
  bet<-model.params[2]
  th_1<-model.params[3]
  xi<-model.params[4]
  pps.alp<-rnorm(1,mean=alp,sd= alp/5)
  pps.bet<-rnorm(1,mean=bet,sd= bet/5)
  pps.th_1<-rnorm(1,mean=th_1,sd= th_1/5)
  pps.xi<-rnorm(1,mean=xi,sd=xi/5)
  pps.pred.params<-rnorm(length(pred.params),mean=pred.params,sd=abs(pred.params/5))
  return(c(abs(pps.alp), abs(pps.bet), abs(pps.th_1), abs(pps.xi), pps.pred.params))
  }


#iterations<-1000
run_metropolis_MCMC<-function(startvalue,iterations){
  chain=array(0,c(iterations+1,length(startvalue) ))
  chain[1,] <-startvalue
  for(i in 1:iterations){
    if(i%%500==0){print(i)}
    proposal<-proposalfunction(model.params=chain[i,1:4], pred.params=chain[i,5:length(startvalue)])
    probab<-exp(posterior(model.params=proposal[1:4], pred.params=proposal[5:length(proposal)], n=n, G=G, distG=distG, one=one,response=response,pred.dataset=pred.dataset,dsm=dsm)-posterior( model.params=chain[i,1:4], pred.params=chain[i,5:length(proposal)], n=n, G=G, distG=distG, one=one,response=response,pred.dataset=pred.dataset,dsm=dsm))
#probab
    if(runif(1)< probab){
      chain[i+1,]<-proposal
    }else{
      chain[i+1,]<-chain[i,]
    }
  #print(round(chain[i,],2))
  }
  return(mcmc(chain))
}
#-----------------------------test--------------------------------------
tree.size <- 7
num.pred <- 3
x.tre<-tree <- pbtree(n = tree.size,scale = 1)
#plot(tree)
#opt.model<-"OU"
#G<-vcv(x.tre)
response <- matrix(rnorm(tree.size),ncol = 1)
pred.dataset <- matrix(rnorm(tree.size*num.pred),ncol = num.pred)
opt.model.params<-0.2 #bet
# x <- c(5,6)
pred.int.dataset <- pred.dataset
pred.names<-colnames(pred.dataset) <- c("a","b","c")
int.vars <- matrix(c("a","b","c","c"),c(2,2))
predictor.dataset.estimate <- full.pred.dataset.estimate(pred.dataset=pred.dataset,int.vars=int.vars,response=response)
predictor.dataset.estimate

G<-vcv(x.tre)
distG<-2*(max(G)-G)
n<-dim(G)[1]
one<-array(1,c(n,1))
solveG<-solve(G)

dsm<-cbind(pred.dataset)#design matrix for initial estimate of b

b_ini<- solve(t(dsm)%*%dsm)%*%t(dsm)%*%response
regb <- b_ini

alpha0<-rexp(1,rate=5)
beta0<-rexp(1,rate=5)
th_10<-rexp(1,rate=50)
xi0<-rinvgamma(n=1,  shape= 1.1)

p0<-c(alpha0, beta0, th_10, xi0)
p0
#model.params <- p0
# c(rnorm(1),rexp(1))
# c(runif(1),runif(1))
print(compute.optimal.sigma(model.params = p0 , response=response, pred.dataset=pred.dataset,G=G,opt.model = NULL))
print(NegLogLike(x=p0,regb=b_ini,n=n,G=G,distG=distG,one=one,response=response,pred.dataset=pred.dataset,dsm=dsm))

startvalue<- c(alpha0,beta0,th_10,xi0,b_ini)
chain<-run_metropolis_MCMC(startvalue,1000)
#setwd("~/Dropbox/FCU/Teaching/Mentoring/2017Spring/TerryLai/R_code/maincode/")
#save.image("test_round2.RData")

#load("test_round2.RData")
burnIn<-100
acceptance <- 1- mean(duplicated(chain[-(1:burnIn),]))
print(acceptance)

par(mfrow=c(2,4))
hist(chain[-(1:burnIn),1],nclass=30, main="Posterior of alp", xlab= " True value = red line")
abline(v=mean(chain[-(1:burnIn),1]),col="red")

hist(chain[-(1:burnIn),2],nclass=30, main="Posterior of bet", xlab= " True value = red line")
abline(v=mean(chain[-(1:burnIn),2]),col="red")

hist(chain[-(1:burnIn),3],nclass=30, main="Posterior of th_1", xlab= " True value = red line")
abline(v=mean(chain[-(1:burnIn),3]),col="red")

hist(chain[-(1:burnIn),4],nclass=30, main="Posterior of xi", xlab= " True value = red line")
abline(v=mean(chain[-(1:burnIn),4]),col="red")


plot(chain[-(1:burnIn),1],type="l",xlab="True value of read line")
abline(h=mean(chain[-(1:burnIn),1]), col="red")

plot(chain[-(1:burnIn),2],type="l",xlab="True value of read line")
abline(h=mean(chain[-(1:burnIn),2]), col="red")

plot(chain[-(1:burnIn),3],type="l",xlab="True value of read line")
abline(h=mean(chain[-(1:burnIn),3]), col="red")

plot(chain[-(1:burnIn),4],type="l",xlab="True value of read line")
abline(h=mean(chain[-(1:burnIn),4]), col="red")
