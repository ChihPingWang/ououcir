####OUBM
rm(list=ls())
library(ape)
library(phytools)
#library(BMhyb)
library(corpcor)
library(stringi)
library(MCMCpack)
##given tree of n taxa,response,predictor.dataset,interaction,model
##output everything
#for(i in 1 : dim(int.vars)[2]){
#------------------------------------------------------------------------------

#oubm<-function(x.tre,response=response, pred.dataset=pred.dataset){

CovRes<-function(x,n=n,G=G,distG=distG,one=one,response=response,pred.dataset=pred.dataset){
    alp<-x[1]
    sig<-x[2]

    rho_y_th<-0
    y_0<-0
    th_0<-0


    #sig1<-compute.optimal.sigma(opt.model.params, pred.int.params = predictor.dataset.estimate$b_ini[1:dim(pred.dataset)[2]], pred.params.mtx = predictor.dataset.estimate$pred.params.mtx,pred.int.dataset=predictor.dataset.estimate$pred.int.dataset, G=G , opt.model = "BM")
    sig1<-compute.optimal.sigma(opt.model.params = 0 , response=response, pred.dataset=pred.dataset,G=G,opt.model="BM")


    cov_yi_yj<-array(0,c(n,n))
    cov_y_th<-array(0,c(n,n))
    cov_thi_thj<-array(0,c(n,n))
    b1_hat<-array(0,c(n,n))
    cov_ri_rj<-array(0,c(n,n))

    for(i in 1:n){
      for(j in 1:n){

        cov_yi_yj[i,j]<-
          sig1^2*G[i,j]*(exp(1/2*alp*distG[i,j]) - 1)^2*exp(-alp*distG[i,j]) - 2*((th_0*(exp(alp*G[i,j]) - 1) + y_0)*th_0*exp(-alp*G[i,j]) + (sig1*rho_y_th*sig - (alp*exp(alp*G[i,j]) - alp)*th_0^2 - alp*th_0*y_0 - sig1^2 - (alp*sig1^2*G[i,j] + sig1*rho_y_th*sig - sig1^2)*exp(alp*G[i,j]))*exp(-alp*G[i,j])/alp)*(exp(1/2*alp*distG[i,j]) - 1)*exp(-alp*distG[i,j]) - 1/2*(2*(th_0*(exp(alp*G[i,j]) - 1) + y_0)^2*exp(-2*alp*G[i,j]) - (2*sig1*rho_y_th*sig + 2*(alp*exp(2*alp*G[i,j]) - 2*alp*exp(alp*G[i,j]) + alp)*th_0^2 + 4*(alp*exp(alp*G[i,j]) - alp)*th_0*y_0 + 2*alp*y_0^2 - sig1^2 - sig^2 + (2*alp*sig1^2*G[i,j] + 2*sig1*rho_y_th*sig - 3*sig1^2 + sig^2)*exp(2*alp*G[i,j]) - 4*(sig1*rho_y_th*sig - sig1^2)*exp(alp*G[i,j]))*exp(-2*alp*G[i,j])/alp)*exp(-alp*distG[i,j])

        cov_y_th[i,j]<-
          -(th_0*(exp(alp*G[i,j]) - 1) + y_0)*th_0*exp(-alp*G[i,j]) - (sig1*rho_y_th*sig - (alp*exp(alp*G[i,j]) - alp)*th_0^2 - alp*th_0*y_0 - sig1^2 - (alp*sig1^2*G[i,j] + sig1*rho_y_th*sig - sig1^2)*exp(alp*G[i,j]))*exp(-alp*G[i,j])/alp

        cov_thi_thj[i,j] <- sig1^2*G[i,j]

        if(distG[i,j]<10^(-10)){b1_hat[i,j]<-0}else{
          b1_hat[i,j]<-   #1-   ((1-exp(-alp*distG[i,j]/2))/ (alp*distG[i,j]/2))
            -((th_0*(exp(alp*distG[i,j]/2) - 1) + y_0)*th_0*exp(-alp*distG[i,j]/2)
              - ((sig1*exp(alp*distG[i,j]/2) - sig1)*rho_y_th*sig +
                   (alp*exp(alp*distG[i,j]/2) - alp)*th_0^2 + alp*th_0*y_0 + sig1^2 +
                   (alp*sig1^2*distG[i,j]/2 - sig1^2)*exp(alp*distG[i,j]/2))*exp(-alp*distG[i,j]/2)/alp)/(sig1^2*distG[i,j]/2) # we got problem here
        }
        cov_ri_rj[i,j]<-
          cov_yi_yj[i,j]-2*b1_hat[i,j]*cov_y_th[i,j]+(b1_hat[i,j])^2*cov_thi_thj[i,j]

      }#end of j loop
    }#end of i loop

    return(cov_ri_rj)
  }#end of CovRes fcn

  #-------------------------------------------------------------------------------
  NegLogLike<-function(x,regb=regb,n=n,G=G,distG=distG,one=one,response=response,pred.dataset=pred.dataset,dsm=dsm){
    badval<-(0.5)*.Machine$double.xmax
    alp<-x[1]
    sig<-x[2]

    p_est<-c(alp,sig)
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
    if(alp<0 || alp>100 || sig<0 ||  !is.finite(negloglike) || negloglike <= -1000 ){
      return(badval)
    }


    return(negloglike[1])
  }#Neglog





#-------------------------------------------------------------------------------
#sim.dsm.oubm
sim.dsm.oubm<-function(model.params,pred.dataset=pred.dataset,x.tre=x.tre){
  alp<-model.params[1]
  sig<-model.params[2]
  #b1<-x[3]

  n<-dim(pred.dataset)[1]
  one<-array(1,c(n,1))
  #print(one)
  second<-phyloCorrectFactor.OUBM(c(alp,sig),pred.dataset = pred.dataset,x.tre=x.tre)*pred.dataset
  return(second)
}
#-------------------------------------------------------------------------------
phyloCorrectFactor.OUBM<-function(model.params, pred.params=pred.params, pred.dataset=pred.dataset,x.tre=x.tre){
  alp<-model.params[1]
  sig<-model.params[2]
  rho_y_th<-0
  y_0<-0
  th_0<-0
  tree.vcv<-vcv(x.tre)
  sig1<-compute.optimal.sigma(opt.model.params, pred.int.params = predictor.dataset.estimate$b_ini[1:dim(pred.dataset)[2]], pred.params.mtx = predictor.dataset.estimate$pred.params.mtx,pred.int.dataset=predictor.dataset.estimate$pred.int.dataset, G=G , opt.model = "BM")
  trh<-max(tree.vcv) #tree height
  cov_y_th<--(th_0*(exp(alp*trh) - 1) + y_0)*th_0*exp(-alp*trh) - (sig1*rho_y_th*sig - (alp*exp(alp*trh) - alp)*th_0^2 - alp*th_0*y_0 - sig1^2 - (alp*sig1^2*trh + sig1*rho_y_th*sig - sig1^2)*exp(alp*trh))*exp(-alp*trh)/alp
  var_th<-sig1^2*trh
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

compute.optimal.sigma <- function(opt.model.params = 0 , response=response, pred.dataset=pred.dataset,G=G,opt.model=opt.model){
  if(opt.model == "OU"){
    bet <- opt.model.params[1]
    V<- (1/(2*bet))*(1-exp(-2*bet*max(G)))*exp(-2*bet*(max(G)-G))
    inv.V<-pseudoinverse(V)
  }
  #prdc <- sum(stri_length(colnames(pred.dataset)) == 1)
  n<-dim(G)[1]
  one<-array(1,c(n,1))
  inv.G<-pseudoinverse(G)
  sig1<-0
  v.pred.array<-NULL
  b_est<- pseudoinverse(t(pred.dataset)%*%pred.dataset)%*%t(pred.dataset)%*%response

  for(bIndex in 1 : dim(pred.dataset)[2]){
    if(opt.model=="BM"){
      temp.m.pred <- c(t(one)%*%inv.G%*%pred.dataset[,bIndex])
      temp.m.pred <- c(temp.m.pred/(t(one)%*%inv.G%*%one))
      temp.v.pred <- t(pred.dataset[,bIndex] - temp.m.pred*one)%*%inv.G%*%(pred.dataset[,bIndex] - temp.m.pred*one)/n
    }
    if(opt.model=="OU"){
      temp.m.pred <- c(t(one)%*%inv.V%*%pred.dataset[,bIndex])
      temp.m.pred <- c(temp.m.pred/(t(one)%*%inv.V%*%one))
      temp.v.pred <- t(pred.dataset[,bIndex] - temp.m.pred*one)%*%inv.V%*%(pred.dataset[,bIndex] - temp.m.pred*one)/n
    }
    sig1<-sig1 + b_est[bIndex]*temp.v.pred
   }
    return(sig1)
  }


#model.params<-startvalue[1:2]
#pred.params<-startvalue[3:length(startvalue)]
prior<-function(model.params, pred.params){
  alp<-model.params[1]
  sig<-model.params[2]
  #alpprior<- -dgamma(x=alp,shape=1,log=TRUE)
  #sigprior<- -dinvgamma(x=sig, shape=1,log=TRUE)
  #predprior<-dnorm(pred.params, mean=0,sd=0.001 ,log=TRUE)
  alpprior<- -dunif(x=alp,min=0,max=100,log=TRUE)
  sigprior<- -dunif(x=sig,min=0, max=100,log=TRUE)
  predprior<- sum(-dunif(pred.params, min=-50,max=50 ,log=TRUE))
  return(alpprior + sigprior + predprior)
  }

posterior <- function(model.params=model.params, pred.params=pred.params, n=n, G=G, distG=distG, one=one,response=response,pred.dataset=pred.dataset,dsm=dsm){

#  print( prior(model.params,pred.params)  )
  return(NegLogLike(x=model.params,regb=b_ini, n=n, G=G,distG=distG,one=one,response=response,pred.dataset=pred.dataset,dsm=dsm)+ prior(model.params,pred.params))
  }


proposalfunction<-function(model.params, pred.params){
  alp<-model.params[1]
  sig<-model.params[2]
  pps.alp<-rnorm(1,mean=alp,sd= alp/5)
  pps.sig<-rnorm(1,mean=sig,sd= sig/5)
  pps.pred.params<-rnorm(length(pred.params),mean=pred.params,sd=abs(pred.params/5))
  return(c(abs(pps.alp), abs(pps.sig), pps.pred.params))
  }


#iterations<-1000
run_metropolis_MCMC<-function(startvalue,iterations){
  chain=array(0,c(iterations+1,length(startvalue) ))
  chain[1,] <-startvalue
  for(i in 1:iterations){
    if(i%%500==0){print(i)}
    proposal<-proposalfunction(model.params=chain[i,1:2], pred.params=chain[i,3:length(startvalue)])
    probab<-exp(posterior(model.params=proposal[1:2], pred.params=proposal[3:length(proposal)], n=n, G=G, distG=distG, one=one,response=response,pred.dataset=pred.dataset,dsm=dsm)-posterior( model.params=chain[i,1:2], pred.params=chain[i,3:length(proposal)], n=n, G=G, distG=distG, one=one,response=response,pred.dataset=pred.dataset,dsm=dsm))
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
opt.model<-"BM"
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
sig0<-rinvgamma(n=1,  shape= 1.1)

p0<-c(alpha0, sig0)
p0
# c(rnorm(1),rexp(1))
# c(runif(1),runif(1))
print(compute.optimal.sigma(opt.model.params = 0 , response=response, pred.dataset=pred.dataset,G=G,opt.model="BM"))
print(NegLogLike(x=p0,regb=b_ini,n=n,G=G,distG=distG,one=one,response=response,pred.dataset=pred.dataset,dsm=dsm))

startvalue<- c(alpha0,sig0,b_ini)
chain<-run_metropolis_MCMC(startvalue,100000)
#setwd("~/Dropbox/FCU/Teaching/Mentoring/2017Spring/TerryLai/R_code/maincode/")
setwd("~/Dropbox/TerryLai/R_code/maincode")
#save.image("test_round2.RData")
system("pwd")
#load("test_round2.RData")
burnIn<-5000
acceptance <- 1- mean(duplicated(chain[-(1:burnIn),]))
print(acceptance)

par(mfrow=c(2,2))
hist(chain[-(1:burnIn),1],nclass=30, main="Posterior of alp", xlab= " True value = red line")
abline(v=mean(chain[-(1:burnIn),1]))

hist(chain[-(1:burnIn),2],nclass=30, main="Posterior of sig", xlab= " True value = red line")
abline(v=mean(chain[-(1:burnIn),2]))

plot(chain[-(1:burnIn),1],type="l",xlab="True value of read line")
abline(h=mean(chain[-(1:burnIn),1]), col="red")

plot(chain[-(1:burnIn),2],type="l",xlab="True value of read line")
abline(h=mean(chain[-(1:burnIn),2]), col="red")





#oubm(x.tre,response=response, pred.dataset=pred.dataset)
