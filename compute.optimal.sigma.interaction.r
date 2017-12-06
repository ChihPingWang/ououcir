rm(list=ls())
library(phangorn)
library(corpcor)
library(stringi)


compute.optimal.sigma<-function(model.params=NULL,pred.int.params = pred.int.params,pred.int.dataset=pred.int.dataset, V=V , opt.model = opt.model){
  
  #given dataset, identify which preditors, interactions 
  #for example r = 5 predictors, s<- 2 interactions
  pred.int.params <- b_ini
  prdc <- sum(stri_length(colnames(pred.int.dataset)) == 1)
  colnames(pred.int.dataset)[1:prdc]
  int.prdc <- dim(pred.int.dataset)[2] - prdc
  # itrc <- 
  #dim(predictor.dataset)
  
  n<-dim(V)[1]
  one<-array(1,c(n,1))
  inv.V<-pseudoinverse(V)
  sig1<-0

    v.pred.array<-NULL

  for(bIndex in 1:prdc){
    temp.m.pred<-c(t(one)%*%inv.V%*%pred.int.dataset[,bIndex])
    temp.m.pred<-c(temp.m.pred/(t(one)%*%inv.V%*%one))
    temp.v.pred <- t(pred.int.dataset[,bIndex] - temp.m.pred*one)%*%inv.V%*%(pred.int.dataset[,bIndex] - temp.m.pred*one)/n
    v.pred.array<-c(v.pred.array,temp.v.pred)
   
    for(jIndex in 1:prdc){
      temp.bij.EX.sq <- 0
      temp.bij.EX <- 0
      if(jIndex != bIndex){
        if(opt.model == "BM"){
          EX.pred <-0 
          EX.sq.pred <- max(V)
          temp.bij.EX.sq <-  pred.params.mtx[bIndex,jIndex]^2*temp.v.pred*EX.sq.pred
          
        }
        if(opt.model == "OU"){
          x_0 <- temp.m.pred
          EX.pred <- x_0*exp(bet*max(V))
          EX.sq.pred <- (1/(2*bet))*(1 - exp(-2*bet*max(V))) + x_0^2*exp(-2*bet*max(V))# + mu*(1-exp(-bet*max(V))))^2
        
        
        
        temp.bij.EX.sq <-  pred.params.mtx[bIndex,jIndex]^2*temp.v.pred*EX.sq.pred
        temp.bij.EX <- pred.params.mtx[bIndex,jIndex]*EX.pred
        
        }#end of ou
      }
    }# jIndex
    sig1<-sig1+pred.int.params[bIndex]^2*temp.v.pred + temp.bij.EX.sq*temp.v.pred + 2*pred.int.params[bIndex]*temp.v.pred*temp.bij.EX
    }# THIS LOOP IS FOR OPTIMAL ASSUMING BM
  
    
    ##WILL DO OPTIMAL ASSUMING OU
    
    
    
    
  
    #     for(bIndex in (prdc+1):(prdc+itrc)){
  #      
  #   pred.params[bIndex]^2*temp.v.pred*max(V)
  #   
  #       
  #       
  #       
  # }
  # 
  # 
#  START FROM HERE: COINSTRCUT ARRAY FIRST, THEN IDENTIFY THE POSITION OF INTERACTION IN THE FORMULA bij....multreg.pdf
 # colnames(predictor.dataset)<-c("X1","X2","X3","X12","X23")
  
  
  
  sig1<-sqrt(sig1)
  return(sig1)
}

size<-7
num.pred<-3
pred.params<-rexp(3, rate=5)
pred.params

predictor.dataset<-matrix(rnorm(size*num.pred),ncol=num.pred)
colnames(predictor.dataset) <- c("a","b","c")
int.vars <- expand.grid("a",c("b","c"))

predictor.dataset<-data.frame(predictor.dataset)
predictor.dataset$a
int_datset_ac <- predictor.dataset$a * predictor.dataset$c
pred.names <- colnames(predictor.dataset)
pred.int.dataset <- predictor.dataset
for(j in 1 : dim(int.vars)[2]){
  #typeof(int.vars)
  #j <- 1
  #int.vars = c("a","c")
  int_j1 <- which(colnames(predictor.dataset) == unlist(int.vars[j,])[1] )
  int_j2 <- which(colnames(predictor.dataset) == unlist(int.vars[j,])[2] )
  colnames(predictor.dataset)[3]
  pred.int_j12 <- matrix(predictor.dataset[,int_j1]*predictor.dataset[,int_j2],ncol = 1)
  pred.int.dataset <- cbind(pred.int.dataset,pred.int_j12)
  temp.names <- paste(colnames(predictor.dataset)[int_j1],colnames(predictor.dataset)[int_j2],sep ="")
  pred.names = c(pred.names,temp.names)
}  
colnames(pred.int.dataset) <- pred.names
typeof(pred.names)
pred.int.dataset <- as.matrix(pred.int.dataset)
response <- matrix(rnorm(size),ncol = 1)
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

V<-diag(abs(rnorm(size)),c(size,size))
print(compute.optimal.sigma(pred.params, predictor.dataset = predictor.dataset, V=V))


