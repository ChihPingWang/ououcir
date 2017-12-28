rm(list=ls())
library(sde)
#?rcCIR

#rcCIR(n=10,Dt=0.1,x0=1,theta=c(6,2,2)

sim.cir<-function(x){
    mu <-x[1]
    alp<-x[2]
    sig<-x[3]
    cir0<-rcCIR(n=n/2,Dt=0.01,x0=1,theta=c(mu,alp,sig))
    cir1<-rcCIR(n=n/2-1,Dt=0.01,x0=20,theta=c(mu,alp,sig))
    cir2<-rcCIR(n=n/2-1,Dt=0.01,x0=100,theta=c(mu,alp,sig))
    result<-list(cir0=cir0,cir1=cir1,cir2=cir2)
    return(result)
    }
      
sim.bm<-function(x){
    sig<-x[1]
    bm00<-cumsum(c(0,sig*rnorm(n)))
    bm0<-cumsum(c(0,sig*rnorm(n/2)))
    bm1<-cumsum(c(bm0[n/2],sig*rnorm(n/2-1)))
    bm2<-cumsum(c(bm0[n/2],sig*rnorm(n/2-1)))
    result<-list(bm0=bm0,bm1=bm1,bm2=bm2)
    return(result)
    }
      
      n<-200
      #out<-sim.bm(1)
      mu<-100;alp<-10;sig<-1
      print(2*mu*alp > sig^2)
      
      out<-sim.cir(c(mu,alp,sig))
      #miny<-min(out$bm0,out$bm1,out$bm2)
      #maxy<-max(out$bm0,out$bm1,out$bm2)
      miny<-min(out$cir0,out$cir1,out$cir2)
      maxy<-max(out$cir0,out$cir1,out$cir2)
      #plot(out$bm0,xlim=c(-n/10,n+n/10),ylim=c(min(out$bm00,out$bm0,out$bm1,out$bm2),1.4*max(out$bm00,out$bm0,out$bm1,out$bm2)),type="l",lwd=2,ylab="",xlab="",axes=F,frame.plot=F, xaxt='n', ann=FALSE, yaxt='n')
      
      # #load("illtree1.RData")
      # 
      # 
      # #setEPS()
      # #postscript("tree.eps",height=8,width=13)
      # par(mfrow=c(1,2))
      # plot(out$bm0,xlim=c(-n/10,n+n/10),ylim=c(1.5*min(out$bm0,out$bm1,out$bm2),1.2*max(out$bm0,out$bm1,out$bm2)),type="l",lwd=4,ylab="",xlab="",frame.plot=F,  xaxt='n', cex=2,yaxt='n')
      # axis(side=2,at=c(miny, (miny+maxy)/2,maxy) ,labels=c("",expression(y[t]),""), cex.axis=2)
      # segments(c(5001:9999),out$bm1[1:4999], c( 5002:10000),out$bm1[2:5000],col="red",lwd=4)
      # segments(c(5001:7499),out$bm2[1:2499], c( 5002:7500),out$bm2[2:2500],col="blue",lwd=4)
      # segments(1,1.25*miny,n,1.25*miny,col="blue",lwd=3,lty=1)
      # segments(n,out$bm1[n/2],  n  ,1.3*miny,col="blue",lwd=3,lty=2)
      # segments(3/4*n,out$bm2[2500], 3/4*n ,1.25*miny,col="blue",lwd=3,lty=2)
      # segments(1,0,1,1.25*miny,col="blue",lwd=3,lty=2)
      # segments(n/2,out$bm0[n/2],n/2,1.25*miny,col="blue",lwd=3,lty=2)
      # 
      # text(-n/20,0,expression(y[0]),cex=2)
      # text(1.05*n,out$bm1[n/2],expression(y[i]),cex=2)
      # text(3*n/4*1.05,out$bm2[n/4],expression(y[j]),cex=2)
      # text(n/2, 1*out$bm0[n/2] + 0.8* abs( out$bm0[n/2]),expression(y[a]),cex=2)
      # 
      # text(5,1.45*miny,expression(paste("t=",0,sep="")),cex=2)
      # text(n/2,1.45*miny,expression(t[a]),cex=2)
      # text(3/4*n,1.45*miny,expression(paste(t[a],"+", t[j],sep="")),cex=2)
      # text(1*n,1.45*miny,expression(paste(t[a],"+", t[i],sep="")),cex=2)
      # 
      # #The tree
      # 
      # plot(NA,xlim=c(-n/10,n+n/10),ylim=c(1.5*min(out$bm0,out$bm1,out$bm2),1.2*max(out$bm0,out$bm1,out$bm2)),type="l",lwd=4,ylab="",xlab="",frame.plot=F,  xaxt='n', cex=2,yaxt='n')
      # segments(c(1:4999),(max(out$bm0,out$bm1,out$bm2)+min(out$bm0,out$bm1,out$bm2))/2, c( 2:5000),(max(out$bm0,out$bm1,out$bm2)+min(out$bm0,out$bm1,out$bm2))/2,lwd=15)
      # segments(c(5001:9999), max(out$bm0,out$bm1,out$bm2) , c( 5002:10000),max(out$bm0,out$bm1,out$bm2),lwd=15)
      # segments(c(5001:7499),min(out$bm0,out$bm1,out$bm2), c( 5002:7500),min(out$bm0,out$bm1,out$bm2),lwd=15)
      # segments(5001,min(out$bm0,out$bm1,out$bm2), 5001,max(out$bm0,out$bm1,out$bm2),lwd=15)
      # segments(1,1.25*miny,n,1.25*miny,col="blue",lwd=3,lty=1)
      # segments(n,1.25*miny,  n  ,1.30*miny,col="blue",lwd=3,lty=2)
      # text(1.05*n,max(out$bm0,out$bm1,out$bm2), expression(i),cex=3)
      # text(3*n/4*1.05,min(out$bm0,out$bm1,out$bm2),expression(j),cex=3)
      # segments(1,1.25*miny,  1  ,1.30*miny,col="blue",lwd=3,lty=2)
      # text(5,1.45*miny,expression(paste("t=",0,sep="")),cex=2)
      # segments(n/2,1.25*miny, n/2  ,1.30*miny,col="blue",lwd=3,lty=2)
      # text(n/2,1.45*miny,expression(t[a]),cex=2)
      # segments(3*n/4,1.25*miny, 3*n/4  ,1.30*miny,col="blue",lwd=3,lty=2)
      # text(3/4*n,1.45*miny,expression(paste(t[a],"+", t[j],sep="")),cex=2)
      # text(1*n,1.45*miny,expression(paste(t[a],"+", t[i],sep="")),cex=2)
      # 
      # #dev.off()
      
      
      par(mfrow=c(1,2))
      plot(out$cir0,xlim=c(-n/10,n+n/10),ylim=c(1.5*min(out$cir0,out$cir1,out$cir2),1.2*max(out$cir0,out$cir1,out$cir2)),type="l",lwd=2,ylab="",xlab="",frame.plot=F,  xaxt='n', cex=2,yaxt='n')
      axis(side=2,at=c(miny, (miny+maxy)/2,maxy) ,labels=c("",expression(y[t]),""), cex.axis=2)
      #segments(c(5001:9999),out$cir1[1:4999], c( 5002:10000),out$cir1[2:5000],col="red",lwd=4)
      #segments(c(5001:7499),out$cir2[1:2499], c( 5002:7500),out$cir2[2:2500],col="blue",lwd=4)
     
      segments(c( (n/2+1):(n-1)),out$cir1[1:(n/2-1)], c( (n/2+2):n),out$cir1[2:(n/2)],col="red",lwd=2)
      segments(c( (n/2+1):(3*n/4-1)),out$cir2[1:(n/4-1)], c( (n/2+2):(3*n/4)),out$cir2[2:(n/4)],col="blue",lwd=2)
      
      
      segments(1,1.25*miny,n,1.25*miny,col="blue",lwd=3,lty=1)
      segments(n,out$cir1[n/2],  n  ,1.3*miny,col="blue",lwd=3,lty=2)
      segments(3/4*n,out$cir2[n/4], 3/4*n ,1.25*miny,col="blue",lwd=3,lty=2)
      segments(1,0,1,1.25*miny,col="blue",lwd=3,lty=2)
      segments(n/2,out$cir0[n/2],n/2,1.25*miny,col="blue",lwd=3,lty=2)
      
      text(-n/20,0,expression(y[0]),cex=2)
      text(1.05*n,out$cir1[n/2],expression(y[i]),cex=2)
      text(3*n/4*1.05,out$cir2[n/4],expression(y[j]),cex=2)
      text(n/2, 1*out$cir0[n/2] + 0.8* abs( out$cir0[n/2]),expression(y[a]),cex=2)
      
      text(5,1.45*miny,expression(paste("t=",0,sep="")),cex=2)
      text(n/2,1.45*miny,expression(t[a]),cex=2)
      text(3/4*n,1.45*miny,expression(paste(t[a],"+", t[j],sep="")),cex=2)
      text(1*n,1.45*miny,expression(paste(t[a],"+", t[i],sep="")),cex=2)
      
      #The tree
      
      plot(NA,xlim=c(-n/10,n+n/10),ylim=c(1.5*min(out$cir0,out$cir1,out$cir2),1.2*max(out$cir0,out$cir1,out$cir2)),type="l",lwd=4,ylab="",xlab="",frame.plot=F,  xaxt='n', cex=2,yaxt='n')
      #segments(c(1:4999),(max(out$cir0,out$cir1,out$cir2)+min(out$cir0,out$cir1,out$cir2))/2, c( 2:5000),(max(out$cir0,out$cir1,out$cir2)+min(out$cir0,out$cir1,out$cir2))/2,lwd=15)
      #segments(c(5001:9999), max(out$cir0,out$cir1,out$cir2) , c( 5002:10000),max(out$cir0,out$cir1,out$cir2),lwd=15)
      #segments(c(5001:7499),min(out$cir0,out$cir1,out$cir2), c( 5002:7500),min(out$cir0,out$cir1,out$cir2),lwd=15)
      #segments(5001,min(out$cir0,out$cir1,out$cir2), 5001,max(out$cir0,out$cir1,out$cir2),lwd=15)
      
      segments(c(1:(n/2-1)),(max(out$cir0,out$cir1,out$cir2)+min(out$cir0,out$cir1,out$cir2))/2, c( 2:n/2),(max(out$cir0,out$cir1,out$cir2)+min(out$cir0,out$cir1,out$cir2))/2,lwd=15)
      segments(c((n/2+1):(n-1)), max(out$cir0,out$cir1,out$cir2) , c( (n/2+2):n),max(out$cir0,out$cir1,out$cir2),lwd=15)
      segments(c((n/2+1):(3*n/4-1)),min(out$cir0,out$cir1,out$cir2), c( (n/2+2):(3*n/4)),min(out$cir0,out$cir1,out$cir2),lwd=15)
      segments((n/2+1),min(out$cir0,out$cir1,out$cir2), (n/2+1),max(out$cir0,out$cir1,out$cir2),lwd=15)
      
      
      segments(1,1.25*miny,n,1.25*miny,col="blue",lwd=3,lty=1)
      segments(n,1.25*miny,  n  ,1.30*miny,col="blue",lwd=3,lty=2)
      text(1.05*n,max(out$cir0,out$cir1,out$cir2), expression(i),cex=3)
      text(3*n/4*1.05,min(out$cir0,out$cir1,out$cir2),expression(j),cex=3)
      segments(1,1.25*miny,  1  ,1.30*miny,col="blue",lwd=3,lty=2)
      text(5,1.45*miny,expression(paste("t=",0,sep="")),cex=2)
      segments(n/2,1.25*miny, n/2  ,1.30*miny,col="blue",lwd=3,lty=2)
      text(n/2,1.45*miny,expression(t[a]),cex=2)
      segments(3*n/4,1.25*miny, 3*n/4  ,1.30*miny,col="blue",lwd=3,lty=2)
      text(3/4*n,1.45*miny,expression(paste(t[a],"+", t[j],sep="")),cex=2)
      text(1*n,1.45*miny,expression(paste(t[a],"+", t[i],sep="")),cex=2)
      
      