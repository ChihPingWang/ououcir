rm(list=ls())
#install.packages("Sim.DiffProc")
library(Sim.DiffProc)
?st.int
a<-2
fexpr<-expression(a*exp(w-0.5*t))
res<-st.int(fexpr,type="ito",M=1,lower=0,upper=1)
ls(res)
median(res$X)

true.alpha.y<-2
ta<-2
ta
fexpr<-expression(true.alpha.y*exp(true.alpha.y*t)*w)
res<-st.int(fexpr,type="ito",M=1,lower=0,upper=1)
median(res$X)
?st.int

integrand<-function(x){1/((x+1)*sqrt(x))}
value<-integrate(integrand,lower=0,upper=Inf)
value$value

integrand<-function(s,alpha=alpha){
  alpha*exp(alpha*s)*rnorm(n=1,mean=0,sd=s)
  }
integrate(integrand,lower=0 ,upper=1, alpha=1)
