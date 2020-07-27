## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  error=TRUE,
  warning=FALSE
)

## ---- include=FALSE------------------------------------------------------
library(mlma)
#source('O:/My Documents/My Research/Research/Multilevel mediation analysis/mlma package/current version/R/mlma.r')

## ------------------------------------------------------------------------
set.seed(1)
n=20       # the number of observations in each group
J<-600/n   # there are 30 groups
level=rep(1:J,each=n)
alpha_211<-0.8     #covariates coefficients
alpha_1111<-0.8
alpha_2111<-0.8
beta_1<-0.4
beta_2<-0.4
beta_3<-0.4
beta_4<-0.4
beta_5<-0.4
v1=5              #the level 1 variance
v2=v1/5           #the level 2 variance

#The exposure variables
x1<-rbinom(600,1,0.5) #binary level 1 exposure, xij
x2<-rep(rnorm(J),each=n) #continuous level 2 exposure

#The mediators
m2<-rep(rbinom(J,1,exp(alpha_211*unique(x2^2))/(1+exp(alpha_211*unique(x2^2)))),each=n)    #level 2 binary mediator
u1<-rep(rnorm(J,0,0.5),each=n) #level 2 variance for mij
e1<-rnorm(n*J)  #level 1 variance for mij
m1<-u1+alpha_1111*x1+alpha_2111*x2+e1 #level 1 continuous mediator

#The response variable
u0<-rep(rnorm(J,0,v2),each=n)
e0<-rnorm(n*J,0,v1)
y<-u0+beta_1*x1+beta_2*x2+beta_3*ifelse(x2<=0,0,log(1+x2))+beta_4*m1+beta_5*m2+e0


## ------------------------------------------------------------------------
example1<-data.org(x=cbind(x1=x1,x2=x2), m=cbind(m1=m1,m2=m2), 
                   f01y=list(2,c("x","ifelse(x>0,log(x+1),0)")),
                   level=level,
                   f01km2=list(matrix(c(2,2),1,2),"x^2")) 

## ------------------------------------------------------------------------
mlma.e1<-mlma(y=y,data1=example1,intercept=F)
mlma.e1

## ------------------------------------------------------------------------
summary(mlma.e1)

## ------------------------------------------------------------------------
mlma.e1$f1   #the full model
mlma.e1$fm1  #models for level 1 mediators
mlma.e1$fm2  #models for level 2 mediators

## ------------------------------------------------------------------------
plot(mlma.e1)

## ------------------------------------------------------------------------
plot(mlma.e1,var="m1")
plot(mlma.e1,var="m2")

## ------------------------------------------------------------------------
boot.e1<-boot.mlma(y=y,data1=example1,echo=F,intercept = F)
summary(boot.e1)

## ------------------------------------------------------------------------
plot(boot.e1)
plot(boot.e1,var="m1")
plot(boot.e1,var="m2")

