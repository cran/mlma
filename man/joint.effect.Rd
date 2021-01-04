\name{joint.effect}
\alias{joint.effect}
\alias{print.joint.effect}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Find the joint effect of a list of mediators.
}
\description{
This function provides the inference information (estimation, standard deviation, etc.) for the joint and relative effect of a list of mediators.
}
\usage{
joint.effect(object,var.list,digits=4,...,alpha=0.05)
\method{print}{joint.effect}(x,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
an "boot.mlma" object.
}
  \item{x}{
an "joint.effect" object from the joint.effect function.
}
 \item{var.list}{
a list of mediators whose joint effects are of interest.
}
  \item{\dots}{
arguments to be passed to methods.
}
  \item{alpha}{
the significance level at which to draw the confidence bands.  
}
  \item{digits}{
number of digits to be printed.  
}
}
\author{
Qingzhao Yu (qyu@lsuhsc.edu), Bin Li (bli@lsu.edu).
}
\examples{
data(sim.111)
data(sim.211)
data3<-data.org(x=cbind(ifelse(sim.211$x>1.9,1,0),sim.111$x), m=sim.211$m, 
                f20ky=list(2,c("x","x^3")), f01km1=list(matrix(c(2,1),1),"sqrt(x)+3"), 
                f01km2=list(matrix(c(1,1),1),c("x^1.2","x^2.3")), level=sim.211$level)
temp3.boot<-boot.mlma(y=sim.211$y, data1=data3,boot=2) 
joint.effect(temp3.boot, var.list=c("m.4"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{joint effect statistics}
