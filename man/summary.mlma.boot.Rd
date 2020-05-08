\name{summary.mlma.boot}
\alias{summary.mlma.boot}
\alias{print.summary.mlma.boot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summary of "mlma.boot" Object
}
\description{
This function provide summary statistics for all mediation effects.
}
\usage{
\method{summary}{mlma.boot}(object, ..., alpha = 0.05,RE=FALSE,digits=4)
\method{print}{summary.mlma.boot}(x, ..., digits=x$digits)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
an "mlma" object.
}
 \item{x}{
a summary.mlma.boot object created initially call to summary.mlma.boot.
}
  \item{\dots}{
arguments to be passed to methods.
}
  \item{alpha}{
the significance level at which to draw the confidence bands.  
}
  \item{RE}{
if true, print the relative effects, otherwise show the mediation effects.  
}
  \item{digits}{
the number of digits to print.
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
summary(temp3.boot, digits=2)
summary(temp3.boot, RE=TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{summary statistics}
