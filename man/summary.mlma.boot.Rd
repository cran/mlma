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
\method{summary}{mlma.boot}(object, ..., alpha = 0.05)
\method{print}{summary.mlma.boot}(x, ...)
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
}
\author{
Qingzhao Yu (qyu@lsuhsc.edu), Bin Li (bli@lsu.edu).
}
\examples{
#1-1-1 model
data(sim.111)
temp111<-boot.mlma(y=sim.111$y, biny=FALSE, sim.111$x, levelx=1, 
                   m=sim.111$m, l1=1:2, c1=3,c1r=1,
                   f01y=c("x","log(x^2)"), f10y=c("x^2","sqrt(x+6)"), 
                   f20ky=list(2,c("x","x^3")), f01km1=list(2,"sqrt(x)+3"), 
                   f10km=list(2,"log(x+2)"), level=sim.111$level, boot=2)
summary(temp111)

#2-1-1 model
data(sim.211)
temp211<-boot.mlma(y=sim.211$y, biny=FALSE, x=sim.211$x, levelx=2, m=sim.211$m, 
                   l1=2,l2=1, c1=3,c1r=1, f01y=c("x","log(x^2)"), 
                   f02ky=list(1,c("x","x^2")), f20ky=list(2,c("x","x^3")),
                   f01km1=list(2,"sqrt(x)+3"), f01km2=list(1,c("x^1.2","x^2.3")), 
                   level=sim.211$level, boot=2)
summary(temp211)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{summary statistics}