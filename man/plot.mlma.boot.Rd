\name{plot.mlma.boot}
\alias{plot.mlma.boot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot the "mlma.boot" Object
}
\description{
For the mediator identified by var, the function draws the level 1 and/or (aggregated) level 2 indirect effects versus the predictor and the confidence bands at alpha significance level.
}
\usage{
\method{plot}{mlma.boot}(x, ..., var = NULL, alpha = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
an "mlma" object.
}
  \item{\dots}{
arguments to be passed to methods.
}
  \item{var}{
the name of the mediator that is to be plotted.  
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
plot(temp111,var="m.2")
plot(temp111,var="m.3")
plot(temp111,var="m.4")
 

#2-1-1 model
data(sim.211)
temp211<-boot.mlma(y=sim.211$y, biny=FALSE, x=sim.211$x, levelx=2, m=sim.211$m, 
                   l1=2,l2=1, c1=3,c1r=1, f01y=c("x","log(x^2)"), 
                   f02ky=list(1,c("x","x^2")), f20ky=list(2,c("x","x^3")),
                   f01km1=list(2,"sqrt(x)+3"), f01km2=list(1,c("x^1.2","x^2.3")), 
                   level=sim.211$level, boot=2)
plot(temp211,var="m.1")
plot(temp211,var="m.4")
plot(temp211,var="m.3")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{plot}
