\name{plot.mlma.boot}
\alias{plot.mlma.boot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot the "mlma.boot" Object
}
\description{
For the mediator identified by var, the function draws the level 1 and/or (aggregated) level 2 indirect effects versus the predictor and the confidence bands at alpha significance level. If var is NULL, draw the relative mediation effects with confidence intervals.
}
\usage{
\method{plot}{mlma.boot}(x, ..., var = NULL, alpha = 0.05, quant=FALSE)
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
  \item{quant}{
if true, confidence interval is calculated using quantil when plot the relative effects. By default, the CIs are calculated using normal approximation.  This argument does nothing for the CIs calculated when var is not null.
}
}
\author{
Qingzhao Yu (qyu@lsuhsc.edu), Bin Li (bli@lsu.edu).
}
\examples{
data(sim.211)
data1<-data.org(x=ifelse(sim.211$x>2.2,1,0), m=sim.211$m, 
                f02ky=list(1,c("x","x^2")),  
                f20ky=list(2,c("x","x^3")), f01km2=list(matrix(c(1,1),1),c("x^1.2","x^2.3")),
                f01km1=list(matrix(c(2,1),1),"sqrt(x)+3"),level=sim.211$level) 
temp1.boot<-boot.mlma(y=sim.211$y, data1=data1,boot=2)  
plot(temp1.boot)
plot(temp1.boot,var="m.1")
plot(temp1.boot,var="m.3")
plot(temp1.boot,var="m.4")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{plot}
