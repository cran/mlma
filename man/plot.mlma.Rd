\name{plot.mlma}
\alias{plot.mlma}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot "mlma" Object
}
\description{
Plot the overall mediation effect or decomposed indirect effect from the selected mediator.
}
\usage{
\method{plot}{mlma}(x, ..., var = NULL, cate = FALSE, w2 = rep(1, length(object$de2)))
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
the name of the mediator that is to be plotted.  If var is NULL, plot the relative mediation effects of all mediators.
}
  \item{cate}{
TRUE when the mediator is a categorical variable.  By default is FALSE.
}
  \item{w2}{
the weight for observations at level 2, which should be the same order as unique(level[!is.na(level)]). The default is rep(1,length(object$de2)).
}
}
\details{
Plot the relative effects of direct effects and indirect effects of mediators at level 1 (if levelx=1) and level 2 respectively if var=NULL.  Otherwise, plot the indirect effect of var, the estimated differential effect of the predictor on var, and the predicted relationship between y and var at individual level and/or (aggregated) group level.
}
\author{
Qingzhao Yu (qyu@lsuhsc.edu), Bin Li (bli@lsu.edu).
}
\examples{
data(sim.111)
data1<-data.org(x=sim.111$x, levelx=1, m=sim.111$m, l1=1:2, c1=3, c1r=1, 
                f01y=c("x","log(x^2)"), f10y=c("x^2","sqrt(x+6)"), 
                f20ky=list(2,c("x","x^3")), f01km1=list(2,"sqrt(x)+3"), 
                f10km=list(2,"log(x+2)"), level=sim.111$level)
temp111<-mlma(y=sim.111$y, biny=FALSE, data1=data1, x=sim.111$x, levelx=1, m=sim.111$m, 
                l1=1:2,c1=3, c1r=1,level=sim.111$level)
plot(temp111)
plot(temp111,var="m.2")
plot(temp111,var="m.3")
plot(temp111,var="m.4", cate=TRUE)

data(sim.211)
data1<-data.org(x=sim.211$x, levelx=2, m=sim.211$m, l1=2,l2=1, c1=3, c1r=1,
                   f01y=c("x","log(x^2)"), f02ky=list(1,c("x","x^2")), 
                   f20ky=list(2,c("x","x^3")), f01km1=list(2,"sqrt(x)+3"), 
                   f01km2=list(1,c("x^1.2","x^2.3")), level=sim.211$level)
temp211<-mlma(y=sim.211$y, biny=FALSE, data1, x=sim.211$x, levelx=2, m=sim.211$m, 
           l1=2, l2=1, c1=3, c1r=1,level=sim.211$level)                               
plot(temp211)
plot(temp211,var="m.1")
plot(temp211,var="m.4", cate=TRUE)
plot(temp211,var="m.3")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{plot}
