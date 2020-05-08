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
\method{plot}{mlma}(x, ..., var = NULL, cate = FALSE, w2 = rep(1,
                   length(object$de2)))
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
an indicator of whether the variable is categorical or not.
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
data(sim.211)
data1<-data.org(x=sim.211$x, m=sim.211$m, 
                   f01y=list(1,c("x","log(x^2)")), f02ky=list(1,c("x","x^2")),  
                   f20ky=list(2,c("x","x^3")), f01km2=list(matrix(c(1,1),1),c("x^1.2","x^2.3")),
                   f01km1=list(matrix(c(2,1),1),"sqrt(x)+3"),level=sim.211$level) 
temp1<-mlma(y=sim.211$y, data1) 
plot(temp1)
plot(temp1,var="m.1")
plot(temp1,var="m.3")
#plot(temp1,var="m.4")

data(sim.111)
data2<-data.org(sim.111$x, m=sim.111$m, 
                   f10y=list(1,c("x^2","sqrt(x+6)")), 
                   f20ky=list(2,c("x","x^3")), 
                   f10km=list(matrix(c(2,1),1),"log(x+2)"), level=sim.111$level)
temp2<-mlma(y=sim.111$y, data1=data2)
plot(temp2)
plot(temp2,var="m.2")
#plot(temp2,var="m.4")
plot(temp2,var="m.3")

data3<-data.org(x=cbind(sim.211$x,sim.111$x), m=sim.211$m, 
                   f01y=list(1,c("x","log(x^2)")), f02ky=list(1,c("x","x^2")), 
                   f20ky=list(2,c("x","x^3")), f01km1=list(matrix(c(2,1),1),"sqrt(x)+3"), 
                   f01km2=list(matrix(c(1,1),1),c("x^1.2","x^2.3")), level=sim.211$level)
temp3<-mlma(y=sim.211$y, data3) 
plot(temp3)
plot(temp3,var="m.1")
#plot(temp3,var="m.4")
plot(temp3,var="m.3")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{plot}
