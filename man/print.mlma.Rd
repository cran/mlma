\name{print.mlma}
\alias{print.mlma}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Print "mlma" Object
}
\description{
print the level 1 and level 2 mediation effecs from the object.
}
\usage{
\method{print}{mlma}(x, ..., w2 = rep(1, length(object$de2)))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
an "mlma" object.
}
  \item{\dots}{
arguments to be passed to methods.
}
  \item{w2}{
the weight for observations at level 2, which should be the same order as unique(level[!is.na(level)]). The default is rep(1,length(object$de2)).
}
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
print(temp111)

data(sim.211)
data1<-data.org(x=sim.211$x, levelx=2, m=sim.211$m, l1=2,l2=1, c1=3, c1r=1,
                   f01y=c("x","log(x^2)"), f02ky=list(1,c("x","x^2")), 
                   f20ky=list(2,c("x","x^3")), f01km1=list(2,"sqrt(x)+3"), 
                   f01km2=list(1,c("x^1.2","x^2.3")), level=sim.211$level)
temp211<-mlma(y=sim.211$y, biny=FALSE, data1, x=sim.211$x, levelx=2, m=sim.211$m, 
           l1=2, l2=1, c1=3, c1r=1,level=sim.211$level)                               
print(temp211)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ print }
