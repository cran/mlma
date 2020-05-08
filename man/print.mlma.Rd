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
\method{print}{mlma}(x, ..., w2 = rep(1, length(object$de2)),digits=2)
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
  \item{digits}{
the number of digits to print.
}
}
\author{
Qingzhao Yu (qyu@lsuhsc.edu), Bin Li (bli@lsu.edu).
}
\examples{
data(sim.111)
#with a level 1 exposure
data(sim.111)
data2<-data.org(sim.111$x, m=sim.111$m, 
                   f10y=list(1,c("x^2","sqrt(x+6)")), 
                   f20ky=list(2,c("x","x^3")), 
                   f10km=list(matrix(c(2,1),1),"log(x+2)"), level=sim.111$level)
temp2<-mlma(y=sim.111$y, data1=data2)
print(temp2)

#with a level 2 exposure
data(sim.211)
data1<-data.org(x=sim.211$x, m=sim.211$m, 
                   f01y=list(1,c("x","log(x^2)")), f02ky=list(1,c("x","x^2")),  
                   f20ky=list(2,c("x","x^3")), f01km2=list(matrix(c(1,1),1),c("x^1.2","x^2.3")),
                   f01km1=list(matrix(c(2,1),1),"sqrt(x)+3"),level=sim.211$level) 
temp1<-mlma(y=sim.211$y, data1) 
temp1

#with both level 1 and 2 exposure
data3<-data.org(x=cbind(sim.211$x,sim.111$x), m=sim.211$m, 
                   f01y=list(1,c("x","log(x^2)")), f02ky=list(1,c("x","x^2")), 
                   f20ky=list(2,c("x","x^3")), f01km1=list(matrix(c(2,1),1),"sqrt(x)+3"), 
                   f01km2=list(matrix(c(1,1),1),c("x^1.2","x^2.3")), level=sim.211$level)
temp3<-mlma(y=sim.211$y, data3) 
temp3
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ print }
