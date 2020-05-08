\name{print.mlma.boot}
\alias{print.mlma.boot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Print "mlma.boot" Object
}
\description{
print the level 1 and level 2 mediation effecs from the object.
}
\usage{
\method{print}{mlma.boot}(x,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
an "mlma.boot" object.
}
  \item{...}{other arguments}
}
\author{
Qingzhao Yu (qyu@lsuhsc.edu), Bin Li (bli@lsu.edu).
}
\examples{
data(sim.111)
data2<-data.org(ifelse(sim.111$x>1.9,1,0), m=sim.111$m, 
                f10y=list(1,c("x^2","sqrt(x+6)")), 
                f20ky=list(2,c("x","x^3")), 
                f10km=list(matrix(c(2,1),1),"log(x+2)"), level=sim.111$level)
temp2.boot<-boot.mlma(y=sim.111$y, data1=data2,boot=2)

print(temp2.boot)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ print }
