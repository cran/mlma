\name{summary.mlma}
\alias{summary.mlma}
\alias{print.summary.mlma}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summary of "mlma" Object
}
\description{
This function provides ANOVA tests on the predictors and mediators in the full model and on the predictors for models in explaining each mediators.
}
\usage{
\method{summary}{mlma}(object,...,type="III")
\method{print}{summary.mlma}(x, ...)
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
  \item{type}{
type of test, "II", "III", 2, or 3.  
}
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
temp2<-mlma(y=sim.111$y, data1=data2)
summary(temp2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{summary statistics}
