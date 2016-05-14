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
temp<-mlma(y=sim.111$y, biny=FALSE, data1=data1, x=sim.111$x, levelx=1, m=sim.111$m, 
                l1=1:2,c1=3, c1r=1,level=sim.111$level,org.data=TRUE,
                f01y=c("x","log(x^2)"), f10y=c("x^2","sqrt(x+6)"), 
                f20ky=list(2,c("x","x^3")), f01km1=list(2,"sqrt(x)+3"), 
                f10km=list(2,"log(x+2)"))
summary(temp)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{summary statistics}
