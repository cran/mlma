\name{data.org}
\alias{data.org}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Transform and Organize Data for Mediation Analysis
}
\description{
To transform variables and generate data sets for mediation analysis.
}
\usage{
data.org(x, levelx=1, levely=1, m, l1 = NULL, l2 = NULL, c1 = NULL, 
         c1r = rep(1, length(c1)), c2 = NULL, c2r = rep(1, length(c2)), 
         f01y = NULL, f10y = NULL, f02ky = NULL, f20ky = NULL, f01km1 = NULL, 
         f01km2 = NULL, f10km = NULL, level, weight = rep(1, length(x)))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
the vector of the predictive variable.
}
  \item{levelx}{
the level of x (1 or 2), 1 by default.
}
  \item{levely}{
the level of y (1 or 2), 1 by default.
}
  \item{m}{
the matrix or vector of mediators.
}
  \item{l1}{
the column numbers of level 1 continuous mediators in m or the list of names of the level 1 continuous mediators.
}
  \item{l2}{
the column numbers of level 2 continuous mediators in m or the list of names of the level 2 continuous mediators.
}
  \item{c1}{
the column numbers of level 1 categorical mediators in m or the list of names of the level 1 categorical mediators.
}
  \item{c1r}{
the reference groups of categorical mediators specified by c1.
}
  \item{c2}{
the column numbers of level 2 categorical mediators in m or the list of names of the level 2 categorical mediators.
}
  \item{c2r}{
the reference groups of categorical mediators specified by c2.
}
  \item{f01y}{
the transformation function expressions on level 2 predictive variable (x.j) in explaining y (eg, c("x^2","log(x)")).
}
  \item{f10y}{
the transformation function expressions on level 1 predictive variable (xij-x.j) in explaining y.
}
  \item{f02ky}{
the transformation-function-expression list on level 2 mediators (m.jk) in explaining y (eg, list(2:3,c("log(x)","sqrt(x)"), "2*x")). The first item lists column numbers/variable names of the level 2 mediators in m, which needs to be transformed. By that order, each of the rest items of f01ky list the transformation functional expressions for each mediator. The mediators not specified in the list will not be transformed in any way.
}
  \item{f20ky}{
the transformation-function-expression list on level 1 mediators (mijk-m.jk) in explaining y.The first item lists column numbers/variable names of the level 1 mediators in m, which needs to be transformed. By that order, each of the rest items of f02ky list the transformation functional expressions for each mediator. The mediators not specified in the list will not be transformed in any way.
}
  \item{f01km1}{
the transformation-function-expression list on level 2 predictor (x.j) in explaining the level 1 mediators. The first item lists column numbers/variable names of the level 1 mediators in m, which should be explained by the transformed predictor(s). By that order, each of the rest items of f01km1 lists the transformation functional expressions for the (aggregated) level 2 predictor in explaining each mediator. The mediators not specified in the list will be explained by the original format of the level 2 predictor only.
}
  \item{f01km2}{
the transformation-function-expression list on level 2 predictor (x.j) in explaining the level 2 mediators. The first item lists column numbers/variable names of the level 2 mediators in m, which should be explained by the transformed level 2 predictor(s). By that order, each of the rest items of f01km2 lists the transformation functional expressions for the predictor in explaining each mediator. The mediators not specified in the list will be explained by the original format of the (aggregated) level 2 predictor only.
}
  \item{f10km}{
the transformation-function-expression list on level 1 predictor (xij-x.j) in explaining the level 1 mediators. The first item lists column numbers/variable names of the level 1 mediators in m, which should be explained by the transformed level 1 predictor(s). By that order, each of the rest items of f10km lists the transformation functional expressions for the level 1 predictor in explaining each mediator. The mediators not specified in the list will be explained by the original format of the level 1 predictor only.
}
  \item{level}{
a vector that record the group number for each observation.
}
  \item{weight}{
the weight of cases in groups.
}
}

\details{
The arguments starting with "f" are used to specify the transformation functions of the predictor or mediators in explaining y, or the transformation functions of the predictor in explaining the mediators. If the name of the argument includes a "k", the transformation is on the mediators.  If the names of the arguments end with "y", the transformation is to explain the outcome. Otherwise, the transformation is on x to predict mediators (the argument ends with "m1" or "m" (for level 1 mediator), or "m2" (for level 2 mediator)). The functions corresponds to the functions in the following multilevel additive models, reading as f+subscript+superscript.  For example, f01y specifies \eqn{{f}_{01}^Y}.
\deqn{Y_{ij} = u_{0j}^Y(X_{.j}, \mathbf{M}_{.j}, \mathbf{Z}_{.j})+u_{1j}^Y(X_{ij}-X_{.j})+\sum_{k=1}^Ku_{2jk}^Y(M_{ijk}-M_{.jk})+u_{3j}^Y(\mathbf{Z}_{ij}-\mathbf{Z}_{.j})+r_{ij}^Y}

where \deqn{u_{0j}^Y(X_{.j}, \mathbf{M}_{.j}, \mathbf{Z}_{.j}) = c_{00}^Y + {\boldsymbol{\beta}_{01}^Y}^T\mathbf{f}_{01}^Y(X_{.j}) + \sum_{k=1}^K{\boldsymbol{\beta}_{02k}^Y}^T\mathbf{f}_{02k}^Y(M_{.jk}) + {\boldsymbol{\beta}_{03}^Y}^T\mathbf{Z}_{.j} + r_{0j}^Y}
\deqn{u_{1j}^Y(X_{ij}-X_{.j}) = {\boldsymbol{\beta}_{10}^Y}^T\mathbf{f}_{10}^Y(X_{ij}-X_{.j}) + r_{1j}^Y}
\deqn{u_{2jk}^Y(M_{ijk}-M_{.jk}) = {\boldsymbol{\beta}_{20k}^Y}^T\mathbf{f}_{20k}^Y(M_{ijk}-M_{.jk}) + r_{2jk}^Y, \mbox{~~~~~~}k=1,\ldots,K}
\deqn{u_{3j}^Y(\mathbf{Z}_{ij}-\mathbf{Z}_{.j}) = {\boldsymbol{\beta}_{30}^Y}^T(\mathbf{Z}_{ij}-\mathbf{Z}_{.j}) + r_{3j}^Y}
For \eqn{k=1,\ldots,K,}
\deqn{M_{.jk} = c_{00k}^M+u_{1jk}^M(X_{ij}-X_{.j})+r_{ijk}^M}
\deqn{u_{0jk}^M(X_{.j}) = c_{00k}^M + {\boldsymbol{\beta}_{01k}^M}^T\mathbf{f}_{01k}^{M1}(X_{.j}) + r_{0jk}^M}
\deqn{u_{1jk}^M(X_{ij}-X_{.j}) = {\boldsymbol{\beta}_{10k}^M}^T\mathbf{f}_{10k}^M(X_{ij}-X_{.j}) + r_{1jk}^M}
If for some k, \eqn{M_k} is level 2 variable,
\deqn{M_{.jk} = c_{00k}^M + {\boldsymbol{\beta}_{01k}^M}^T\mathbf{f}_{01k}^{M2}(X_{.j}) + r_{0jk}^M} 

The transformation function can be any function that is differentiable by the function deriv(), or the ifelse function with those functions. The transformation function can also be the ns() and bs() functions for natural and b spline basis.
}
\value{
The function returns a list with transformed and organized data with the following items:
  \item{x1}{the level 1 and 2 transformed predictor variable matrix in explaining y (eg, f01y(x.j) & f10y(xij-x.j)).}
  \item{l1x}{the column numbers of level 1 predictors in x1.}
  \item{l2x}{the column numbers of level 2 predictors in x1.}
  \item{m1y}{the level 1 mediator matrix in explaining y (eg, f20ky(mijk-m.jk) & mijk or binarized mijk for categorical mediators).}
  \item{m1}{a list where the first item identify column numbers of level 1 mediators in m (eg l1 and c1). For every mediator identified by m1[[1]] and by that order, each of the rest item identify the column number(s) in m1y the (transformed) value(s) of the mediator in explaining y.}
  \item{m2y}{the level 2 mediator (original or aggregated) matrix in explaining y (eg, f02ky(m.jk) & m.jk).}
  \item{m2}{a list where the first item identify column numbers of level 2 mediators in m (eg l2 and c2). For every mediator identified by m2[[1]] and by that order, each of the rest item identify the column number(s) in m2y the (transformed) value(s) of the mediator in explaining y.}
  \item{m12}{a list where the first item identify column numbers of aggregated level 2 mediators in m (eg, aggregated f20ky). For every mediator identified by m12[[1]] and by that order, each of the rest item identify the column number(s) in m2y the aggregated value(s) of the mediator in explaining y.}
  \item{xm1}{the (transformed) level 1 and level 2 predictor(s) in explaining level 1 mediators.}
  \item{fm11}{a list where the first item identify column numbers of level 1 mediators in m. For every mediator identified by fm11[[1]] and by that order, each of the rest item identify the column number(s) in xm1 the (transformed) level 1 predictor(s) in explaining the mediator.}
  \item{fm12}{a list where the first item identify column numbers of level 1 mediators in m. For every mediator identified by fm12[[1]] and by that order, each of the rest item identify the column number(s) in xm1 the (transformed/aggregated) level 2 predictor(s) in explaining the mediator.}
  \item{m.2}{a matrix of level 2 mediators (one row for each group).}
  \item{xm2}{the (transformed/aggregated) level 2 predictor(s) in explaining level 2 mediators (one row for each group).}
  \item{fm22}{a list where the first item identify column numbers of level 2 mediators in m. For every mediator identified by fm22[[1]] and by that order, each of the rest item identify the column number(s) in xm2 the (transformed) level 2 predictor(s) in explaining the mediator.}
  \item{x1.der, m2y.der, m1y.der, xm2.der, xm1.der}{the derivative of x1, m2y, m1y, xm2, and xm1 respectively.}
  \item{parameter}{The list of arguments used.}
}
\author{
Qingzhao Yu (qyu@lsuhsc.edu), Bin Li (bli@lsu.edu).
}
\examples{
data(sim.211)
example1<-data.org(x=sim.211$x, levelx=2, m=sim.211$m, l1=2,l2=1, c1=3, c1r=1,
                   f01y=c("x","log(x^2)"), f02ky=list(1,c("x","x^2")), 
                   f20ky=list(2,c("x","x^3")), f01km1=list(2,"sqrt(x)+3"), 
                   f01km2=list(1,c("x^1.2","x^2.3")), level=sim.211$level)

data(sim.111)
example2<-data.org(sim.111$x, levelx=1, m=sim.111$m, l1=1:2, c1=3, c1r=1,
                   f01y=c("x","log(x^2)"), f10y=c("x^2","sqrt(x+6)"), 
                   f20ky=list(2,c("x","x^3")), f01km1=list(2,"sqrt(x)+3"), 
                   f10km=list(2,"log(x+2)"), level=sim.111$level)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Data Transformation }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
